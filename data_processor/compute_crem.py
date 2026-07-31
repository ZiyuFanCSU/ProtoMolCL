import os
import json
import argparse
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import copy
import lmdb
import json
import pickle
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
from rdkit import Chem
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')

from functools import partial
from crem.crem import mutate_mol2  # Generate perturbations using CReM
from multiprocessing import Pool

# Import your own compute function here.
# For example, if compute is defined in mutation.py:
# from mutation import compute
# If compute is defined in this file, no import is required.

from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem.Scaffolds.MurckoScaffold import (
    GetScaffoldForMol,
    MakeScaffoldGeneric,
    MurckoScaffoldSmiles,
)


def get_scaffold(smi, generic=False, return_mol=False):
    if generic:
        assert not return_mol

    if isinstance(smi, str):
        mol = Chem.MolFromSmiles(smi)
    else:
        mol = smi

    if generic:
        if return_mol:
            return MakeScaffoldGeneric(mol)
        else:
            return Chem.MolToSmiles(MakeScaffoldGeneric(mol))
    else:
        if return_mol:
            return GetScaffoldForMol(mol)
        else:
            return Chem.MolToSmiles(GetScaffoldForMol(mol))


GLOBAL_FRAGMENT_PATH = None


def init_worker(fragment_path):
    """
    Store fragment_path when each worker process is initialized
    to avoid passing it repeatedly.
    """
    global GLOBAL_FRAGMENT_PATH
    GLOBAL_FRAGMENT_PATH = fragment_path


def filter_valid_smiles(smiles_list, max_amount=None, canonical=True):
    """
    Filter valid SMILES strings.

    Parameters
    ----------
    smiles_list : list
        A list of candidate molecular SMILES strings.

    max_amount : int or None
        Maximum number of valid molecules to retain.
        If None, all valid molecules are retained.

    canonical : bool
        Whether to convert SMILES strings into RDKit canonical form.
        True is recommended.

    Returns
    -------
    valid_smiles : list
        A list of valid molecular SMILES strings.
    """

    valid_smiles = []
    seen = set()

    for smi in smiles_list:

        # Skip empty values and non-string objects.
        if smi is None or not isinstance(smi, str):
            continue

        # Attempt to convert the SMILES string into an RDKit Mol object.
        mol = Chem.MolFromSmiles(smi)

        # If conversion fails, the SMILES string is invalid.
        if mol is None:
            continue

        # Perform additional chemical validity checks,
        # such as valence and aromaticity validation.
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue

        # Optionally convert the molecule to canonical SMILES for deduplication.
        if canonical:
            smi_new = Chem.MolToSmiles(mol, canonical=True)
        else:
            smi_new = smi

        # Remove duplicates.
        if smi_new in seen:
            continue

        seen.add(smi_new)
        valid_smiles.append(smi_new)

        # Stop after reaching the requested maximum number.
        if max_amount is not None and len(valid_smiles) >= max_amount:
            break

    return valid_smiles


def compute(smi, db_name, max_amount=100, max_atom=200):

    # Determine mutation indices.
    mol = Chem.MolFromSmiles(smi)
    if len(mol.GetAtoms()) > max_atom:
        return {}

    scaff_smi = get_scaffold(mol)

    mapped_inds = mol.GetSubstructMatch(Chem.MolFromSmiles(scaff_smi))
    replace_inds = set(np.arange(len(mol.GetAtoms()))) - set(mapped_inds)

    for idx in copy.deepcopy(replace_inds):
        for n_atom in mol.GetAtomWithIdx(int(idx)).GetNeighbors():
            replace_inds.add(n_atom.GetIdx())

    replace_ids = [int(i) for i in replace_inds]

    # Mutate the molecule.
    random_mutations = set()

    if replace_inds:
        for radius in [2, 3]:
            random_mutations.update(
                set(
                    mutate_mol2(
                        mol,
                        db_name=db_name,
                        max_size=5,
                        radius=radius,
                        replace_ids=replace_ids,
                    )
                )
            )

        # mutations = filter_by_scaffold(
        #     scaff_smi,
        #     list(random_mutations),
        #     identical=True,
        #     max_amount=max_amount
        # )

        # The scaffold is no longer required to remain identical.
        # Only chemically valid generated SMILES strings are retained.
        mutations = filter_valid_smiles(
            list(random_mutations),
            max_amount=max_amount,
            canonical=True,
        )

    output = {
        'smi': smi,
        'mutate': mutations,
    }

    return output


def safe_compute(smi):
    """
    Call compute safely so that an error in one molecule
    does not terminate the entire program.
    """
    try:
        if pd.isna(smi):
            return {
                "smi": smi,
                "mutate": [],
                "success": False,
                "error": "NaN SMILES",
            }

        smi = str(smi).strip()

        if smi == "":
            return {
                "smi": smi,
                "mutate": [],
                "success": False,
                "error": "Empty SMILES",
            }

        result = compute(smi, GLOBAL_FRAGMENT_PATH)

        # Normally, result has the form:
        # {'smi': smi, 'mutate': mutations}
        mutate = result.get("mutate", [])

        return {
            "smi": result.get("smi", smi),
            "mutate": mutate,
            "success": True,
            "error": "",
        }

    except Exception as e:
        return {
            "smi": smi,
            "mutate": [],
            "success": False,
            "error": str(e),
        }


def process_chunk(smiles_list, num_workers):
    """
    Process one chunk of SMILES strings using multiple processes.
    """
    results = []

    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=init_worker,
        initargs=(GLOBAL_FRAGMENT_PATH,),
    ) as executor:

        futures = [executor.submit(safe_compute, smi) for smi in smiles_list]

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Processing chunk",
        ):
            results.append(future.result())

    return results


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input_csv",
        type=str,
        default="./ProtoMolCL/data_process2/zinc_200w.csv",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        default="./ProtoMolCL/data_process3/zinc_crem.csv",
    )
    parser.add_argument(
        "--fragment_path",
        type=str,
        default="./ProtoMolCL/data/dataset.db",
    )

    parser.add_argument("--smiles_col", type=str, default="smiles")
    parser.add_argument("--chunksize", type=int, default=10000)
    parser.add_argument("--num_workers", type=int, default=16)

    args = parser.parse_args()

    global GLOBAL_FRAGMENT_PATH
    GLOBAL_FRAGMENT_PATH = args.fragment_path

    # Remove the existing output file to avoid duplicate appends.
    if os.path.exists(args.output_csv):
        os.remove(args.output_csv)

    first_write = True
    total_processed = 0

    reader = pd.read_csv(args.input_csv, chunksize=args.chunksize)

    for chunk_id, chunk in enumerate(reader):
        print(f"\nProcessing chunk {chunk_id}, size = {len(chunk)}")

        if args.smiles_col not in chunk.columns:
            raise ValueError(f"Column not found in CSV: {args.smiles_col}")

        smiles_list = chunk[args.smiles_col].tolist()

        results = process_chunk(
            smiles_list=smiles_list,
            num_workers=args.num_workers,
        )

        out_rows = []

        for res in results:
            out_rows.append({
                "smi": res["smi"],

                # mutate may be a list, so save it as a JSON string
                # for convenient downstream loading.
                "mutate": json.dumps(res["mutate"], ensure_ascii=False),

                "success": res["success"],
                "error": res["error"],
            })

        out_df = pd.DataFrame(out_rows)

        out_df.to_csv(
            args.output_csv,
            mode="w" if first_write else "a",
            index=False,
            header=first_write,
        )

        first_write = False
        total_processed += len(chunk)

        print(f"Finished chunk {chunk_id}, total processed = {total_processed}")

    print("All done!")
    print(f"Result saved to: {args.output_csv}")


if __name__ == "__main__":
    compute(
        "CCCCCN",
        "./ProtoMolCL/data/dataset.db",
    )
    main()
