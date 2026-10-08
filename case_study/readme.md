# ProtoMolCL Case Studies

This directory contains four groups of analyses accompanying our KDD rebuttal, illustrating complementary information in low-, mid-, and high-granularity molecular representations.

## 1. MoleculeNet BACE cases

`case_study_BACE_1.png`–`case_study_BACE_3.png` use molecules from the MoleculeNet BACE dataset. Regions receiving high attention scores closely correspond to molecular regions implicated by protein-pocket simulation, providing context for interpreting the attention patterns.

## 2. Activity-cliff molecular pairs

`case_study_cliffpairs_1.png`–`case_study_cliffpairs_6.png` compare structurally similar molecules with large activity differences. Selected cases show complementary attention preferences: low and mid highlight local modifications and connecting fragments, while high emphasizes additional shared-scaffold regions.

## 3. MoleculeNet molecular cases

`case_study_moleculenet_1.png`–`case_study_moleculenet_6.png` show differences in attended regions across low, mid, and high granularities. **Our linear probes also show differences in the information carried by these representations: low favors functional-group information, mid BRICS-fragment information, and high global molecular descriptors. Even if all three branches attend to the same atom, their encoded information can differ because that atom has a different representation at each granularity.**

## 4. Semantic organization and cluster interpretation

Our selected cases are presented in `tsne_semantic.png`, which compares representations colored by functional-group, BRICS-fragment, and global-property clusters, and in `tsne_low_representatives.png`, `tsne_mid_representatives.png`, and `tsne_high_representatives.png`, which show representative molecules.

`cluster_interpretation_report1.md`–`cluster_interpretation_report3.md` contain cluster interpretations for the original training dataset, providing background for the chemical semantics used in training.

## Reading the attention visualizations

The attention analysis contains two components:

- **GPS attention:** atom-to-atom attention within the GraphGPS encoder, describing how atom representations exchange information across the molecular graph.
- **One-prompt attention:** attention from a single learned prompt to atom representations during aggregation, describing how atoms are weighted when forming a molecular representation.

GPS attention describes information exchange; one-prompt attention describes readout weighting. Highlighted atoms can encode broader molecular context, so attention locations alone do not establish semantic granularity or causal explanations of activity.


