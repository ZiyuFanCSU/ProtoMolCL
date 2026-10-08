# Functional-group prototype cluster interpretation
Total valid molecules: 2000000
Number of clusters: 34

## How to read this report
- `Top RDKit fragments` is usually the most intuitive chemical interpretation.
- `cluster_presence` means the fraction of molecules in this cluster containing that fragment.
- `presence_enrichment` = cluster presence - global presence. Larger means more characteristic of this cluster.
- Representative SMILES are the molecules closest to the cluster center in the RDKit2DNormalized feature space.

## Cluster 0
- Count: 28721
- Ratio: 0.0144

### Top RDKit fragments
- fr_aryl_methyl: presence=0.940, global=0.398, enrich=0.542
- fr_furan: presence=0.395, global=0.056, enrich=0.339
- fr_thiophene: presence=0.154, global=0.074, enrich=0.080
- fr_sulfone: presence=0.089, global=0.025, enrich=0.064
- fr_amide: presence=0.754, global=0.701, enrich=0.054
- fr_C_O_noCOO: presence=0.818, global=0.779, enrich=0.039
- fr_Ndealkylation1: presence=0.099, global=0.061, enrich=0.038
- fr_piperzine: presence=0.086, global=0.062, enrich=0.025
- fr_C_O: presence=0.818, global=0.797, enrich=0.021
- fr_phenol_noOrthoHbond: presence=0.044, global=0.024, enrich=0.020
- fr_phenol: presence=0.044, global=0.025, enrich=0.019
- fr_C_S: presence=0.048, global=0.031, enrich=0.017
- fr_hdrzone: presence=0.032, global=0.020, enrich=0.012
- fr_Ar_OH: presence=0.046, global=0.034, enrich=0.012
- fr_benzene: presence=0.844, global=0.834, enrich=0.010
- fr_unbrch_alkane: presence=0.034, global=0.025, enrich=0.008
- fr_NH2: presence=0.064, global=0.057, enrich=0.008
- fr_priamide: presence=0.026, global=0.019, enrich=0.007
- fr_oxime: presence=0.011, global=0.005, enrich=0.006
- fr_aldehyde: presence=0.009, global=0.003, enrich=0.006

### Top high RDKit2DNormalized features
- ('fr_furan', <class 'numpy.float64'>): z_diff=1.488
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.106
- ('fr_sulfone', <class 'numpy.float64'>): z_diff=0.417
- ('fr_thiophene', <class 'numpy.float64'>): z_diff=0.303
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.248
- ('fr_Ndealkylation1', <class 'numpy.float64'>): z_diff=0.171
- ('fr_phenol_noOrthoHbond', <class 'numpy.float64'>): z_diff=0.130
- ('fr_phenol', <class 'numpy.float64'>): z_diff=0.127

### Representative SMILES
- `CCCN(C(=O)c1ccc(CC)cc1)[C@H]1CCNC1`
- `CCCN(C(=O)c1ccc(CC)cc1)[C@@H]1CCNC1`
- `CCN(Cc1ccccc1C)C(=O)[C@H]1CCNC1`
- `CCc1ccc(C(=O)N[C@H]2CC[C@@H](N(C)C)C2)cc1`
- `Cc1cccc(CN(C(=O)[C@@H]2CCN2)C2CCCCC2)c1`

## Cluster 1
- Count: 51320
- Ratio: 0.0257

### Top RDKit fragments
- fr_bicyclic: presence=0.947, global=0.397, enrich=0.550
- fr_ether: presence=0.998, global=0.514, enrich=0.484
- fr_methoxy: presence=0.709, global=0.230, enrich=0.479
- fr_ester: presence=0.524, global=0.107, enrich=0.416
- fr_Ar_N: presence=0.928, global=0.531, enrich=0.397
- fr_thiazole: presence=0.343, global=0.062, enrich=0.281
- fr_allylic_oxid: presence=0.247, global=0.055, enrich=0.192
- fr_para_hydroxylation: presence=0.401, global=0.216, enrich=0.186
- fr_halogen: presence=0.534, global=0.384, enrich=0.150
- fr_benzene: presence=0.964, global=0.834, enrich=0.129
- fr_NH0: presence=0.978, global=0.857, enrich=0.121
- fr_Ar_OH: presence=0.080, global=0.034, enrich=0.045
- fr_phenol_noOrthoHbond: presence=0.056, global=0.024, enrich=0.032
- fr_phenol: presence=0.057, global=0.025, enrich=0.032
- fr_furan: presence=0.080, global=0.056, enrich=0.024
- fr_imidazole: presence=0.069, global=0.054, enrich=0.015
- fr_hdrzone: presence=0.030, global=0.020, enrich=0.010
- fr_oxazole: presence=0.022, global=0.014, enrich=0.008
- fr_lactone: presence=0.011, global=0.006, enrich=0.005
- fr_term_acetylene: presence=0.008, global=0.003, enrich=0.005

### Top high RDKit2DNormalized features
- ('fr_ester', <class 'numpy.float64'>): z_diff=1.351
- ('fr_thiazole', <class 'numpy.float64'>): z_diff=1.180
- ('fr_methoxy', <class 'numpy.float64'>): z_diff=1.158
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=1.132
- ('fr_ether', <class 'numpy.float64'>): z_diff=1.043
- ('fr_allylic_oxid', <class 'numpy.float64'>): z_diff=0.866
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.761
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.676

### Representative SMILES
- `COC(=O)c1ccccc1Cn1cnc(=O)c2cc(OC)c(F)cc21`
- `COC(=O)c1cc2c(ccn2-c2ccc(F)cc2)n1Cc1cccc(OC)c1`
- `COC(=O)c1cccc(Cn2cnc3cc(Cl)ccc3c2=O)c1`
- `COC(=O)c1cccc(Cn2cnc3cc(F)ccc3c2=O)c1`
- `COC(=O)Cc1nn(Cc2ccccc2Cl)c(=O)c2ccccc12`

## Cluster 2
- Count: 89626
- Ratio: 0.0448

### Top RDKit fragments
- fr_bicyclic: presence=1.000, global=0.397, enrich=0.603
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_Ar_N: presence=0.999, global=0.531, enrich=0.468
- fr_amide: presence=0.997, global=0.701, enrich=0.296
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_pyridine: presence=0.259, global=0.138, enrich=0.121
- fr_NH1: presence=0.741, global=0.642, enrich=0.100
- fr_imidazole: presence=0.138, global=0.054, enrich=0.084
- fr_para_hydroxylation: presence=0.293, global=0.216, enrich=0.077
- fr_aniline: presence=0.531, global=0.456, enrich=0.075
- fr_thiophene: presence=0.147, global=0.074, enrich=0.073
- fr_thiazole: presence=0.112, global=0.062, enrich=0.050
- fr_sulfide: presence=0.155, global=0.122, enrich=0.033
- fr_oxazole: presence=0.036, global=0.014, enrich=0.021
- fr_piperzine: presence=0.075, global=0.062, enrich=0.013
- fr_furan: presence=0.068, global=0.056, enrich=0.012
- fr_priamide: presence=0.021, global=0.019, enrich=0.002
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=1.231
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.969
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.774
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.586
- ('fr_imidazole', <class 'numpy.float64'>): z_diff=0.386
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.345
- ('fr_thiophene', <class 'numpy.float64'>): z_diff=0.277

### Representative SMILES
- `Cc1nc([C@H]2CCNC2)nc2c1CCC(=O)N2Cc1ccccc1`
- `C=Cc1ccc(C(=O)N[C@H]2CCCc3nc(N(C)C)ncc32)cc1`
- `Cc1cc(C)nc(NCCCC(=O)N2CCc3ccccc3C2)n1`
- `Cc1cc(C)c(C(=O)N[C@H]2CCCc3nc(N(C)C)ncc32)cc1C`
- `O=C(Cn1ncc(N2CCCC2)cc1=O)N[C@@H]1CCc2ccccc2C1`

## Cluster 3
- Count: 71580
- Ratio: 0.0358

### Top RDKit fragments
- fr_bicyclic: presence=1.000, global=0.397, enrich=0.603
- fr_Ar_N: presence=0.983, global=0.531, enrich=0.452
- fr_aryl_methyl: presence=0.649, global=0.398, enrich=0.251
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_pyridine: presence=0.276, global=0.138, enrich=0.138
- fr_imidazole: presence=0.175, global=0.054, enrich=0.121
- fr_para_hydroxylation: presence=0.294, global=0.216, enrich=0.078
- fr_thiophene: presence=0.134, global=0.074, enrich=0.060
- fr_Ar_OH: presence=0.068, global=0.034, enrich=0.034
- fr_thiazole: presence=0.089, global=0.062, enrich=0.028
- fr_oxazole: presence=0.038, global=0.014, enrich=0.024
- fr_nitrile: presence=0.047, global=0.030, enrich=0.017
- fr_sulfide: presence=0.135, global=0.122, enrich=0.012
- fr_phenol_noOrthoHbond: presence=0.036, global=0.024, enrich=0.012
- fr_phenol: presence=0.036, global=0.025, enrich=0.012
- fr_tetrazole: presence=0.018, global=0.012, enrich=0.007
- fr_unbrch_alkane: presence=0.032, global=0.025, enrich=0.006
- fr_hdrzone: presence=0.026, global=0.020, enrich=0.006
- fr_Imine: presence=0.035, global=0.031, enrich=0.005
- fr_SH: presence=0.005, global=0.003, enrich=0.003

### Top high RDKit2DNormalized features
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=1.246
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.970
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.912
- ('fr_imidazole', <class 'numpy.float64'>): z_diff=0.550
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=0.512
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.395
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.232
- ('fr_thiophene', <class 'numpy.float64'>): z_diff=0.227

### Representative SMILES
- `CCc1ccc2oc(=O)cc(CN(C)[C@@H](C)c3ccncn3)c2c1`
- `CCCN(Cc1cnc(C)cn1)C1Cc2ccccc2C1`
- `CCc1nnc(CN2CCc3ccccc3C2(C)C)o1`
- `CCCN(Cc1cc(=O)n(C)c(=O)n1C)[C@@H]1CCc2ccccc2C1`
- `CN(Cc1cc(=O)n(C)c(=O)n1C)[C@@H]1CCCc2ccccc21`

## Cluster 4
- Count: 54557
- Ratio: 0.0273

### Top RDKit fragments
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_Ar_N: presence=0.998, global=0.531, enrich=0.467
- fr_amide: presence=0.995, global=0.701, enrich=0.294
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_thiazole: presence=0.131, global=0.062, enrich=0.069
- fr_pyridine: presence=0.194, global=0.138, enrich=0.056
- fr_Ndealkylation2: presence=0.197, global=0.161, enrich=0.037
- fr_imidazole: presence=0.088, global=0.054, enrich=0.033
- fr_morpholine: presence=0.069, global=0.042, enrich=0.027
- fr_urea: presence=0.092, global=0.069, enrich=0.023
- fr_Ndealkylation1: presence=0.082, global=0.061, enrich=0.020
- fr_oxazole: presence=0.033, global=0.014, enrich=0.019
- fr_piperzine: presence=0.076, global=0.062, enrich=0.015
- fr_tetrazole: presence=0.025, global=0.012, enrich=0.013
- fr_priamide: presence=0.030, global=0.019, enrich=0.011
- fr_SH: presence=0.003, global=0.003, enrich=0.001
- fr_furan: presence=0.056, global=0.056, enrich=0.000
- fr_sulfone: presence=0.025, global=0.025, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.933
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.657
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.579
- ('fr_thiazole', <class 'numpy.float64'>): z_diff=0.304
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.226
- ('fr_oxazole', <class 'numpy.float64'>): z_diff=0.174
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.159

### Representative SMILES
- `CCn1nccc1CCN1CCNC(=O)C[C@H]1c1ccccc1`
- `Cc1ccnc([C@H](NC(=O)c2ccc(CN(C)C)cc2)C2CC2)n1`
- `Cc1cccc(-c2ncc(CN3CCC(=O)NC[C@@H]3C)cn2)c1`
- `CCc1ncc(C(=O)N[C@@H](CCN(C)C)c2ccccc2)cn1`
- `Cc1ccn(CC(=O)N[C@H](CCN(C)C)c2ccccc2)n1`

## Cluster 5
- Count: 70236
- Ratio: 0.0351

### Top RDKit fragments
- fr_halogen: presence=1.000, global=0.384, enrich=0.616
- fr_ether: presence=1.000, global=0.514, enrich=0.486
- fr_amide: presence=0.972, global=0.701, enrich=0.271
- fr_methoxy: presence=0.454, global=0.230, enrich=0.224
- fr_C_O_noCOO: presence=0.999, global=0.779, enrich=0.220
- fr_C_O: presence=0.999, global=0.797, enrich=0.202
- fr_NH1: presence=0.801, global=0.642, enrich=0.160
- fr_benzene: presence=0.983, global=0.834, enrich=0.149
- fr_ester: presence=0.249, global=0.107, enrich=0.142
- fr_urea: presence=0.179, global=0.069, enrich=0.110
- fr_alkyl_halide: presence=0.161, global=0.054, enrich=0.107
- fr_aniline: presence=0.549, global=0.456, enrich=0.093
- fr_imide: presence=0.141, global=0.050, enrich=0.092
- fr_barbitur: presence=0.042, global=0.006, enrich=0.037
- fr_morpholine: presence=0.070, global=0.042, enrich=0.028
- fr_hdrzone: presence=0.042, global=0.020, enrich=0.022
- fr_C_S: presence=0.050, global=0.031, enrich=0.019
- fr_Ndealkylation1: presence=0.073, global=0.061, enrich=0.012
- fr_lactone: presence=0.015, global=0.006, enrich=0.009
- fr_priamide: presence=0.027, global=0.019, enrich=0.008

### Top high RDKit2DNormalized features
- ('fr_halogen', <class 'numpy.float64'>): z_diff=1.272
- ('fr_ether', <class 'numpy.float64'>): z_diff=0.984
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.719
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.683
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.664
- ('fr_methoxy', <class 'numpy.float64'>): z_diff=0.548
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.507
- ('fr_alkyl_halide', <class 'numpy.float64'>): z_diff=0.494

### Representative SMILES
- `CN(C)C(=O)c1cccc(NC(=O)COc2cccc(F)c2)c1`
- `CCOc1ccc(CNC(=O)[C@H]2CC(=O)N(c3ccc(Br)cc3)C2)cc1`
- `O=C(COc1ccc(N2CCCC2=O)cc1)N[C@H]1C[C@H]1c1cccc(F)c1`
- `C[C@@H](Oc1ccc(Cl)cc1)C(=O)NCc1ccc(N2CCCC2=O)cc1`
- `CC(=O)Nc1ccc(S(=O)(=O)Oc2ccc(Cl)cc2CN(C[C@H]2CCCO2)C(=O)C(C)C)cc1`

## Cluster 6
- Count: 57256
- Ratio: 0.0286

### Top RDKit fragments
- fr_halogen: presence=0.715, global=0.384, enrich=0.331
- fr_amide: presence=0.996, global=0.701, enrich=0.295
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_NH1: presence=0.816, global=0.642, enrich=0.174
- fr_urea: presence=0.222, global=0.069, enrich=0.153
- fr_imide: presence=0.153, global=0.050, enrich=0.104
- fr_piperzine: presence=0.144, global=0.062, enrich=0.082
- fr_alkyl_halide: presence=0.125, global=0.054, enrich=0.071
- fr_Ndealkylation1: presence=0.125, global=0.061, enrich=0.064
- fr_benzene: presence=0.880, global=0.834, enrich=0.046
- fr_sulfone: presence=0.070, global=0.025, enrich=0.046
- fr_priamide: presence=0.060, global=0.019, enrich=0.041
- fr_thiophene: presence=0.103, global=0.074, enrich=0.029
- fr_NH2: presence=0.085, global=0.057, enrich=0.028
- fr_furan: presence=0.079, global=0.056, enrich=0.023
- fr_nitrile: presence=0.045, global=0.030, enrich=0.015
- fr_C_S: presence=0.042, global=0.031, enrich=0.011
- fr_hdrzone: presence=0.031, global=0.020, enrich=0.011
- fr_unbrch_alkane: presence=0.035, global=0.025, enrich=0.009

### Top high RDKit2DNormalized features
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.780
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.703
- ('fr_halogen', <class 'numpy.float64'>): z_diff=0.689
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.656
- ('fr_urea', <class 'numpy.float64'>): z_diff=0.619
- ('fr_imide', <class 'numpy.float64'>): z_diff=0.480
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.435
- ('fr_piperzine', <class 'numpy.float64'>): z_diff=0.345

### Representative SMILES
- `CC(C)[C@H](NC(=O)c1ccc(Cl)cc1)C(=O)N1CC[S@@](=O)C(C)(C)C1`
- `CN(C)C(=O)CCNC(=O)c1ccc(Br)cc1`
- `CN(C)C(=O)CCNC(=O)[C@H]1C[C@H]1c1cccc(Cl)c1`
- `CC(C)CNC(=O)[C@H](C)N(Cc1ccc(F)cc1)C(=O)C(C)C`
- `CC[C@@H](C(=O)N1CCCNC(=O)[C@@H]1C)c1ccc(F)cc1`

## Cluster 7
- Count: 72730
- Ratio: 0.0364

### Top RDKit fragments
- fr_sulfonamd: presence=1.000, global=0.132, enrich=0.868
- fr_aniline: presence=0.811, global=0.456, enrich=0.355
- fr_halogen: presence=0.672, global=0.384, enrich=0.288
- fr_NH1: presence=0.912, global=0.642, enrich=0.271
- fr_amide: presence=0.931, global=0.701, enrich=0.230
- fr_C_O_noCOO: presence=0.950, global=0.779, enrich=0.171
- fr_benzene: presence=0.989, global=0.834, enrich=0.155
- fr_C_O: presence=0.950, global=0.797, enrich=0.153
- fr_methoxy: presence=0.305, global=0.230, enrich=0.075
- fr_ether: presence=0.570, global=0.514, enrich=0.056
- fr_Ndealkylation1: presence=0.095, global=0.061, enrich=0.034
- fr_morpholine: presence=0.058, global=0.042, enrich=0.016
- fr_piperzine: presence=0.073, global=0.062, enrich=0.011
- fr_alkyl_halide: presence=0.064, global=0.054, enrich=0.010
- fr_unbrch_alkane: presence=0.032, global=0.025, enrich=0.007
- fr_para_hydroxylation: presence=0.216, global=0.216, enrich=0.000
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000
- fr_isocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_thiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isothiocyan: presence=0.000, global=0.000, enrich=-0.000

### Top high RDKit2DNormalized features
- ('fr_sulfonamd', <class 'numpy.float64'>): z_diff=2.581
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.781
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.709
- ('fr_halogen', <class 'numpy.float64'>): z_diff=0.597
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.531
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.490
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.304
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.225

### Representative SMILES
- `CCOc1cccc(CNC(=O)CN(c2ccc(F)cc2)S(C)(=O)=O)c1`
- `CCN(CC)S(=O)(=O)c1ccc(OCC(=O)Nc2ccc(Br)cc2)cc1`
- `CN(C)S(=O)(=O)c1cccc(NC(=O)COc2ccc(Br)cc2)c1`
- `CCOc1ccc(N(CC(=O)NC23CC4CC(CC(C4)C2)C3)S(=O)(=O)c2ccc(Br)cc2)cc1`
- `CCOc1ccc([C@@H](C)NC(=O)CN(c2ccc(F)cc2)S(C)(=O)=O)cc1`

## Cluster 8
- Count: 60537
- Ratio: 0.0303

### Top RDKit fragments
- fr_methoxy: presence=0.849, global=0.230, enrich=0.619
- fr_ether: presence=1.000, global=0.514, enrich=0.486
- fr_aryl_methyl: presence=0.852, global=0.398, enrich=0.454
- fr_Ar_N: presence=0.949, global=0.531, enrich=0.417
- fr_bicyclic: presence=0.775, global=0.397, enrich=0.378
- fr_aniline: presence=0.826, global=0.456, enrich=0.370
- fr_amide: presence=0.928, global=0.701, enrich=0.227
- fr_NH1: presence=0.851, global=0.642, enrich=0.209
- fr_C_O_noCOO: presence=0.954, global=0.779, enrich=0.175
- fr_para_hydroxylation: presence=0.380, global=0.216, enrich=0.164
- fr_C_O: presence=0.954, global=0.797, enrich=0.157
- fr_NH0: presence=0.984, global=0.857, enrich=0.127
- fr_benzene: presence=0.949, global=0.834, enrich=0.115
- fr_pyridine: presence=0.226, global=0.138, enrich=0.088
- fr_thiazole: presence=0.123, global=0.062, enrich=0.061
- fr_ester: presence=0.151, global=0.107, enrich=0.043
- fr_sulfide: presence=0.158, global=0.122, enrich=0.036
- fr_imidazole: presence=0.080, global=0.054, enrich=0.026
- fr_thiophene: presence=0.094, global=0.074, enrich=0.021
- fr_oxazole: presence=0.026, global=0.014, enrich=0.011

### Top high RDKit2DNormalized features
- ('fr_methoxy', <class 'numpy.float64'>): z_diff=1.495
- ('fr_ether', <class 'numpy.float64'>): z_diff=0.982
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=0.927
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.840
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=0.772
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.736
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.483
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.448

### Representative SMILES
- `COc1cccc([C@@H]2c3ccccc3CCN2CC(=O)Nc2cc(C)on2)c1`
- `COc1cc2c(cc1OC)[C@H](c1ccccc1)N(CC(=O)Nc1cc(C)on1)CC2`
- `COc1ccc(C(C)(C)C)cc1NC(=O)c1nn(C)c(=O)c2ccccc21`
- `COc1ccc(-n2nc3c(c2NC(=O)c2ccccc2C)C[S@@](=O)C3)cc1`
- `COc1ccc(C(=O)Nc2c3nsnc3ccc2C)cc1`

## Cluster 9
- Count: 78423
- Ratio: 0.0392

### Top RDKit fragments
- fr_Ar_N: presence=1.000, global=0.531, enrich=0.469
- fr_pyridine: presence=0.469, global=0.138, enrich=0.331
- fr_amide: presence=0.998, global=0.701, enrich=0.297
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_NH1: presence=0.709, global=0.642, enrich=0.067
- fr_urea: presence=0.119, global=0.069, enrich=0.050
- fr_piperzine: presence=0.106, global=0.062, enrich=0.044
- fr_thiazole: presence=0.104, global=0.062, enrich=0.042
- fr_morpholine: presence=0.073, global=0.042, enrich=0.031
- fr_Ndealkylation1: presence=0.091, global=0.061, enrich=0.029
- fr_thiophene: presence=0.102, global=0.074, enrich=0.029
- fr_tetrazole: presence=0.031, global=0.012, enrich=0.020
- fr_priamide: presence=0.033, global=0.019, enrich=0.014
- fr_imidazole: presence=0.065, global=0.054, enrich=0.011
- fr_Ndealkylation2: presence=0.169, global=0.161, enrich=0.008
- fr_furan: presence=0.063, global=0.056, enrich=0.007
- fr_nitrile: presence=0.034, global=0.030, enrich=0.004
- fr_oxazole: presence=0.018, global=0.014, enrich=0.003

### Top high RDKit2DNormalized features
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.952
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.923
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.743
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.625
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.227
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.212
- ('fr_urea', <class 'numpy.float64'>): z_diff=0.212
- ('fr_tetrazole', <class 'numpy.float64'>): z_diff=0.198

### Representative SMILES
- `C[C@@H](Cc1cccnc1)NC(=O)c1ccc(-c2ccc(CN(C)C)cn2)cc1`
- `O=C(NCc1cccnc1)c1cnn(Cc2ccccc2)c1`
- `O=C(Cn1cnc(-c2ccccc2)cc1=O)NCc1ccccn1`
- `O=C(NCc1cccnc1)c1nnc(-c2ccccc2)o1`
- `CC(C)c1c(C(=O)NCC2(c3ccccc3)CCC2)cnn1-c1ccccn1`

## Cluster 10
- Count: 50611
- Ratio: 0.0253

### Top RDKit fragments
- fr_bicyclic: presence=1.000, global=0.397, enrich=0.603
- fr_ether: presence=0.873, global=0.514, enrich=0.359
- fr_ester: presence=0.350, global=0.107, enrich=0.242
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_amide: presence=0.877, global=0.701, enrich=0.176
- fr_methoxy: presence=0.359, global=0.230, enrich=0.129
- fr_imide: presence=0.173, global=0.050, enrich=0.124
- fr_thiophene: presence=0.153, global=0.074, enrich=0.079
- fr_allylic_oxid: presence=0.121, global=0.055, enrich=0.066
- fr_benzene: presence=0.865, global=0.834, enrich=0.030
- fr_urea: presence=0.095, global=0.069, enrich=0.026
- fr_lactone: presence=0.025, global=0.006, enrich=0.019
- fr_Ndealkylation1: presence=0.078, global=0.061, enrich=0.017
- fr_morpholine: presence=0.054, global=0.042, enrich=0.012
- fr_priamide: presence=0.031, global=0.019, enrich=0.012
- fr_alkyl_carbamate: presence=0.013, global=0.004, enrich=0.009
- fr_amidine: presence=0.026, global=0.019, enrich=0.007
- fr_unbrch_alkane: presence=0.031, global=0.025, enrich=0.006
- fr_sulfone: presence=0.030, global=0.025, enrich=0.005

### Top high RDKit2DNormalized features
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=1.242
- ('fr_ester', <class 'numpy.float64'>): z_diff=0.790
- ('fr_ether', <class 'numpy.float64'>): z_diff=0.771
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.727
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.706
- ('fr_imide', <class 'numpy.float64'>): z_diff=0.571
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.451
- ('fr_methoxy', <class 'numpy.float64'>): z_diff=0.321

### Representative SMILES
- `CC(C)NC(=O)[C@H]1c2ccccc2C(=O)N(C[C@H]2CCCO2)C12CCCC2`
- `CC[C@@H](C)NC(=O)[C@H]1c2ccccc2C(=O)N(C[C@H]2CCCO2)C12CCCC2`
- `CCOc1cccc(CNC(=O)CC2=C3CCCC[C@@H]3N(C(C)C)C2=O)c1`
- `CC(C)NC(=O)[C@@H]1c2ccccc2C(=O)N(C[C@@H]2CCCO2)C12CCCC2`
- `C=C1c2ccccc2C(=O)N1CC(=O)NCCO[C@@H]1CCCC[C@H]1C`

## Cluster 11
- Count: 76874
- Ratio: 0.0384

### Top RDKit fragments
- fr_piperdine: presence=1.000, global=0.136, enrich=0.864
- fr_Ndealkylation2: presence=0.974, global=0.161, enrich=0.813
- fr_Ar_N: presence=0.993, global=0.531, enrich=0.462
- fr_aryl_methyl: presence=0.630, global=0.398, enrich=0.232
- fr_pyridine: presence=0.321, global=0.138, enrich=0.183
- fr_amide: presence=0.857, global=0.701, enrich=0.156
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_C_O_noCOO: presence=0.868, global=0.779, enrich=0.090
- fr_imidazole: presence=0.129, global=0.054, enrich=0.075
- fr_C_O: presence=0.869, global=0.797, enrich=0.072
- fr_thiazole: presence=0.089, global=0.062, enrich=0.027
- fr_Nhpyrrole: presence=0.096, global=0.078, enrich=0.018
- fr_Ar_NH: presence=0.096, global=0.078, enrich=0.018
- fr_oxazole: presence=0.028, global=0.014, enrich=0.013
- fr_priamide: presence=0.025, global=0.019, enrich=0.006
- fr_tetrazole: presence=0.016, global=0.012, enrich=0.004
- fr_HOCCN: presence=0.013, global=0.011, enrich=0.001
- fr_lactam: presence=0.002, global=0.001, enrich=0.001
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000
- fr_isocyan: presence=0.000, global=0.000, enrich=-0.000

### Top high RDKit2DNormalized features
- ('fr_piperdine', <class 'numpy.float64'>): z_diff=2.539
- ('fr_Ndealkylation2', <class 'numpy.float64'>): z_diff=2.230
- ('fr_NH0', <class 'numpy.float64'>): z_diff=1.040
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.940
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.525
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=0.473
- ('fr_imidazole', <class 'numpy.float64'>): z_diff=0.345
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.320

### Representative SMILES
- `CC(=O)N1CCC[C@@H](c2ccnc(-c3cccc(C)c3)n2)C1`
- `Cc1cnc(C(=O)N2CCC(Cc3ccccc3)CC2)cn1`
- `Cc1ccc(-c2noc([C@H]3CCCCN3C(=O)C34CC5CC(CC(C5)C3)C4)n2)cc1`
- `Cc1ccccc1-c1noc([C@@H]2CCC(=O)N(C3CC3)C2)n1`
- `CCc1noc([C@@H]2CCC(=O)N(CCc3ccccc3)C2)n1`

## Cluster 12
- Count: 56413
- Ratio: 0.0282

### Top RDKit fragments
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_aniline: presence=0.997, global=0.456, enrich=0.540
- fr_Ar_N: presence=0.989, global=0.531, enrich=0.458
- fr_amide: presence=0.997, global=0.701, enrich=0.296
- fr_NH1: presence=0.908, global=0.642, enrich=0.266
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_thiazole: presence=0.141, global=0.062, enrich=0.079
- fr_pyridine: presence=0.212, global=0.138, enrich=0.074
- fr_sulfide: presence=0.193, global=0.122, enrich=0.071
- fr_urea: presence=0.129, global=0.069, enrich=0.060
- fr_piperzine: presence=0.093, global=0.062, enrich=0.032
- fr_tetrazole: presence=0.029, global=0.012, enrich=0.017
- fr_oxazole: presence=0.021, global=0.014, enrich=0.007
- fr_guanido: presence=0.005, global=0.003, enrich=0.003
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000
- fr_thiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isothiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isocyan: presence=0.000, global=0.000, enrich=-0.000

### Top high RDKit2DNormalized features
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_aniline', <class 'numpy.float64'>): z_diff=1.078
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.928
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.630
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.628
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.547
- ('fr_thiazole', <class 'numpy.float64'>): z_diff=0.344
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.264

### Representative SMILES
- `CCn1nccc1CCN[C@H](C)c1cccc(N2CCCC2=O)c1`
- `CCc1nc(CN[C@@H](C)c2ccc(N3CCCC3=O)cc2)no1`
- `Cc1ncncc1C(=O)NCCc1ccc(N(C)C)cc1`
- `CNC(=O)[C@H]1CCN(c2cc(C)nc(-c3ccccc3)n2)C1`
- `Cc1ccc(C(=O)NCc2nccc(N3CCCC3)n2)cc1`

## Cluster 13
- Count: 61116
- Ratio: 0.0306

### Top RDKit fragments
- fr_halogen: presence=1.000, global=0.384, enrich=0.616
- fr_Ar_N: presence=1.000, global=0.531, enrich=0.469
- fr_amide: presence=0.988, global=0.701, enrich=0.287
- fr_pyridine: presence=0.372, global=0.138, enrich=0.234
- fr_C_O_noCOO: presence=0.993, global=0.779, enrich=0.215
- fr_C_O: presence=0.993, global=0.797, enrich=0.196
- fr_alkyl_halide: presence=0.211, global=0.054, enrich=0.157
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_aniline: presence=0.590, global=0.456, enrich=0.134
- fr_NH1: presence=0.771, global=0.642, enrich=0.129
- fr_benzene: presence=0.889, global=0.834, enrich=0.054
- fr_thiazole: presence=0.107, global=0.062, enrich=0.045
- fr_piperzine: presence=0.096, global=0.062, enrich=0.034
- fr_oxazole: presence=0.034, global=0.014, enrich=0.020
- fr_tetrazole: presence=0.025, global=0.012, enrich=0.013
- fr_imidazole: presence=0.061, global=0.054, enrich=0.007
- fr_urea: presence=0.076, global=0.069, enrich=0.007
- fr_priamide: presence=0.023, global=0.019, enrich=0.004
- fr_SH: presence=0.004, global=0.003, enrich=0.001
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_halogen', <class 'numpy.float64'>): z_diff=1.270
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.910
- ('fr_alkyl_halide', <class 'numpy.float64'>): z_diff=0.717
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.672
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.569
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.518
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.260
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.218

### Representative SMILES
- `O=C(Cc1ccccc1Cl)N[C@@H]1CCN(c2ccc(=O)n(CC3CC3)n2)C1`
- `CN(C)c1ccnc(CNC(=O)c2ccc(F)cc2)n1`
- `CNC(=O)CCCN(C)c1nnc(-c2ccccc2Cl)s1`
- `CC[C@@H](C)NC(=O)c1cnc(-c2ccc(F)cc2)nc1N(C)C`
- `C[C@H](CNC(=O)c1ccc(N(C)C)cc1F)Cn1cccn1`

## Cluster 14
- Count: 38331
- Ratio: 0.0192

### Top RDKit fragments
- fr_piperdine: presence=0.977, global=0.136, enrich=0.841
- fr_Ndealkylation2: presence=0.965, global=0.161, enrich=0.805
- fr_halogen: presence=1.000, global=0.384, enrich=0.616
- fr_amide: presence=0.901, global=0.701, enrich=0.200
- fr_C_O_noCOO: presence=0.917, global=0.779, enrich=0.138
- fr_NH0: presence=0.995, global=0.857, enrich=0.138
- fr_alkyl_halide: presence=0.180, global=0.054, enrich=0.126
- fr_C_O: presence=0.917, global=0.797, enrich=0.120
- fr_benzene: presence=0.906, global=0.834, enrich=0.071
- fr_urea: presence=0.119, global=0.069, enrich=0.050
- fr_priamide: presence=0.037, global=0.019, enrich=0.019
- fr_NH2: presence=0.068, global=0.057, enrich=0.012
- fr_HOCCN: presence=0.018, global=0.011, enrich=0.006
- fr_lactam: presence=0.002, global=0.001, enrich=0.001
- fr_Al_OH: presence=0.077, global=0.076, enrich=0.001
- fr_oxazole: presence=0.015, global=0.014, enrich=0.001
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000
- fr_thiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isothiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isocyan: presence=0.000, global=0.000, enrich=-0.000

### Top high RDKit2DNormalized features
- ('fr_piperdine', <class 'numpy.float64'>): z_diff=2.472
- ('fr_Ndealkylation2', <class 'numpy.float64'>): z_diff=2.205
- ('fr_halogen', <class 'numpy.float64'>): z_diff=1.257
- ('fr_alkyl_halide', <class 'numpy.float64'>): z_diff=0.576
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.466
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.289
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.250
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.225

### Representative SMILES
- `CN(C)S(=O)(=O)NC[C@@H]1CCCN(C(=O)/C=C/c2cccc(Cl)c2)C1`
- `CN(C)CC(C)(C)CNC(=O)C1CCN(Cc2ccc(Cl)cc2)CC1`
- `CCCN1CC2(C[C@@H]1C(=O)NCc1cccc(F)c1)CCN(C)CC2`
- `C[C@@H]1CCCCN1C1CN(CC(=O)NCCc2ccc(F)cc2)C1`
- `O=C(CCN1CCCCC1)NC[C@H](c1ccc(Cl)cc1)N1CCCC1`

## Cluster 15
- Count: 58507
- Ratio: 0.0293

### Top RDKit fragments
- fr_Al_COO: presence=1.000, global=0.031, enrich=0.969
- fr_COO: presence=1.000, global=0.046, enrich=0.954
- fr_COO2: presence=1.000, global=0.046, enrich=0.954
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_allylic_oxid: presence=0.097, global=0.055, enrich=0.041
- fr_Ndealkylation2: presence=0.200, global=0.161, enrich=0.039
- fr_alkyl_carbamate: presence=0.025, global=0.004, enrich=0.021
- fr_unbrch_alkane: presence=0.040, global=0.025, enrich=0.015
- fr_lactam: presence=0.007, global=0.001, enrich=0.006
- fr_sulfide: presence=0.126, global=0.122, enrich=0.004
- fr_oxime: presence=0.008, global=0.005, enrich=0.003
- fr_guanido: presence=0.004, global=0.003, enrich=0.001
- fr_phenol_noOrthoHbond: presence=0.025, global=0.024, enrich=0.000
- fr_nitroso: presence=0.000, global=0.000, enrich=0.000
- fr_phenol: presence=0.025, global=0.025, enrich=0.000
- fr_epoxide: presence=0.001, global=0.001, enrich=0.000
- fr_phos_ester: presence=0.000, global=0.000, enrich=0.000
- fr_phos_acid: presence=0.000, global=0.000, enrich=0.000
- fr_azide: presence=0.000, global=0.000, enrich=0.000
- fr_term_acetylene: presence=0.003, global=0.003, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_Al_COO', <class 'numpy.float64'>): z_diff=5.587
- ('fr_COO', <class 'numpy.float64'>): z_diff=4.578
- ('fr_COO2', <class 'numpy.float64'>): z_diff=4.577
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.835
- ('fr_alkyl_carbamate', <class 'numpy.float64'>): z_diff=0.325
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.240
- ('fr_allylic_oxid', <class 'numpy.float64'>): z_diff=0.198
- ('fr_lactam', <class 'numpy.float64'>): z_diff=0.175

### Representative SMILES
- `O=C(O)CN(Cc1ccccc1)C(=O)[C@@H]1C[C@@H]1C1CC1`
- `CC(=O)N(CC(=O)O)Cc1ccccc1`
- `CC(C)CC(=O)N(CC(=O)O)Cc1ccccc1`
- `C[C@H](C(=O)O)N(Cc1ccccc1)C(=O)C1CCC1`
- `C[C@@H](C1CC1)N(Cc1ccccc1)C(=O)CCC(=O)O`

## Cluster 16
- Count: 59372
- Ratio: 0.0297

### Top RDKit fragments
- fr_piperdine: presence=0.916, global=0.136, enrich=0.780
- fr_Ndealkylation2: presence=0.939, global=0.161, enrich=0.778
- fr_amide: presence=0.912, global=0.701, enrich=0.211
- fr_C_O_noCOO: presence=0.933, global=0.779, enrich=0.154
- fr_C_O: presence=0.933, global=0.797, enrich=0.136
- fr_NH0: presence=0.992, global=0.857, enrich=0.135
- fr_urea: presence=0.156, global=0.069, enrich=0.087
- fr_ether: presence=0.567, global=0.514, enrich=0.053
- fr_NH2: presence=0.094, global=0.057, enrich=0.037
- fr_priamide: presence=0.053, global=0.019, enrich=0.034
- fr_Ndealkylation1: presence=0.095, global=0.061, enrich=0.034
- fr_methoxy: presence=0.245, global=0.230, enrich=0.014
- fr_alkyl_carbamate: presence=0.010, global=0.004, enrich=0.005
- fr_furan: presence=0.061, global=0.056, enrich=0.005
- fr_thiophene: presence=0.078, global=0.074, enrich=0.005
- fr_lactam: presence=0.004, global=0.001, enrich=0.003
- fr_morpholine: presence=0.044, global=0.042, enrich=0.002
- fr_quatN: presence=0.003, global=0.001, enrich=0.002
- fr_term_acetylene: presence=0.005, global=0.003, enrich=0.002
- fr_epoxide: presence=0.002, global=0.001, enrich=0.001

### Top high RDKit2DNormalized features
- ('fr_piperdine', <class 'numpy.float64'>): z_diff=2.292
- ('fr_Ndealkylation2', <class 'numpy.float64'>): z_diff=2.134
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.538
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.476
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.420
- ('fr_urea', <class 'numpy.float64'>): z_diff=0.360
- ('fr_priamide', <class 'numpy.float64'>): z_diff=0.253
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.241

### Representative SMILES
- `CN1CCC(N(CCc2ccccc2)CC(=O)NC[C@H]2CCCO2)CC1`
- `O=C(CN1CCCC1)NC[C@@H]1CCC2(CCN(Cc3ccccc3)CC2)O1`
- `CN1CCC(OCCCNC(=O)C2(c3ccccc3)CCN(C)CC2)CC1`
- `CCOc1ccc(CN2CCC(N3CCC[C@H](C(=O)NC4CC4)C3)CC2)cc1`
- `C[C@H]1CCCN(CCCNC(=O)C2(c3ccccc3)CCOCC2)C1`

## Cluster 17
- Count: 62564
- Ratio: 0.0313

### Top RDKit fragments
- fr_Ar_NH: presence=1.000, global=0.078, enrich=0.922
- fr_Nhpyrrole: presence=1.000, global=0.078, enrich=0.922
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_Ar_N: presence=1.000, global=0.531, enrich=0.469
- fr_NH1: presence=1.000, global=0.642, enrich=0.358
- fr_bicyclic: presence=0.656, global=0.397, enrich=0.259
- fr_imidazole: presence=0.205, global=0.054, enrich=0.151
- fr_pyridine: presence=0.225, global=0.138, enrich=0.087
- fr_NH0: presence=0.927, global=0.857, enrich=0.071
- fr_para_hydroxylation: presence=0.262, global=0.216, enrich=0.046
- fr_tetrazole: presence=0.041, global=0.012, enrich=0.029
- fr_thiophene: presence=0.081, global=0.074, enrich=0.007
- fr_Ar_OH: presence=0.036, global=0.034, enrich=0.002
- fr_guanido: presence=0.003, global=0.003, enrich=0.000
- fr_sulfide: presence=0.122, global=0.122, enrich=0.000
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000
- fr_isothiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_thiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_diazo: presence=0.000, global=0.000, enrich=-0.000

### Top high RDKit2DNormalized features
- ('fr_Ar_NH', <class 'numpy.float64'>): z_diff=3.414
- ('fr_Nhpyrrole', <class 'numpy.float64'>): z_diff=3.414
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.997
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.932
- ('fr_imidazole', <class 'numpy.float64'>): z_diff=0.684
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=0.530
- ('fr_tetrazole', <class 'numpy.float64'>): z_diff=0.287

### Representative SMILES
- `Cc1cc(=O)[nH]c(CCNC(=O)[C@@H](C)N2CCc3ccccc3C2)n1`
- `Cc1n[nH]c(C2CC2)c1C(=O)N[C@H](C)CN1CCc2ccccc2C1`
- `CCCc1cc(C(=O)N[C@H](C)CN2CCc3ccccc3C2)[nH]n1`
- `Cc1n[nH]c(C2CC2)c1C(=O)N[C@@H](C)CN1CCc2ccccc2C1`
- `Cc1cc(=O)[nH]c(CCNC(=O)[C@H](C)N2CCc3ccccc3C2)n1`

## Cluster 18
- Count: 57135
- Ratio: 0.0286

### Top RDKit fragments
- fr_para_hydroxylation: presence=0.943, global=0.216, enrich=0.727
- fr_bicyclic: presence=0.881, global=0.397, enrich=0.485
- fr_aniline: presence=0.864, global=0.456, enrich=0.408
- fr_amide: presence=0.977, global=0.701, enrich=0.277
- fr_C_O_noCOO: presence=0.985, global=0.779, enrich=0.206
- fr_C_O: presence=0.985, global=0.797, enrich=0.188
- fr_benzene: presence=1.000, global=0.834, enrich=0.165
- fr_NH1: presence=0.776, global=0.642, enrich=0.135
- fr_imide: presence=0.126, global=0.050, enrich=0.076
- fr_thiazole: presence=0.106, global=0.062, enrich=0.044
- fr_urea: presence=0.110, global=0.069, enrich=0.041
- fr_furan: presence=0.093, global=0.056, enrich=0.037
- fr_NH0: presence=0.892, global=0.857, enrich=0.035
- fr_Ndealkylation1: presence=0.088, global=0.061, enrich=0.026
- fr_sulfide: presence=0.145, global=0.122, enrich=0.023
- fr_N_O: presence=0.018, global=0.003, enrich=0.015
- fr_piperzine: presence=0.073, global=0.062, enrich=0.012
- fr_priamide: presence=0.030, global=0.019, enrich=0.012
- fr_oxazole: presence=0.026, global=0.014, enrich=0.012
- fr_hdrzone: presence=0.026, global=0.020, enrich=0.006

### Top high RDKit2DNormalized features
- ('fr_para_hydroxylation', <class 'numpy.float64'>): z_diff=1.752
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=1.001
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.824
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.685
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.650
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.628
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.586
- ('fr_imide', <class 'numpy.float64'>): z_diff=0.353

### Representative SMILES
- `CCN(CC)Cc1ccc(CNC(=O)[C@@H]2CN(C(C)=O)c3ccccc3O2)cc1`
- `CC(C)N(CCCNC(=O)CCC(=O)N1C[C@@H](C)Oc2ccccc21)Cc1ccccc1`
- `CN(CCCNC(=O)CCC(=O)N1CCOc2ccccc21)Cc1ccccc1`
- `O=C(CN1CCOc2ccccc2C1)N[C@@H]1CC(=O)N(c2ccccc2)C1`
- `CN(C)c1cccc(CNC(=O)CCN2C(=O)COc3ccccc32)c1`

## Cluster 19
- Count: 62076
- Ratio: 0.0310

### Top RDKit fragments
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_amide: presence=0.955, global=0.701, enrich=0.254
- fr_NH1: presence=0.893, global=0.642, enrich=0.252
- fr_aniline: presence=0.697, global=0.456, enrich=0.241
- fr_ether: presence=0.755, global=0.514, enrich=0.240
- fr_C_O_noCOO: presence=0.974, global=0.779, enrich=0.196
- fr_C_O: presence=0.974, global=0.797, enrich=0.177
- fr_benzene: presence=0.994, global=0.834, enrich=0.159
- fr_halogen: presence=0.535, global=0.384, enrich=0.151
- fr_imide: presence=0.161, global=0.050, enrich=0.111
- fr_urea: presence=0.165, global=0.069, enrich=0.096
- fr_methoxy: presence=0.309, global=0.230, enrich=0.079
- fr_ester: presence=0.176, global=0.107, enrich=0.068
- fr_barbitur: presence=0.053, global=0.006, enrich=0.048
- fr_C_S: presence=0.073, global=0.031, enrich=0.043
- fr_hdrzine: presence=0.011, global=0.005, enrich=0.005
- fr_lactone: presence=0.010, global=0.006, enrich=0.005
- fr_lactam: presence=0.003, global=0.001, enrich=0.002
- fr_unbrch_alkane: presence=0.027, global=0.025, enrich=0.002
- fr_guanido: presence=0.003, global=0.003, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.756
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.731
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.674
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.671
- ('fr_barbitur', <class 'numpy.float64'>): z_diff=0.605
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.525
- ('fr_imide', <class 'numpy.float64'>): z_diff=0.515

### Representative SMILES
- `Cc1ccc(C(C)C)cc1OCC(=O)Nc1ccc(F)cc1C(=O)N(C)C`
- `Cc1ccc(OCC(=O)NC[C@@H]2CC(=O)N(c3ccc(F)cc3)C2)cc1`
- `Cc1ccc(F)cc1OCC(=O)Nc1ccc(C(=O)N(C)C)cc1`
- `Cc1cccc(OCCNC(=O)[C@H]2CC(=O)N(c3ccc(Cl)cc3)C2)c1`
- `CCOc1ccc(N2C[C@H](C(=O)NCc3ccc(F)c(C)c3)CC2=O)cc1`

## Cluster 20
- Count: 92329
- Ratio: 0.0462

### Top RDKit fragments
- fr_Al_OH_noTert: presence=0.996, global=0.064, enrich=0.932
- fr_Al_OH: presence=1.000, global=0.076, enrich=0.924
- fr_HOCCN: presence=0.177, global=0.011, enrich=0.166
- fr_Ndealkylation2: presence=0.209, global=0.161, enrich=0.048
- fr_Imine: presence=0.070, global=0.031, enrich=0.039
- fr_piperzine: presence=0.084, global=0.062, enrich=0.022
- fr_Ndealkylation1: presence=0.079, global=0.061, enrich=0.017
- fr_piperdine: presence=0.150, global=0.136, enrich=0.014
- fr_lactone: presence=0.009, global=0.006, enrich=0.003
- fr_quatN: presence=0.003, global=0.001, enrich=0.003
- fr_phos_acid: presence=0.002, global=0.000, enrich=0.002
- fr_phos_ester: presence=0.002, global=0.000, enrich=0.002
- fr_phenol: presence=0.026, global=0.025, enrich=0.001
- fr_phenol_noOrthoHbond: presence=0.026, global=0.024, enrich=0.001
- fr_unbrch_alkane: presence=0.027, global=0.025, enrich=0.001
- fr_epoxide: presence=0.002, global=0.001, enrich=0.001
- fr_term_acetylene: presence=0.005, global=0.003, enrich=0.001
- fr_sulfone: presence=0.026, global=0.025, enrich=0.001
- fr_oxime: presence=0.006, global=0.005, enrich=0.001
- fr_azide: presence=0.001, global=0.000, enrich=0.001

### Top high RDKit2DNormalized features
- ('fr_Al_OH_noTert', <class 'numpy.float64'>): z_diff=3.762
- ('fr_Al_OH', <class 'numpy.float64'>): z_diff=3.487
- ('fr_HOCCN', <class 'numpy.float64'>): z_diff=1.582
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.268
- ('fr_Imine', <class 'numpy.float64'>): z_diff=0.242
- ('fr_Ndealkylation2', <class 'numpy.float64'>): z_diff=0.136
- ('fr_phos_acid', <class 'numpy.float64'>): z_diff=0.110
- ('fr_phos_ester', <class 'numpy.float64'>): z_diff=0.110

### Representative SMILES
- `CN(Cc1ccccc1)S(=O)(=O)N(CCO)CC12CC3CC(CC(C3)C1)C2`
- `CN([C@@H](CCO)c1ccccc1)S(=O)(=O)N1CCCCCC1`
- `CC(C)N(CC#CC[C@](O)(c1ccccc1)C1CCCC1)CCO`
- `C=C(C)CN(CCO)Cc1ccccc1`
- `C[C@@H](CO)CN(Cc1ccccc1)CC1CC1`

## Cluster 21
- Count: 34160
- Ratio: 0.0171

### Top RDKit fragments
- fr_bicyclic: presence=1.000, global=0.397, enrich=0.603
- fr_aryl_methyl: presence=0.958, global=0.398, enrich=0.560
- fr_sulfonamd: presence=0.544, global=0.132, enrich=0.412
- fr_aniline: presence=0.861, global=0.456, enrich=0.404
- fr_NH1: presence=0.846, global=0.642, enrich=0.204
- fr_amide: presence=0.866, global=0.701, enrich=0.165
- fr_benzene: presence=0.990, global=0.834, enrich=0.155
- fr_C_O_noCOO: presence=0.892, global=0.779, enrich=0.113
- fr_C_O: presence=0.892, global=0.797, enrich=0.095
- fr_thiophene: presence=0.120, global=0.074, enrich=0.046
- fr_para_hydroxylation: presence=0.256, global=0.216, enrich=0.040
- fr_imide: presence=0.071, global=0.050, enrich=0.021
- fr_amidine: presence=0.022, global=0.019, enrich=0.003
- fr_sulfone: presence=0.027, global=0.025, enrich=0.002
- fr_priamide: presence=0.021, global=0.019, enrich=0.002
- fr_N_O: presence=0.005, global=0.003, enrich=0.001
- fr_Imine: presence=0.032, global=0.031, enrich=0.001
- fr_benzodiazepine: presence=0.001, global=0.000, enrich=0.001
- fr_lactam: presence=0.001, global=0.001, enrich=0.000
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=1.228
- ('fr_sulfonamd', <class 'numpy.float64'>): z_diff=1.225
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.141
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.819
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.726
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.407
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.371
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.245

### Representative SMILES
- `CCNC(=O)c1ccc2c(c1)CCN2S(=O)(=O)c1ccc(C)cc1`
- `CCC[C@@H](C)NC(=O)c1ccc2c(c1)CCN2S(=O)(=O)c1ccc(C)cc1`
- `Cc1ccccc1CNC(=O)c1ccc2c(c1)CCN2S(C)(=O)=O`
- `CC(C)(C)NS(=O)(=O)c1ccc2c(c1)CCCN2C(=O)c1ccccc1`
- `CCS(=O)(=O)N1c2ccc(C(=O)N[C@@H](C)CCc3ccccc3)cc2C[C@H]1C`

## Cluster 22
- Count: 51826
- Ratio: 0.0259

### Top RDKit fragments
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_Ar_N: presence=0.973, global=0.531, enrich=0.441
- fr_NH0: presence=0.999, global=0.857, enrich=0.142
- fr_pyridine: presence=0.248, global=0.138, enrich=0.110
- fr_sulfonamd: presence=0.205, global=0.132, enrich=0.073
- fr_C_S: presence=0.092, global=0.031, enrich=0.061
- fr_thiazole: presence=0.092, global=0.062, enrich=0.030
- fr_nitrile: presence=0.057, global=0.030, enrich=0.027
- fr_imidazole: presence=0.081, global=0.054, enrich=0.026
- fr_tetrazole: presence=0.033, global=0.012, enrich=0.021
- fr_oxazole: presence=0.033, global=0.014, enrich=0.019
- fr_morpholine: presence=0.058, global=0.042, enrich=0.016
- fr_sulfone: presence=0.033, global=0.025, enrich=0.008
- fr_piperzine: presence=0.069, global=0.062, enrich=0.007
- fr_SH: presence=0.009, global=0.003, enrich=0.007
- fr_alkyl_halide: presence=0.057, global=0.054, enrich=0.003
- fr_aldehyde: presence=0.006, global=0.003, enrich=0.003
- fr_Ar_OH: presence=0.037, global=0.034, enrich=0.002
- fr_guanido: presence=0.004, global=0.003, enrich=0.001
- fr_azo: presence=0.002, global=0.001, enrich=0.001

### Top high RDKit2DNormalized features
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.914
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.745
- ('fr_C_S', <class 'numpy.float64'>): z_diff=0.359
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.316
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.241
- ('fr_tetrazole', <class 'numpy.float64'>): z_diff=0.216
- ('fr_sulfonamd', <class 'numpy.float64'>): z_diff=0.216

### Representative SMILES
- `CCN(Cc1nc(C2CC2)no1)Cc1ccccc1C`
- `Cc1ccc(-c2ccc([C@@H](C)N(C)C)cc2)nn1`
- `C=Cc1cccc(CN(C)Cc2cnn(C)c2C)c1`
- `CC[C@@H]1C=CCN1Cc1cn(C)nc1-c1ccccc1`
- `CCc1ccc(CN(Cc2nnc(C(C)C)o2)C2CC2)cc1`

## Cluster 23
- Count: 29039
- Ratio: 0.0145

### Top RDKit fragments
- fr_Ar_COO: presence=1.000, global=0.015, enrich=0.985
- fr_COO: presence=1.000, global=0.046, enrich=0.954
- fr_COO2: presence=1.000, global=0.046, enrich=0.954
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_furan: presence=0.108, global=0.056, enrich=0.052
- fr_benzene: presence=0.882, global=0.834, enrich=0.047
- fr_C_S: presence=0.069, global=0.031, enrich=0.038
- fr_Imine: presence=0.063, global=0.031, enrich=0.032
- fr_pyridine: presence=0.166, global=0.138, enrich=0.027
- fr_phenol: presence=0.042, global=0.025, enrich=0.017
- fr_amidine: presence=0.036, global=0.019, enrich=0.017
- fr_Ar_OH: presence=0.050, global=0.034, enrich=0.016
- fr_aniline: presence=0.466, global=0.456, enrich=0.010
- fr_hdrzone: presence=0.029, global=0.020, enrich=0.009
- fr_imide: presence=0.058, global=0.050, enrich=0.009
- fr_barbitur: presence=0.010, global=0.006, enrich=0.005
- fr_azo: presence=0.005, global=0.001, enrich=0.003
- fr_guanido: presence=0.004, global=0.003, enrich=0.001
- fr_aldehyde: presence=0.003, global=0.003, enrich=0.000
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_Ar_COO', <class 'numpy.float64'>): z_diff=8.044
- ('fr_COO', <class 'numpy.float64'>): z_diff=4.578
- ('fr_COO2', <class 'numpy.float64'>): z_diff=4.577
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.525
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.257
- ('fr_furan', <class 'numpy.float64'>): z_diff=0.226
- ('fr_C_S', <class 'numpy.float64'>): z_diff=0.223
- ('fr_Imine', <class 'numpy.float64'>): z_diff=0.202

### Representative SMILES
- `Cc1c(CNC[C@H](C)c2ccccc2)c(C(=O)O)c(C)n1Cc1ccccc1`
- `Cc1c(CNCc2ccccc2)c(C(=O)O)c(C)n1Cc1ccccc1`
- `Cc1c(CNC2CCCCC2)c(C(=O)O)c(C)n1Cc1ccccc1`
- `Cc1c(CN[C@H](C)C23CC4CC(CC(C4)C2)C3)c(C(=O)O)c(C)n1Cc1ccccc1`
- `Cc1c(CN[C@@H]2CCC[C@H](C)[C@@H]2C)c(C(=O)O)c(C)n1Cc1ccccc1`

## Cluster 24
- Count: 77714
- Ratio: 0.0389

### Top RDKit fragments
- fr_ketone: presence=1.000, global=0.062, enrich=0.938
- fr_ketone_Topliss: presence=0.990, global=0.054, enrich=0.936
- fr_allylic_oxid: presence=0.297, global=0.055, enrich=0.241
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_ester: presence=0.302, global=0.107, enrich=0.194
- fr_ether: presence=0.705, global=0.514, enrich=0.191
- fr_Al_OH: presence=0.206, global=0.076, enrich=0.130
- fr_Al_OH_noTert: presence=0.173, global=0.064, enrich=0.109
- fr_bicyclic: presence=0.497, global=0.397, enrich=0.101
- fr_methoxy: presence=0.329, global=0.230, enrich=0.099
- fr_benzene: presence=0.895, global=0.834, enrich=0.061
- fr_dihydropyridine: presence=0.056, global=0.003, enrich=0.053
- fr_phenol_noOrthoHbond: presence=0.061, global=0.024, enrich=0.037
- fr_phenol: presence=0.061, global=0.025, enrich=0.037
- fr_Ar_OH: presence=0.065, global=0.034, enrich=0.030
- fr_imide: presence=0.064, global=0.050, enrich=0.015
- fr_unbrch_alkane: presence=0.034, global=0.025, enrich=0.008
- fr_lactone: presence=0.011, global=0.006, enrich=0.005
- fr_hdrzine: presence=0.011, global=0.005, enrich=0.005

### Top high RDKit2DNormalized features
- ('fr_ketone_Topliss', <class 'numpy.float64'>): z_diff=4.176
- ('fr_ketone', <class 'numpy.float64'>): z_diff=3.912
- ('fr_allylic_oxid', <class 'numpy.float64'>): z_diff=1.085
- ('fr_dihydropyridine', <class 'numpy.float64'>): z_diff=1.066
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.950
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.943
- ('fr_ester', <class 'numpy.float64'>): z_diff=0.635
- ('fr_Al_OH', <class 'numpy.float64'>): z_diff=0.491

### Representative SMILES
- `CC(C)N(C(=O)c1ccccc1-c1ccc2c(c1)OCCC2=O)C(C)C`
- `CCOc1cccc(CN(C(=O)c2cccc(C(C)=O)c2)C2CC2)c1`
- `CC(=O)c1ccccc1-c1cccc(C(=O)N(C)CCC[C@@H]2CCCO2)c1`
- `CC(=O)c1ccc2c(c1)CN(C(=O)OC(C)(C)C)CC2`
- `CCN(Cc1ccc2c(c1)OCCO2)C(=O)CCCOc1ccc(C(C)=O)cc1`

## Cluster 25
- Count: 35349
- Ratio: 0.0177

### Top RDKit fragments
- fr_ArN: presence=1.000, global=0.020, enrich=0.979
- fr_NH2: presence=1.000, global=0.057, enrich=0.943
- fr_aniline: presence=0.970, global=0.456, enrich=0.514
- fr_Ar_N: presence=0.848, global=0.531, enrich=0.317
- fr_nitrile: presence=0.082, global=0.030, enrich=0.052
- fr_NH0: presence=0.906, global=0.857, enrich=0.049
- fr_sulfide: presence=0.169, global=0.122, enrich=0.047
- fr_pyridine: presence=0.179, global=0.138, enrich=0.041
- fr_Nhpyrrole: presence=0.109, global=0.078, enrich=0.031
- fr_Ar_NH: presence=0.109, global=0.078, enrich=0.031
- fr_thiophene: presence=0.097, global=0.074, enrich=0.023
- fr_priamide: presence=0.033, global=0.019, enrich=0.014
- fr_hdrzine: presence=0.016, global=0.005, enrich=0.011
- fr_unbrch_alkane: presence=0.034, global=0.025, enrich=0.009
- fr_azo: presence=0.009, global=0.001, enrich=0.008
- fr_thiazole: presence=0.069, global=0.062, enrich=0.007
- fr_ketone: presence=0.065, global=0.062, enrich=0.004
- fr_phos_ester: presence=0.004, global=0.000, enrich=0.003
- fr_phos_acid: presence=0.004, global=0.000, enrich=0.003
- fr_nitroso: presence=0.001, global=0.000, enrich=0.001

### Top high RDKit2DNormalized features
- ('fr_ArN', <class 'numpy.float64'>): z_diff=6.971
- ('fr_NH2', <class 'numpy.float64'>): z_diff=4.073
- ('fr_aniline', <class 'numpy.float64'>): z_diff=1.052
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.663
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.365
- ('fr_nitrile', <class 'numpy.float64'>): z_diff=0.310
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.280
- ('fr_azo', <class 'numpy.float64'>): z_diff=0.230

### Representative SMILES
- `Nc1nc(-c2ccccc2)nn1C(=O)C1CC1`
- `CC[C@@H](C(=O)n1cnnc1N)c1ccccc1`
- `C[C@H](NCc1cnc(N)nc1)c1ccccc1`
- `CC1(C)C[C@@]1(CNc1nccc(N)n1)c1ccccc1`
- `CC1(C)C[C@]1(CNc1nccc(N)n1)c1ccccc1`

## Cluster 26
- Count: 29761
- Ratio: 0.0149

### Top RDKit fragments
- fr_nitro_arom: presence=1.000, global=0.042, enrich=0.958
- fr_nitro: presence=1.000, global=0.048, enrich=0.952
- fr_nitro_arom_nonortho: presence=0.656, global=0.026, enrich=0.630
- fr_benzene: presence=1.000, global=0.834, enrich=0.166
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_hdrzone: presence=0.098, global=0.020, enrich=0.078
- fr_ester: presence=0.181, global=0.107, enrich=0.073
- fr_allylic_oxid: presence=0.120, global=0.055, enrich=0.065
- fr_Imine: presence=0.085, global=0.031, enrich=0.055
- fr_Ar_OH: presence=0.077, global=0.034, enrich=0.043
- fr_phenol: presence=0.064, global=0.025, enrich=0.039
- fr_phenol_noOrthoHbond: presence=0.063, global=0.024, enrich=0.039
- fr_thiazole: presence=0.096, global=0.062, enrich=0.034
- fr_sulfonamd: presence=0.154, global=0.132, enrich=0.022
- fr_bicyclic: presence=0.415, global=0.397, enrich=0.018
- fr_lactone: presence=0.019, global=0.006, enrich=0.014
- fr_oxime: presence=0.016, global=0.005, enrich=0.011
- fr_ether: presence=0.525, global=0.514, enrich=0.011
- fr_nitrile: presence=0.041, global=0.030, enrich=0.011
- fr_ketone: presence=0.072, global=0.062, enrich=0.010

### Top high RDKit2DNormalized features
- ('fr_nitro_arom', <class 'numpy.float64'>): z_diff=4.829
- ('fr_nitro', <class 'numpy.float64'>): z_diff=4.516
- ('fr_nitro_arom_nonortho', <class 'numpy.float64'>): z_diff=4.009
- ('fr_hdrzone', <class 'numpy.float64'>): z_diff=0.577
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.568
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.553
- ('fr_Imine', <class 'numpy.float64'>): z_diff=0.331
- ('fr_allylic_oxid', <class 'numpy.float64'>): z_diff=0.303

### Representative SMILES
- `O=[N+]([O-])c1cccc(Oc2cccc([N+](=O)[O-])c2)c1`
- `CCCOc1ccc(/C=C/c2ccnc(-c3cccc([N+](=O)[O-])c3)[n+]2[O-])cc1`
- `CC(C)(C)c1ccc(OCc2nc(-c3cccc([N+](=O)[O-])c3)no2)cc1`
- `CC(C)c1ccc(OCc2nc(-c3ccc([N+](=O)[O-])cc3)no2)cc1`
- `CCCOc1ccc(/C=C\c2ccnc(-c3cccc([N+](=O)[O-])c3)[n+]2[O-])cc1`

## Cluster 27
- Count: 66851
- Ratio: 0.0334

### Top RDKit fragments
- fr_sulfide: presence=0.999, global=0.122, enrich=0.877
- fr_halogen: presence=0.727, global=0.384, enrich=0.343
- fr_amide: presence=0.977, global=0.701, enrich=0.276
- fr_amidine: presence=0.279, global=0.019, enrich=0.260
- fr_Imine: presence=0.262, global=0.031, enrich=0.231
- fr_C_O_noCOO: presence=0.986, global=0.779, enrich=0.207
- fr_C_O: presence=0.986, global=0.797, enrich=0.189
- fr_benzene: presence=0.968, global=0.834, enrich=0.134
- fr_imide: presence=0.135, global=0.050, enrich=0.086
- fr_NH0: presence=0.910, global=0.857, enrich=0.053
- fr_C_S: presence=0.084, global=0.031, enrich=0.053
- fr_aniline: presence=0.507, global=0.456, enrich=0.051
- fr_para_hydroxylation: presence=0.266, global=0.216, enrich=0.050
- fr_NH1: presence=0.671, global=0.642, enrich=0.030
- fr_sulfone: presence=0.052, global=0.025, enrich=0.027
- fr_hdrzone: presence=0.038, global=0.020, enrich=0.018
- fr_alkyl_halide: presence=0.069, global=0.054, enrich=0.016
- fr_methoxy: presence=0.242, global=0.230, enrich=0.012
- fr_hdrzine: presence=0.014, global=0.005, enrich=0.009
- fr_phenol_noOrthoHbond: presence=0.032, global=0.024, enrich=0.007

### Top high RDKit2DNormalized features
- ('fr_sulfide', <class 'numpy.float64'>): z_diff=2.672
- ('fr_amidine', <class 'numpy.float64'>): z_diff=1.907
- ('fr_Imine', <class 'numpy.float64'>): z_diff=1.358
- ('fr_halogen', <class 'numpy.float64'>): z_diff=0.706
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.625
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.489
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.476
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.420

### Representative SMILES
- `CN(C)S(=O)(=O)N(C)c1ccc(C(=O)NCCSCc2ccccc2F)cc1`
- `O=C(COc1ccc(F)cc1)NCc1ccc(N2CCSCC2)cc1`
- `C[C@@H](CNC(=O)c1cccc(F)c1)Oc1cccc(CN2CCSCC2)c1`
- `CCN(CC)Cc1ccc(NC(=O)CSCc2ccc(F)cc2)cc1`
- `CCN(CC)c1ccc(CNC(=O)CSc2ccc(Cl)cc2)cc1`

## Cluster 28
- Count: 51754
- Ratio: 0.0259

### Top RDKit fragments
- fr_Ar_N: presence=1.000, global=0.531, enrich=0.469
- fr_pyridine: presence=0.372, global=0.138, enrich=0.234
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_halogen: presence=0.479, global=0.384, enrich=0.095
- fr_nitrile: presence=0.087, global=0.030, enrich=0.057
- fr_C_S: presence=0.081, global=0.031, enrich=0.050
- fr_alkyl_halide: presence=0.097, global=0.054, enrich=0.044
- fr_morpholine: presence=0.074, global=0.042, enrich=0.031
- fr_tetrazole: presence=0.041, global=0.012, enrich=0.030
- fr_sulfonamd: presence=0.148, global=0.132, enrich=0.016
- fr_oxazole: presence=0.029, global=0.014, enrich=0.015
- fr_hdrzone: presence=0.034, global=0.020, enrich=0.014
- fr_sulfone: presence=0.038, global=0.025, enrich=0.013
- fr_SH: presence=0.014, global=0.003, enrich=0.012
- fr_thiazole: presence=0.072, global=0.062, enrich=0.010
- fr_piperzine: presence=0.070, global=0.062, enrich=0.009
- fr_aldehyde: presence=0.010, global=0.003, enrich=0.007
- fr_Ar_OH: presence=0.041, global=0.034, enrich=0.006
- fr_oxime: presence=0.007, global=0.005, enrich=0.002
- fr_guanido: presence=0.004, global=0.003, enrich=0.001

### Top high RDKit2DNormalized features
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.936
- ('fr_NH0', <class 'numpy.float64'>): z_diff=0.691
- ('fr_pyridine', <class 'numpy.float64'>): z_diff=0.673
- ('fr_nitrile', <class 'numpy.float64'>): z_diff=0.337
- ('fr_C_S', <class 'numpy.float64'>): z_diff=0.292
- ('fr_tetrazole', <class 'numpy.float64'>): z_diff=0.292
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.246
- ('fr_SH', <class 'numpy.float64'>): z_diff=0.235

### Representative SMILES
- `O=c1ccc(-c2ccccc2)nn1CN(C[C@H]1CCOC1)C1CC1`
- `CN(Cc1ccc(Cn2cccn2)cc1)C[C@@H]1CCOC1`
- `C=Cn1ncc(CN(Cc2ccccc2)C[C@H]2CCCO2)c1C`
- `CCN(Cc1ccc(-n2cccn2)cc1)CC1CCOCC1`
- `CC[C@@H](C)Cc1nc(C2CCOCC2)nn1Cc1ccccc1`

## Cluster 29
- Count: 65141
- Ratio: 0.0326

### Top RDKit fragments
- fr_ether: presence=1.000, global=0.514, enrich=0.486
- fr_methoxy: presence=0.653, global=0.230, enrich=0.422
- fr_ester: presence=0.331, global=0.107, enrich=0.223
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_amide: presence=0.913, global=0.701, enrich=0.212
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_NH1: presence=0.749, global=0.642, enrich=0.107
- fr_urea: presence=0.148, global=0.069, enrich=0.079
- fr_morpholine: presence=0.094, global=0.042, enrich=0.052
- fr_imide: presence=0.081, global=0.050, enrich=0.031
- fr_Ndealkylation1: presence=0.090, global=0.061, enrich=0.028
- fr_alkyl_carbamate: presence=0.031, global=0.004, enrich=0.027
- fr_priamide: presence=0.039, global=0.019, enrich=0.020
- fr_piperzine: presence=0.081, global=0.062, enrich=0.020
- fr_unbrch_alkane: presence=0.043, global=0.025, enrich=0.018
- fr_lactone: presence=0.020, global=0.006, enrich=0.015
- fr_benzene: presence=0.847, global=0.834, enrich=0.012
- fr_C_S: presence=0.042, global=0.031, enrich=0.011
- fr_nitrile: presence=0.040, global=0.030, enrich=0.010
- fr_allylic_oxid: presence=0.065, global=0.055, enrich=0.009

### Top high RDKit2DNormalized features
- ('fr_methoxy', <class 'numpy.float64'>): z_diff=1.024
- ('fr_ether', <class 'numpy.float64'>): z_diff=1.014
- ('fr_ester', <class 'numpy.float64'>): z_diff=0.728
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.709
- ('fr_C_O', <class 'numpy.float64'>): z_diff=0.672
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.519
- ('fr_alkyl_carbamate', <class 'numpy.float64'>): z_diff=0.414
- ('fr_urea', <class 'numpy.float64'>): z_diff=0.326

### Representative SMILES
- `COCCOCCC(=O)NCCC(=O)N(C)Cc1ccccc1`
- `CCOc1cc(CN(C(=O)[C@H]2CNC(=O)C2)C2CCCCC2)ccc1OC`
- `CCCC(=O)N(Cc1ccc(OC)cc1)[C@@H](C)C(=O)NCCCOCC`
- `COc1cccc(CC(=O)NCCN(C(C)=O)C2CCOCC2)c1`
- `CCCC(=O)N(Cc1ccc(OC)cc1)[C@@H](C)C(=O)NC[C@@H]1CCCO1`

## Cluster 30
- Count: 67692
- Ratio: 0.0338

### Top RDKit fragments
- fr_Nhpyrrole: presence=1.000, global=0.078, enrich=0.922
- fr_Ar_NH: presence=1.000, global=0.078, enrich=0.922
- fr_Ar_N: presence=1.000, global=0.531, enrich=0.469
- fr_NH1: presence=1.000, global=0.642, enrich=0.358
- fr_bicyclic: presence=0.693, global=0.397, enrich=0.297
- fr_para_hydroxylation: presence=0.372, global=0.216, enrich=0.156
- fr_imidazole: presence=0.207, global=0.054, enrich=0.153
- fr_NH0: presence=0.913, global=0.857, enrich=0.056
- fr_pyridine: presence=0.192, global=0.138, enrich=0.054
- fr_benzene: presence=0.860, global=0.834, enrich=0.026
- fr_methoxy: presence=0.249, global=0.230, enrich=0.018
- fr_Ar_OH: presence=0.049, global=0.034, enrich=0.014
- fr_tetrazole: presence=0.026, global=0.012, enrich=0.014
- fr_sulfide: presence=0.134, global=0.122, enrich=0.012
- fr_phenol_noOrthoHbond: presence=0.025, global=0.024, enrich=0.001
- fr_SH: presence=0.003, global=0.003, enrich=0.000
- fr_phenol: presence=0.025, global=0.025, enrich=0.000
- fr_aldehyde: presence=0.003, global=0.003, enrich=0.000
- fr_diazo: presence=0.000, global=0.000, enrich=0.000
- fr_alkyl_carbamate: presence=0.004, global=0.004, enrich=0.000

### Top high RDKit2DNormalized features
- ('fr_Ar_NH', <class 'numpy.float64'>): z_diff=3.414
- ('fr_Nhpyrrole', <class 'numpy.float64'>): z_diff=3.414
- ('fr_Ar_N', <class 'numpy.float64'>): z_diff=0.963
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.956
- ('fr_imidazole', <class 'numpy.float64'>): z_diff=0.693
- ('fr_bicyclic', <class 'numpy.float64'>): z_diff=0.607
- ('fr_para_hydroxylation', <class 'numpy.float64'>): z_diff=0.364
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.221

### Representative SMILES
- `CCCOC1(CNC(=O)c2ccc3[nH]nnc3c2)CCCCC1`
- `Cn1c(CCNC(=O)c2ccc3c(c2)CCO3)n[nH]c1=S`
- `O=C(NCC[C@@H]1CCCO1)c1ccc2[nH]nnc2c1`
- `O=C(N[C@H]1CC[C@H](Oc2cnccn2)CC1)c1ccc2[nH]ccc2c1`
- `C[C@H]1Cc2cc(C(=O)NCCc3cnn[nH]3)ccc2O1`

## Cluster 31
- Count: 77835
- Ratio: 0.0389

### Top RDKit fragments
- fr_ether: presence=0.599, global=0.514, enrich=0.085
- fr_halogen: presence=0.469, global=0.384, enrich=0.085
- fr_Imine: presence=0.086, global=0.031, enrich=0.055
- fr_phenol_noOrthoHbond: presence=0.070, global=0.024, enrich=0.045
- fr_phenol: presence=0.070, global=0.025, enrich=0.045
- fr_nitrile: presence=0.072, global=0.030, enrich=0.042
- fr_NH2: presence=0.098, global=0.057, enrich=0.042
- fr_Ar_OH: presence=0.073, global=0.034, enrich=0.038
- fr_C_S: presence=0.069, global=0.031, enrich=0.038
- fr_hdrzone: presence=0.054, global=0.020, enrich=0.034
- fr_sulfone: presence=0.058, global=0.025, enrich=0.034
- fr_sulfonamd: presence=0.165, global=0.132, enrich=0.033
- fr_allylic_oxid: presence=0.085, global=0.055, enrich=0.030
- fr_alkyl_halide: presence=0.082, global=0.054, enrich=0.028
- fr_unbrch_alkane: presence=0.051, global=0.025, enrich=0.025
- fr_oxime: presence=0.024, global=0.005, enrich=0.019
- fr_aldehyde: presence=0.017, global=0.003, enrich=0.014
- fr_morpholine: presence=0.056, global=0.042, enrich=0.014
- fr_lactone: presence=0.018, global=0.006, enrich=0.012
- fr_guanido: presence=0.014, global=0.003, enrich=0.011

### Top high RDKit2DNormalized features
- ('fr_Imine', <class 'numpy.float64'>): z_diff=0.335
- ('fr_phenol_noOrthoHbond', <class 'numpy.float64'>): z_diff=0.297
- ('fr_phenol', <class 'numpy.float64'>): z_diff=0.293
- ('fr_azide', <class 'numpy.float64'>): z_diff=0.284
- ('fr_oxime', <class 'numpy.float64'>): z_diff=0.278
- ('fr_aldehyde', <class 'numpy.float64'>): z_diff=0.270
- ('fr_hdrzone', <class 'numpy.float64'>): z_diff=0.256
- ('fr_nitrile', <class 'numpy.float64'>): z_diff=0.250

### Representative SMILES
- `CCCN(CCc1ccccc1)CC[C@@H]1CCOC1`
- `CCN(CC)CCCOc1ccc(C(C)C)cc1`
- `CC(C)[C@H]1C[C@@H](CN(C)Cc2ccccc2)CCO1`
- `CCCN(CCc1ccccc1)CC[C@H]1CCOC1`
- `C[C@@H]1CN(C(C)(C)C)[C@H](c2ccccc2)O1`

## Cluster 32
- Count: 44617
- Ratio: 0.0223

### Top RDKit fragments
- fr_nitro_arom: presence=1.000, global=0.042, enrich=0.958
- fr_nitro: presence=1.000, global=0.048, enrich=0.952
- fr_nitro_arom_nonortho: presence=0.669, global=0.026, enrich=0.643
- fr_amide: presence=0.991, global=0.701, enrich=0.290
- fr_aniline: presence=0.729, global=0.456, enrich=0.273
- fr_C_O_noCOO: presence=1.000, global=0.779, enrich=0.221
- fr_C_O: presence=1.000, global=0.797, enrich=0.203
- fr_benzene: presence=1.000, global=0.834, enrich=0.166
- fr_NH0: presence=1.000, global=0.857, enrich=0.143
- fr_imide: presence=0.188, global=0.050, enrich=0.139
- fr_NH1: presence=0.723, global=0.642, enrich=0.082
- fr_sulfide: presence=0.171, global=0.122, enrich=0.049
- fr_C_S: presence=0.069, global=0.031, enrich=0.038
- fr_ester: presence=0.135, global=0.107, enrich=0.028
- fr_hdrzone: presence=0.047, global=0.020, enrich=0.027
- fr_para_hydroxylation: presence=0.237, global=0.216, enrich=0.022
- fr_barbitur: presence=0.027, global=0.006, enrich=0.021
- fr_hdrzine: presence=0.026, global=0.005, enrich=0.021
- fr_ketone_Topliss: presence=0.073, global=0.054, enrich=0.019
- fr_ketone: presence=0.079, global=0.062, enrich=0.018

### Top high RDKit2DNormalized features
- ('fr_nitro_arom', <class 'numpy.float64'>): z_diff=4.829
- ('fr_nitro', <class 'numpy.float64'>): z_diff=4.516
- ('fr_nitro_arom_nonortho', <class 'numpy.float64'>): z_diff=4.096
- ('fr_benzene', <class 'numpy.float64'>): z_diff=0.691
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.648
- ('fr_imide', <class 'numpy.float64'>): z_diff=0.642
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.564
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.543

### Representative SMILES
- `CCOc1ccc(NC(=O)c2ccc([N+](=O)[O-])cc2[N+](=O)[O-])cc1`
- `C=COc1ccc(NC(=O)c2cc([N+](=O)[O-])cc([N+](=O)[O-])c2)cc1`
- `CCN(CC)CCCOc1ccc(NC(=O)C2(c3ccc([N+](=O)[O-])cc3)CCCCC2)cc1`
- `CCOc1ccc(NC(=O)/C=C/c2ccc([N+](=O)[O-])cc2)c([N+](=O)[O-])c1`
- `O=C(COc1cccc([N+](=O)[O-])c1)Nc1cccc(N2CCCC2)c1`

## Cluster 33
- Count: 57947
- Ratio: 0.0290

### Top RDKit fragments
- fr_sulfonamd: presence=1.000, global=0.132, enrich=0.868
- fr_aryl_methyl: presence=1.000, global=0.398, enrich=0.602
- fr_aniline: presence=0.847, global=0.456, enrich=0.391
- fr_NH1: presence=0.933, global=0.642, enrich=0.292
- fr_amide: presence=0.920, global=0.701, enrich=0.219
- fr_benzene: presence=0.999, global=0.834, enrich=0.165
- fr_C_O_noCOO: presence=0.939, global=0.779, enrich=0.160
- fr_halogen: presence=0.526, global=0.384, enrich=0.142
- fr_C_O: presence=0.939, global=0.797, enrich=0.142
- fr_methoxy: presence=0.267, global=0.230, enrich=0.037
- fr_Ndealkylation1: presence=0.089, global=0.061, enrich=0.027
- fr_unbrch_alkane: presence=0.039, global=0.025, enrich=0.014
- fr_prisulfonamd: presence=0.000, global=0.000, enrich=0.000
- fr_isothiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_thiocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_isocyan: presence=0.000, global=0.000, enrich=-0.000
- fr_diazo: presence=0.000, global=0.000, enrich=-0.000
- fr_benzodiazepine: presence=0.000, global=0.000, enrich=-0.000
- fr_nitroso: presence=0.000, global=0.000, enrich=-0.000
- fr_azide: presence=0.000, global=0.000, enrich=-0.000

### Top high RDKit2DNormalized features
- ('fr_sulfonamd', <class 'numpy.float64'>): z_diff=2.581
- ('fr_aryl_methyl', <class 'numpy.float64'>): z_diff=1.228
- ('fr_benzene', <class 'numpy.float64'>): z_diff=1.056
- ('fr_aniline', <class 'numpy.float64'>): z_diff=0.780
- ('fr_NH1', <class 'numpy.float64'>): z_diff=0.542
- ('fr_amide', <class 'numpy.float64'>): z_diff=0.481
- ('fr_C_O_noCOO', <class 'numpy.float64'>): z_diff=0.344
- ('fr_halogen', <class 'numpy.float64'>): z_diff=0.285

### Representative SMILES
- `CC[C@H](C(=O)N[C@H](C)c1ccc(C)cc1)N(c1ccc(Cl)cc1)S(C)(=O)=O`
- `Cc1cccc(CNC(=O)CN(c2ccc(Br)cc2)S(C)(=O)=O)c1`
- `C[C@H](CCc1ccccc1)NC(=O)c1ccc(Cl)c(N(C)S(C)(=O)=O)c1`
- `CC(=O)Nc1ccc(C)cc1S(=O)(=O)N1CCC[C@@H]1c1ccc(Cl)cc1`
- `Cc1cccc(NC(=O)CCCN(C)S(=O)(=O)c2ccc(Cl)cc2)c1C`
