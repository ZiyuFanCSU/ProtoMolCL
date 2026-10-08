# Layer-5 Property-Structural Prototype Cluster Interpretation

This report analyzes each property-structural molecular prototype produced by KMeans. The main diagnostic value is `z_diff_raw`, defined as:

`z_diff_raw = (cluster_mean_raw - global_mean_raw) / global_std_raw`

A positive value means this prototype has a higher descriptor value than the global average; a negative value means it is lower than the global average.

## Overall cluster-size distribution

| Prototype | Count | Ratio |
|---:|---:|---:|
| 0 | 492278 | 0.2461 |
| 1 | 349792 | 0.1749 |
| 2 | 191470 | 0.0957 |
| 3 | 490968 | 0.2455 |
| 4 | 475492 | 0.2377 |

## Prototype 0: sp3-rich, heterocycle-rich

- Count: 492278
- Ratio: 0.2461

### Dominant feature groups

| Feature group | Mean abs z-diff | Top high features | Top low features |
|---|---:|---|---|
| ring_scaffold | 0.514 | NumSaturatedHeterocycles(0.89); NumSaturatedRings(0.88); NumAliphaticRings(0.87) | NumAromaticCarbocycles(-0.27); NumAromaticRings(-0.26); NumAromaticHeterocycles(-0.02) |
| topology_shape | 0.297 | Chi4n(0.64); Chi3n(0.60); Chi2n(0.49) | BalabanJ(-0.60); Kappa3(-0.03); Kappa2(0.05) |
| hbond_flexibility | 0.227 | FractionCSP3(0.49); NumHAcceptors(0.07); NumRotatableBonds(-0.16) | NumHDonors(-0.19); NumRotatableBonds(-0.16); NumHAcceptors(0.07) |
| size_proxy | 0.188 | HeavyAtomCount(0.19) | HeavyAtomCount(0.19) |
| physicochemical | 0.188 | QED(0.38); MolMR(0.15); TPSA(-0.03) | MolLogP(-0.18); TPSA(-0.03); MolMR(0.15) |

### Top enriched descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| NumSaturatedHeterocycles | 1.1107 | 0.5161 | 0.885 |
| NumSaturatedRings | 1.4437 | 0.6854 | 0.885 |
| NumAliphaticRings | 1.7881 | 0.9612 | 0.871 |
| NumAliphaticHeterocycles | 1.3592 | 0.7235 | 0.847 |
| Chi4n | 3.9610 | 3.2080 | 0.645 |
| Chi3n | 5.4793 | 4.6215 | 0.598 |
| FractionCSP3 | 0.4514 | 0.3553 | 0.490 |
| Chi2n | 7.4864 | 6.6081 | 0.489 |
| RingCount | 3.7993 | 3.2513 | 0.475 |
| QED | 0.7129 | 0.6361 | 0.380 |

### Top depleted descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| BalabanJ | 1.5502 | 1.7829 | -0.599 |
| NumAromaticCarbocycles | 1.1740 | 1.4361 | -0.272 |
| NumAromaticRings | 2.0112 | 2.2901 | -0.261 |
| NumHDonors | 0.9060 | 1.0710 | -0.189 |
| MolLogP | 3.0712 | 3.3690 | -0.185 |
| NumRotatableBonds | 4.9103 | 5.3082 | -0.161 |
| TPSA | 70.4028 | 71.3084 | -0.034 |
| Kappa3 | 4.1150 | 4.2314 | -0.032 |
| NumAromaticHeterocycles | 0.8372 | 0.8540 | -0.020 |
| Kappa2 | 8.0621 | 7.9426 | 0.053 |

### Representative molecules closest to the prototype center

- Rank 1, distance=1.5563: `COC(=O)C[C@@H](NC(=O)[C@H]1CC(=O)N(c2ccc3c(c2)CCC3)C1)c1cccs1`
- Rank 2, distance=1.5988: `CC[C@@H](C)NC(=O)[C@@H]1CSCN1C(=O)c1nn(-c2ccccc2)c2c1CCC2`
- Rank 3, distance=1.6123: `Cc1ccccc1C(=O)c1c(NC(=O)CN2CCOCC2)sc2c1CCCC2`
- Rank 4, distance=1.6192: `O=C(NCCc1cccs1)[C@@H]1Cc2cc([N+](=O)[O-])ccc2N2CCCC[C@@H]12`
- Rank 5, distance=1.6418: `CSc1ncccc1C(=O)N1CCc2ccc(NC(=O)[C@@H]3CCCO3)cc2C1`
- Rank 6, distance=1.6431: `COC(=O)c1c(NC(=O)c2ccc(N3CCCC3=O)cc2)sc2c1CCCC2`
- Rank 7, distance=1.6503: `O=C(CN1C(=O)[C@H](C(=O)N2CCCCC2)Sc2ccccc21)NCc1cccs1`
- Rank 8, distance=1.6527: `CCOc1ccc(NC(=O)[C@H]2CCc3sc(C(=O)N4CCOCC4)cc3C2)cc1`

### Brief interpretation

Prototype 0 is mainly characterized by higher NumSaturatedHeterocycles, NumSaturatedRings, NumAliphaticRings, NumAliphaticHeterocycles, Chi4n and lower BalabanJ, NumAromaticCarbocycles, NumAromaticRings, NumHDonors, MolLogP. This suggests a `sp3-rich, heterocycle-rich` molecular property regime.

## Prototype 1: drug-like, sp3-rich, heterocycle-rich, smaller-size

- Count: 349792
- Ratio: 0.1749

### Dominant feature groups

| Feature group | Mean abs z-diff | Top high features | Top low features |
|---|---:|---|---|
| size_proxy | 0.946 | HeavyAtomCount(-0.95) | HeavyAtomCount(-0.95) |
| physicochemical | 0.760 | QED(0.79); TPSA(-0.49); MolLogP(-0.83) | MolMR(-0.94); MolLogP(-0.83); TPSA(-0.49) |
| topology_shape | 0.650 | HallKierAlpha(1.07); BalabanJ(0.27); Kappa3(-0.22) | Chi1(-0.95); Chi0(-0.93); Kappa1(-0.83) |
| hbond_flexibility | 0.569 | FractionCSP3(1.08); NumHDonors(-0.08); NumHAcceptors(-0.51) | NumRotatableBonds(-0.60); NumHAcceptors(-0.51); NumHDonors(-0.08) |
| ring_scaffold | 0.538 | NumSaturatedRings(0.60); NumSaturatedHeterocycles(0.54); NumAliphaticRings(0.42) | NumAromaticRings(-1.08); NumAromaticCarbocycles(-0.92); RingCount(-0.65) |

### Top enriched descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| FractionCSP3 | 0.5678 | 0.3553 | 1.083 |
| HallKierAlpha | -1.5161 | -2.5167 | 1.066 |
| QED | 0.7952 | 0.6361 | 0.787 |
| NumSaturatedRings | 1.2011 | 0.6854 | 0.602 |
| NumSaturatedHeterocycles | 0.8761 | 0.5161 | 0.536 |
| NumAliphaticRings | 1.3578 | 0.9612 | 0.418 |
| NumAliphaticHeterocycles | 0.9776 | 0.7235 | 0.338 |
| NumSaturatedCarbocycles | 0.3250 | 0.1692 | 0.297 |
| BalabanJ | 1.8862 | 1.7829 | 0.266 |
| NumAliphaticCarbocycles | 0.3802 | 0.2377 | 0.230 |

### Top depleted descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| NumAromaticRings | 1.1387 | 2.2901 | -1.079 |
| Chi1 | 9.8181 | 12.8073 | -0.953 |
| HeavyAtomCount | 20.4981 | 26.6259 | -0.946 |
| MolMR | 78.7708 | 102.7736 | -0.939 |
| Chi0 | 14.7427 | 18.9985 | -0.933 |
| NumAromaticCarbocycles | 0.5476 | 1.4361 | -0.923 |
| Kappa1 | 14.6951 | 18.5993 | -0.833 |
| MolLogP | 2.0356 | 3.3690 | -0.827 |
| Kappa2 | 6.2194 | 7.9426 | -0.757 |
| RingCount | 2.4965 | 3.2513 | -0.654 |

### Representative molecules closest to the prototype center

- Rank 1, distance=1.4906: `CNC[C@H]1CCN(C(=O)/C=C\c2c(C)nn(C)c2Cl)C1`
- Rank 2, distance=1.5312: `COC(=O)[C@H](O)C1CCN(Cc2ccc(C)cc2)CC1`
- Rank 3, distance=1.5348: `C=CCNC(=O)[C@H]1CCCN(c2nc(C)cc(C)n2)C1`
- Rank 4, distance=1.5366: `Cc1cc(=O)cc(C(=O)NCC(C)(C)N2CCCC2)o1`
- Rank 5, distance=1.5519: `CNCC1CCN(C(=O)/C=C\c2c(C)nn(C)c2Cl)CC1`
- Rank 6, distance=1.5535: `COCCC[C@H]1CCCN(C(=O)c2ccc(N)nc2)C1`
- Rank 7, distance=1.5580: `CN(C)c1ccc(CN2CCC[C@@H](C(=O)O)C2)cn1`
- Rank 8, distance=1.5596: `COC(=O)[C@@H]1C[C@H](O)CN1C[C@@H](C)c1ccc(F)cc1`

### Brief interpretation

Prototype 1 is mainly characterized by higher FractionCSP3, HallKierAlpha, QED, NumSaturatedRings, NumSaturatedHeterocycles and lower NumAromaticRings, Chi1, HeavyAtomCount, MolMR, Chi0. This suggests a `drug-like, sp3-rich, heterocycle-rich, smaller-size` molecular property regime.

## Prototype 2: larger-size

- Count: 191470
- Ratio: 0.0957

### Dominant feature groups

| Feature group | Mean abs z-diff | Top high features | Top low features |
|---|---:|---|---|
| size_proxy | 1.923 | HeavyAtomCount(1.92) | HeavyAtomCount(1.92) |
| topology_shape | 1.466 | Chi0(1.94); Kappa1(1.93); Chi1(1.91) | HallKierAlpha(-1.35); BalabanJ(-0.49); Kappa3(0.58) |
| physicochemical | 1.419 | MolMR(1.92); MolLogP(1.22); TPSA(0.86) | QED(-1.67); TPSA(0.86); MolLogP(1.22) |
| hbond_flexibility | 0.647 | NumRotatableBonds(1.40); NumHAcceptors(0.84); NumHDonors(0.08) | FractionCSP3(-0.27); NumHDonors(0.08); NumHAcceptors(0.84) |
| ring_scaffold | 0.435 | NumAromaticCarbocycles(1.28); RingCount(1.11); NumAromaticRings(0.97) | NumAromaticHeterocycles(-0.23); NumSaturatedHeterocycles(-0.06); NumSaturatedRings(-0.01) |

### Top enriched descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| Chi0 | 27.8536 | 18.9985 | 1.941 |
| Kappa1 | 27.6469 | 18.5993 | 1.932 |
| HeavyAtomCount | 39.0874 | 26.6259 | 1.923 |
| MolMR | 151.8233 | 102.7736 | 1.918 |
| Chi1 | 18.7957 | 12.8073 | 1.910 |
| Kappa2 | 12.0096 | 7.9426 | 1.787 |
| Chi2n | 9.6897 | 6.6081 | 1.715 |
| Chi3n | 6.8611 | 4.6215 | 1.562 |
| NumRotatableBonds | 8.7675 | 5.3082 | 1.401 |
| Chi4n | 4.8421 | 3.2080 | 1.399 |

### Top depleted descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| QED | 0.2978 | 0.6361 | -1.675 |
| HallKierAlpha | -3.7801 | -2.5167 | -1.345 |
| BalabanJ | 1.5920 | 1.7829 | -0.491 |
| FractionCSP3 | 0.3017 | 0.3553 | -0.273 |
| NumAromaticHeterocycles | 0.6587 | 0.8540 | -0.229 |
| NumSaturatedHeterocycles | 0.4759 | 0.5161 | -0.060 |
| NumSaturatedRings | 0.6786 | 0.6854 | -0.008 |
| NumSaturatedCarbocycles | 0.2028 | 0.1692 | 0.064 |
| NumHDonors | 1.1381 | 1.0710 | 0.077 |
| NumAliphaticCarbocycles | 0.3465 | 0.2377 | 0.176 |

### Representative molecules closest to the prototype center

- Rank 1, distance=1.7011: `Cc1ccc(-n2cc(-c3ccccc3)nc2NC(=O)CN(C[C@@H]2CCCO2)S(=O)(=O)c2ccc(Cl)cc2)cc1F`
- Rank 2, distance=1.7342: `COc1ccc(S(=O)(=O)N(CC(=O)Nc2ccc(N3CCCCC3)cc2)c2cc(C)cc(C)c2)cc1OC`
- Rank 3, distance=1.7388: `CCN(CC)S(=O)(=O)c1ccc(N2CCOCC2)c(NC(=O)c2cc(-c3ccccc3)nn2-c2ccccc2)c1`
- Rank 4, distance=1.7405: `CC(C)c1ccc(N(C(=O)Cn2nnc3ccccc32)[C@@H](C(=O)NC[C@@H]2CCCO2)c2ccccc2F)cc1`
- Rank 5, distance=1.7409: `CC(C)c1ccc([C@H](C(=O)NC[C@@H]2CCCO2)N(C(=O)Cn2nnc3ccccc32)c2ccccc2F)cc1`
- Rank 6, distance=1.7416: `CC(C)c1ccc([C@H](C(=O)NC[C@H]2CCCO2)N(C(=O)Cn2nnc3ccccc32)c2ccc(F)cc2)cc1`
- Rank 7, distance=1.7416: `CC(C)c1ccc([C@H](C(=O)NC[C@@H]2CCCO2)N(C(=O)Cn2nnc3ccccc32)c2ccc(F)cc2)cc1`
- Rank 8, distance=1.7570: `O=C(C[C@H](c1cn(Cc2ccccc2)c2ccc([N+](=O)[O-])cc12)c1ccc(Cl)cc1)NCCN1CCOCC1`

### Brief interpretation

Prototype 2 is mainly characterized by higher Chi0, Kappa1, HeavyAtomCount, MolMR, Chi1 and lower QED, HallKierAlpha, BalabanJ, FractionCSP3, NumAromaticHeterocycles. This suggests a `larger-size` molecular property regime.

## Prototype 3: polar/H-bond-rich, drug-like, aromatic-ring-rich, heterocycle-rich

- Count: 490968
- Ratio: 0.2455

### Dominant feature groups

| Feature group | Mean abs z-diff | Top high features | Top low features |
|---|---:|---|---|
| size_proxy | 0.794 | HeavyAtomCount(-0.79) | HeavyAtomCount(-0.79) |
| topology_shape | 0.719 | BalabanJ(0.85); HallKierAlpha(0.35); Kappa3(-0.15) | Chi3n(-1.00); Chi4n(-0.97); Chi2n(-0.97) |
| ring_scaffold | 0.482 | NumAromaticHeterocycles(0.03); NumAromaticRings(-0.17); NumAromaticCarbocycles(-0.22) | RingCount(-0.80); NumAliphaticRings(-0.78); NumSaturatedRings(-0.75) |
| physicochemical | 0.432 | QED(0.43); TPSA(-0.24); MolLogP(-0.27) | MolMR(-0.77); MolLogP(-0.27); TPSA(-0.24) |
| hbond_flexibility | 0.336 | NumHDonors(0.14); NumRotatableBonds(-0.32); NumHAcceptors(-0.38) | FractionCSP3(-0.50); NumHAcceptors(-0.38); NumRotatableBonds(-0.32) |

### Top enriched descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| BalabanJ | 2.1140 | 1.7829 | 0.851 |
| QED | 0.7240 | 0.6361 | 0.435 |
| HallKierAlpha | -2.1860 | -2.5167 | 0.352 |
| NumHDonors | 1.1929 | 1.0710 | 0.139 |
| NumAromaticHeterocycles | 0.8802 | 0.8540 | 0.031 |
| Kappa3 | 3.6807 | 4.2314 | -0.150 |
| NumAromaticRings | 2.1083 | 2.2901 | -0.170 |
| NumAromaticCarbocycles | 1.2281 | 1.4361 | -0.216 |
| TPSA | 64.8796 | 71.3084 | -0.245 |
| MolLogP | 2.9265 | 3.3690 | -0.274 |

### Top depleted descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| Chi3n | 3.1897 | 4.6215 | -0.999 |
| Chi4n | 2.0776 | 3.2080 | -0.968 |
| Chi2n | 4.8711 | 6.6081 | -0.967 |
| Chi1 | 10.2717 | 12.8073 | -0.809 |
| RingCount | 2.3274 | 3.2513 | -0.800 |
| HeavyAtomCount | 21.4844 | 26.6259 | -0.794 |
| NumAliphaticRings | 0.2191 | 0.9612 | -0.782 |
| MolMR | 83.0185 | 102.7736 | -0.773 |
| NumSaturatedRings | 0.0392 | 0.6854 | -0.754 |
| Chi0 | 15.5654 | 18.9985 | -0.753 |

### Representative molecules closest to the prototype center

- Rank 1, distance=0.7999: `C/C(=N/NS(=O)(=O)c1ccc(C(C)C)cc1)c1ccncc1`
- Rank 2, distance=0.8554: `COC(=O)c1cccc(C(=O)NCc2cc(C)cc(Cl)n2)c1`
- Rank 3, distance=0.8647: `C[C@@H](N[C@@H](C)c1ccccc1[N+](=O)[O-])c1cncc(F)c1`
- Rank 4, distance=0.8731: `C[C@@H](NC(=O)[C@@H](C)c1ccc([N+](=O)[O-])cc1)c1cccs1`
- Rank 5, distance=0.8954: `CC(=O)c1cccc(NC(=O)/C=C/c2c(C)nn(C)c2C)c1`
- Rank 6, distance=0.9081: `Cc1ccc(C(=O)[C@@H](C#N)C(=O)Nc2cccc(C)c2)s1`
- Rank 7, distance=0.9739: `Cc1ccc(NC(=O)[C@@H](C)OC(=O)c2cccnc2Cl)cc1`
- Rank 8, distance=0.9954: `CC(=O)Nc1ccc(/C=C/C(=O)c2cn(C)nc2C)cc1`

### Brief interpretation

Prototype 3 is mainly characterized by higher BalabanJ, QED, HallKierAlpha, NumHDonors, NumAromaticHeterocycles and lower Chi3n, Chi4n, Chi2n, Chi1, RingCount. This suggests a `polar/H-bond-rich, drug-like, aromatic-ring-rich, heterocycle-rich` molecular property regime.

## Prototype 4: lipophilic, aromatic-ring-rich, larger-size

- Count: 475492
- Ratio: 0.2377

### Dominant feature groups

| Feature group | Mean abs z-diff | Top high features | Top low features |
|---|---:|---|---|
| physicochemical | 0.550 | MolLogP(0.59); MolMR(0.56); TPSA(0.30) | QED(-0.75); TPSA(0.30); MolMR(0.56) |
| size_proxy | 0.547 | HeavyAtomCount(0.55) | HeavyAtomCount(0.55) |
| ring_scaffold | 0.476 | NumAromaticRings(0.85); NumAromaticCarbocycles(0.67); RingCount(0.37) | NumSaturatedRings(-0.58); NumSaturatedHeterocycles(-0.53); NumAliphaticRings(-0.51) |
| hbond_flexibility | 0.378 | NumRotatableBonds(0.38); NumHAcceptors(0.37); NumHDonors(0.08) | FractionCSP3(-0.69); NumHDonors(0.08); NumHAcceptors(0.37) |
| topology_shape | 0.350 | Chi1(0.55); Chi0(0.54); Kappa1(0.46) | HallKierAlpha(-0.72); BalabanJ(-0.26); Chi4n(0.06) |

### Top enriched descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| NumAromaticRings | 3.1974 | 2.2901 | 0.850 |
| NumAromaticCarbocycles | 2.0813 | 1.4361 | 0.670 |
| MolLogP | 4.3200 | 3.3690 | 0.590 |
| MolMR | 117.0201 | 102.7736 | 0.557 |
| Chi1 | 14.5261 | 12.8073 | 0.548 |
| HeavyAtomCount | 30.1673 | 26.6259 | 0.547 |
| Chi0 | 21.4560 | 18.9985 | 0.539 |
| Kappa1 | 20.7568 | 18.5993 | 0.461 |
| Kappa2 | 8.9100 | 7.9426 | 0.425 |
| NumRotatableBonds | 6.2346 | 5.3082 | 0.375 |

### Top depleted descriptors

| Feature | Cluster mean | Global mean | z-diff |
|---|---:|---:|---:|
| QED | 0.4850 | 0.6361 | -0.748 |
| HallKierAlpha | -3.1884 | -2.5167 | -0.715 |
| FractionCSP3 | 0.2205 | 0.3553 | -0.687 |
| NumSaturatedRings | 0.1908 | 0.6854 | -0.577 |
| NumSaturatedHeterocycles | 0.1613 | 0.5161 | -0.528 |
| NumAliphaticRings | 0.4775 | 0.9612 | -0.509 |
| NumAliphaticHeterocycles | 0.3974 | 0.7235 | -0.434 |
| NumSaturatedCarbocycles | 0.0294 | 0.1692 | -0.267 |
| BalabanJ | 1.6837 | 1.7829 | -0.255 |
| NumAliphaticCarbocycles | 0.0801 | 0.2377 | -0.254 |

### Representative molecules closest to the prototype center

- Rank 1, distance=1.2507: `CCCCOC(=O)c1ccc(NC(=O)Cn2ccc(=O)c3c(C)cc(C)cc32)cc1`
- Rank 2, distance=1.3040: `CC(=O)Nc1ccc(C(=O)[C@H](C)OC(=O)Cc2coc3c2ccc(C)c3C)cc1`
- Rank 3, distance=1.3397: `CCc1ccc(NC(=O)[C@H](C)OC(=O)c2c(C)nn(Cc3ccccc3)c2Cl)cc1`
- Rank 4, distance=1.3572: `C=CCN1COc2c3sc(NC(=O)c4ccccc4)c(C(=O)OCC)c3ccc2C1`
- Rank 5, distance=1.3635: `Cc1oc(-c2ccccc2)nc1COC(=O)[C@@H](NC(=O)c1ccccc1Cl)C(C)C`
- Rank 6, distance=1.3646: `CCS(=O)(=O)Nc1ccc(C2=NN(C(=O)c3cccs3)[C@@H](c3ccc(F)cc3)C2)cc1`
- Rank 7, distance=1.3651: `COc1ccc(CN/C=C2/C(=O)c3sccc3N(Cc3ccc(C)cc3)S2(=O)=O)cc1`
- Rank 8, distance=1.3661: `Cc1cc(C(=O)COC(=O)c2ccc(O)cc2)c(C)n1CCc1ccc(F)cc1`

### Brief interpretation

Prototype 4 is mainly characterized by higher NumAromaticRings, NumAromaticCarbocycles, MolLogP, MolMR, Chi1 and lower QED, HallKierAlpha, FractionCSP3, NumSaturatedRings, NumSaturatedHeterocycles. This suggests a `lipophilic, aromatic-ring-rich, larger-size` molecular property regime.

