# Descriptors

Modules
-------
* `autocorrelation.py` - Autocorrelation descriptors including MoreauBroto, Moran and Geary Autocorrelations.
* `composition.py` - Composition descriptors including Amino Acid, Dipeptide, Tripeptide, Pseudo, Amphipathic Compositions.
* `conjoint_triad.py` - Conjoint Triad descriptor.
* `ctd.py` - Composition, Transition and Distribution (CTD) descriptors.
* `sequence_order.py` - Sequence Order related descriptors including Quasi Sequence Order and Sequence Order Coupling Number.

Usage
-----
## Import `protpy` after installation: 
```python
import protpy as protpy
```

# Import protein sequence from fasta:
```python
from Bio import SeqIO

with open("test_fasta.fasta") as pro:
    protein_seq = str(next(SeqIO.parse(pro,'fasta')).seq)
```

## Composition Descriptors
Calculate Amino Acid Composition:
```python
amino_acid_comp = protpy.amino_acid_composition(protein_seq)
# A      C      D      E      F ...
# 6.693  3.108  5.817  3.347  6.614 ...
```

Calculate Dipeptide Composition:
```python
dipeptide_comp = protpy.dipeptide_composition(protein_seq)
# AA    AC    AD   AE    AF ...
# 0.72  0.16  0.48  0.4  0.24 ...
```

Calculate Tripeptide Composition:
```python
tripeptide_comp = protpy.tripeptide_composition(protein_seq)
# AAA  AAC  AAD  AAE  AAF ...
# 1    0    0    2    0 ...
```

Calculate Pseudo Composition:
```python
pseudo_comp = protpy.pseudo_amino_acid_composition(protein_seq) 
#using default parameters: lamda=30, weight=0.05, properties=[]

# PAAC_1  PAAC_2  PAAC_3  PAAC_4  PAAC_5 ...
# 0.127        0.059        0.111        0.064        0.126 ...
```

Calculate Amphiphilic Composition:
```python
amphiphilic_comp = protpy.amphiphilic_amino_acid_composition(protein_seq)
#using default parameters: lamda=30, weight=0.5, properties=[hydrophobicity_, hydrophilicity_]

# APAAC_1  APAAC_2  APAAC_3  APAAC_4  APAAC_5 ...
# 6.06    2.814    5.267     3.03    5.988 ...
```

## Autocorrelation Descriptors
Calculate MoreauBroto Autocorrelation:
```python
moreaubroto_autocorrelation = protpy.moreaubroto_autocorrelation(protein_seq)
#using default parameters: lag=30, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"], normalize=True

# MBAuto_CIDH920105_1  MBAuto_CIDH920105_2  MBAuto_CIDH920105_3  MBAuto_CIDH920105_4  MBAuto_CIDH920105_5 ...  
# -0.052               -0.104               -0.156               -0.208               0.246 ...
```

Calculate Moran Autocorrelation:
```python
moran_autocorrelation = protpy.moran_autocorrelation(protein_seq)
#using default parameters: lag=30, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"], normalize=True

# MAuto_CIDH920105_1  MAuto_CIDH920105_2  MAuto_CIDH920105_3  MAuto_CIDH920105_4  MAuto_CIDH920105_5 ...
# -0.07786            -0.07879            -0.07906            -0.08001            0.14911 ...
```

Calculate Geary Autocorrelation:
```python
geary_autocorrelation = protpy.geary_autocorrelation(protein_seq)
#using default parameters: lag=30, properties=["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"], normalize=True

# GAuto_CIDH920105_1  GAuto_CIDH920105_2  GAuto_CIDH920105_3  GAuto_CIDH920105_4  GAuto_CIDH920105_5 ...
# 1.057               1.077               1.04                1.02                1.013 ...
```

## Conjoint Triad Descriptors
Calculate Conjoint Triad:
```python
conj_triad = protpy.conjoint_triad(protein_seq)
# 111  112  113  114  115 ...
# 7    17   11   3    6 ...
```

## CTD Descriptors
Calculate CTD:
```python
ctd = protpy.ctd(protein_seq)
#using default parameters: property="hydrophobicity", all_ctd=True

# hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  hydrophobicity_CTD_C_03  normalized_vdwv_CTD_C_01 ...
# 0.279                    0.386                    0.335                    0.389 ...                   
```

## Sequence Order Descriptors 
Calculate Sequence Order Coupling Number (SOCN):
```python
socn = protpy.sequence_order_coupling_number_(protein_seq)
#using default parameters: d=1, distance_matrix="schneider-wrede-physiochemical-distance-matrix"

#401.387        
```

Calculate All SOCN per distance matrix:
```python
socn_all = protpy.sequence_order_coupling_number(protein_seq)
#using default parameters: lag=30, distance_matrix="schneider-wrede-physiochemical-distance-matrix.json"

# SOCN_SW_1  SOCN_SW_2  SOCN_SW_3  SOCN_SW_4  SOCN_SW_5 ...
# 401.387    409.243    376.946    393.042    396.196 ...        
```

Calculate Quasi Sequence Order (QSO):
```python
qso = protpy.quasi_sequence_order(protein_seq)
#using default parameters: lag=30, weight=0.1, distance_matrix="schneider-wrede-physiochemical-distance-matrix.json"

# QSO_SW1   QSO_SW2   QSO_SW3   QSO_SW4   QSO_SW5 ...
# 0.005692  0.002643  0.004947  0.002846  0.005625 ...        
```