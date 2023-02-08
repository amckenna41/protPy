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
#
```
Calculate Dipeptide Composition:
```python
dipeptide_comp = protpy.dipeptide_composition(protein_seq)
#
```
Calculate Tripeptide Composition:
```python
tripeptide_comp = protpy.tripeptide_composition(protein_seq)
#
```
Calculate Pseudo Composition:
```python
pseudo_comp = protpy.pseudo_amino_acid_composition(protein_seq)
#
```

Calculate Amphiphilic Composition:
```python
amphiphilic_comp = protpy.amphiphilic_amino_acid_composition(protein_seq)
#
```

## Autocorrelation Descriptors
Calculate MoreauBroto Autocorrelation:
```python
moreaubroto_autocorrelation = protpy.moreaubroto_autocorrelation(protein_seq)
#
```
Calculate Moran Autocorrelation:
```python
moran_autocorrelation = protpy.moran_autocorrelation(protein_seq)
#
```

Calculate Geary Autocorrelation:
```python
geary_autocorrelation = protpy.geary_autocorrelation(protein_seq)
#
```

## Conjoint Triad Descriptors

## CTD

## Sequence Order Correlation Factor

## Sequence Order Coupling Number

## Quasi Sequence Order
