
# protpy - Package for generating protein physicochemical, biochemical and structural descriptors using their constituent amino acids. #
[![PyPI](https://img.shields.io/pypi/v/protpy)](https://pypi.org/project/protpy/)
[![pytest](https://github.com/amckenna41/protpy/workflows/Building%20and%20Testing/badge.svg)](https://github.com/amckenna41/protpy/actions?query=workflowBuilding%20and%20Testing)
[![Platforms](https://img.shields.io/badge/platforms-linux%2C%20macOS%2C%20Windows-green)](https://pypi.org/project/protpy/)
[![PythonV](https://img.shields.io/pypi/pyversions/protpy?logo=2)](https://pypi.org/project/protpy/)
[![codecov](https://codecov.io/gh/amckenna41/protPy/graph/badge.svg?token=H508M9J13W)](https://codecov.io/gh/amckenna41/protPy)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)
[![Issues](https://img.shields.io/github/issues/amckenna41/protpy)](https://github.com/amckenna41/protpy/issues)
<!-- [![Size](https://img.shields.io/github/repo-size/amckenna41/protpy)](https://github.com/amckenna41/protpy)
[![Commits](https://img.shields.io/github/commit-activity/w/amckenna41/protpy)](https://github.com/amckenna41/protpy) -->
<!-- <p align="center">
<img src="https://images.newscientist.com/wp-content/uploads/2021/07/22155326/22-july_deepmind-proteome.jpg?width=300" alt="protpyLogo" height="200"/>
</p> -->

* 🧬 A <b>demo</b> of the software is available [here][demo]
* 📝 A <b>Medium</b> article about `protPy` and its background is available [here][article]

<p align="center">
  <img src="https://raw.githubusercontent.com/amckenna41/protpy/main/images/protpy.png" alt="protpy" width="600"/>
</p>

Table of Contents
-----------------

  * [Introduction](#introduction)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Documentation](#documentation)
  * [Directories](#directories)
  * [Tests](#tests)
  * [Issues](#Issues)
  * [Contact](#contact)
  * [References](#references)

Introduction
------------
`protpy` is a Python software package for generating a variety of physicochemical, biochemical and structural descriptors for proteins. All of these descriptors are calculated using sequence-derived or physicochemical features of the amino acids that make up the proteins. These descriptors have been highly studied and used in a series of Bioinformatic applications including protein engineering, SAR (sequence-activity-relationships), predicting protein structure & function, subcellular localization, protein-protein interactions, drug-target interactions etc. 

This software is aimed at any researcher or developer using protein sequence/structural data, and was mainly created to use in my own project [`pySAR`](https://github.com/amckenna41/pySAR) which uses protein sequence data to identify Sequence Activity Relationships (SAR) using Machine Learning [[1]](#references). `protpy` is built and developed in Python 3.10.

The descriptors available in `protpy` include:

<details>
<summary><b>Composition Descriptors (22)</b></summary>

* **Amino Acid Composition (AAComp)**
* **Dipeptide Composition (DPComp)**
* **Tripeptide Composition (TPComp)**
* **Grand Average of Hydropathy (GRAVY)**
* **Aromaticity**
* **Instability Index**
* **Isoelectric Point**
* **Molecular Weight**
* **Charge Distribution**
* **Hydrophobic/Polar/Charged Composition (HPC)**
* **Secondary Structure Propensity (SSP)**
* **k-mer Composition**
* **Reduced Alphabet Composition**
* **Motif Composition**
* **Amino Acid Pair Composition**
* **Aliphatic Index**
* **Extinction Coefficient**
* **Boman Index**
* **Aggregation Propensity**
* **Hydrophobic Moment**
* **Shannon Entropy**
* **Pseudo Amino Acid Composition (PAAComp)**
* **Amphiphilic Amino Acid Composition (APAAComp)**

</details>

<details>
<summary><b>Autocorrelation Descriptors (3)</b></summary>

* **Moreaubroto Autocorrelation (MBAuto)**
* **Moran Autocorrelation (MAuto)**
* **Geary Autocorrelation (GAuto)**

</details>

<details>
<summary><b>Conjoint Triad (1)</b></summary>

* **Conjoint Triad (CTriad)**

</details>

<details>
<summary><b>CTD Descriptors (4)</b></summary>

* **CTD Composition**
* **CTD Transition**
* **CTD Distribution**
* **CTD Combined**

</details>

<details>
<summary><b>Sequence Order Descriptors (5)</b></summary>

* **Sequence Order Coupling Number — single (SOCN)**
* **Sequence Order Coupling Number — series**
* **Sequence Order Coupling Number — all matrices**
* **Quasi Sequence Order (QSO)**
* **Quasi Sequence Order — all matrices**

</details>

More detail of each descriptor is listed in the markdown file: [`DESCRIPTORS.md`](https://github.com/amckenna41/protpy/DESCRIPTORS.md)


Requirements
------------
* [Python][python] >= 3.9
* [aaindex][aaindex] >= 1.2.0
* [numpy][numpy] >= 2.4.4
* [pandas][pandas] >= 3.0.2
* [varname][varname] >= 0.15.1
* [biopython][biopython] >= 1.87 (only required for testing) 

Installation 
------------
Install the latest version of `protpy` using pip:

```bash
pip3 install protpy --upgrade
```

Install by cloning repository:

```bash
git clone https://github.com/amckenna41/protpy.git
python3 setup.py install
```

Usage
-----
### Import `protpy` after installation: 
```python
import protpy as protpy
```

### Import protein sequence from fasta:
```python
from Bio import SeqIO

with open("test_fasta.fasta") as pro:
    protein_seq = str(next(SeqIO.parse(pro,'fasta')).seq)
```


<details>
<summary><b>Composition Descriptors Usage Examples</b></summary>

Calculate Amino Acid Composition:
```python
amino_acid_composition = protpy.amino_acid_composition(protein_seq)
# A      C      D      E      F ...
# 6.693  3.108  5.817  3.347  6.614 ...
```

Calculate Dipeptide Composition:
```python
dipeptide_composition = protpy.dipeptide_composition(protein_seq)
# AA    AC    AD   AE    AF ...
# 0.72  0.16  0.48  0.4  0.24 ...
```

Calculate Tripeptide Composition:
```python
tripeptide_composition = protpy.tripeptide_composition(protein_seq)
# AAA  AAC  AAD  AAE  AAF ...
# 1    0    0    2    0 ...
```

Calculate GRAVY (Grand Average of Hydropathy):
```python
gravy = protpy.gravy(protein_seq)
# GRAVY
# -0.045
```

Calculate Aromaticity:
```python
aromaticity = protpy.aromaticity(protein_seq)
# Aromaticity
# 0.118
```

Calculate Instability Index:
```python
instability = protpy.instability_index(protein_seq)
# InstabilityIndex
# 31.836
```

Calculate Isoelectric Point:
```python
pi = protpy.isoelectric_point(protein_seq)
# IsoelectricPoint
# 5.412
```

Calculate Molecular Weight:
```python
mw = protpy.molecular_weight(protein_seq)
# MolecularWeight (Da)
# 139122.355
```

Calculate Charge Distribution:
```python
charge = protpy.charge_distribution(protein_seq)
#using default parameters: ph=7.4

# PositiveCharge  NegativeCharge  NetCharge
# 99.526          114.956         -15.43
```

Calculate Hydrophobic/Polar/Charged Composition:
```python
hpc = protpy.hydrophobic_polar_charged_composition(protein_seq)
# Hydrophobic  Polar   Charged
# 44.542       32.669  18.247
```

Calculate Secondary Structure Propensity:
```python
ssp = protpy.secondary_structure_propensity(protein_seq)
# Helix  Sheet  Coil
# 0.983  1.05   1.043
```

Calculate k-mer Composition:
```python
kmer = protpy.kmer_composition(protein_seq)
#using default parameters: k=2

# AA     AC     AD  ...
# 0.797  0.159  ... ...
```

Calculate Reduced Alphabet Composition:
```python
reduced = protpy.reduced_alphabet_composition(protein_seq)
#using default parameters: alphabet_size=6

# Group_1  Group_2  Group_3  Group_4  Group_5  Group_6
# 25.339   34.741   9.163    9.084    10.837   10.837
```

Calculate Motif Composition:
```python
motif = protpy.motif_composition(protein_seq)
# NxST_glycosylation  RGD_integrin  KDEL_retention  ...
# 23                  0             0               ...
```

Calculate Amino Acid Pair Composition:
```python
aapair = protpy.amino_acid_pair_composition(protein_seq)
# AA_Hydrophobic-Hydrophobic  AA_Hydrophobic-Polar  ...
# 0.797                       0.159                 ...
```

Calculate Aliphatic Index:
```python
aliphatic = protpy.aliphatic_index(protein_seq)
# AliphaticIndex
# 82.725
```

Calculate Extinction Coefficient:
```python
extinction = protpy.extinction_coefficient(protein_seq)
# ExtCoeff_Reduced  ExtCoeff_Oxidized
# 140960            143335
```

Calculate Boman Index:
```python
boman = protpy.boman_index(protein_seq)
# BomanIndex
# 0.119
```

Calculate Aggregation Propensity:
```python
aggregation = protpy.aggregation_propensity(protein_seq)
# AggregProneRegions  AggregProneFraction
# 58                  11.793
```

Calculate Hydrophobic Moment:
```python
hm = protpy.hydrophobic_moment(protein_seq)
#using default parameters: window=11, angle=100

# HydrophobicMoment_Mean  HydrophobicMoment_Max
# 0.272                   0.813
```

Calculate Shannon Entropy:
```python
se = protpy.shannon_entropy(protein_seq)
# ShannonEntropy
# 4.163
```

Calculate Pseudo Composition:
```python
pseudo_composition = protpy.pseudo_amino_acid_composition(protein_seq)
#using default parameters: lamda=30, weight=0.05, properties=[]

# PAAC_1  PAAC_2  PAAC_3  PAAC_4  PAAC_5 ...
# 0.127   0.059   0.111   0.064   0.126 ...
```

Calculate Amphiphilic Composition:
```python
amphiphilic_composition = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq)
#using default parameters: lamda=30, weight=0.5, properties=[hydrophobicity_, hydrophilicity_]

# APAAC_1  APAAC_2  APAAC_3  APAAC_4  APAAC_5 ...
# 6.624    3.076    5.757    3.032    5.988 ...
```

</details>

<details>
<summary><b>Autocorrelation Descriptors Usage Examples</b></summary>

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

</details>

<details>
<summary><b>Conjoint Triad Descriptors Usage Examples</b></summary>

Calculate Conjoint Triad:
```python
conjoint_triad = protpy.conjoint_triad(protein_seq)
# 111  112  113  114  115 ...
# 7    17   11   3    6 ...
```

</details>

<details>
<summary><b>CTD Descriptors Usage Examples</b></summary>

Calculate CTD:
```python
ctd = protpy.ctd(protein_seq)
#using default parameters: property="hydrophobicity", all_ctd=True

# hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  hydrophobicity_CTD_C_03  normalized_vdwv_CTD_C_01 ...
# 0.279                    0.386                    0.335                    0.389 ...                   
```

</details>

<details>
<summary><b>Sequence Order Descriptors Usage Examples</b></summary>

Calculate Sequence Order Coupling Number (SOCN):
```python
socn = protpy.sequence_order_coupling_number_(protein_seq)
#using default parameters: d=1, distance_matrix="schneider-wrede"

#401.387        
```

Calculate all SOCN's per distance matrix:
```python
#using default parameters: lag=30, distance_matrix="schneider-wrede"
socn_all = protpy.sequence_order_coupling_number(protein_seq)

# SOCN_SW1  SOCN_SW2  SOCN_SW3  SOCN_SW4  SOCN_SW5 ...
# 401.387    409.243    376.946    393.042    396.196 ...  

#using custom parameters: lag=10, distance_matrix="grantham"
socn_all = protpy.sequence_order_coupling_number(protein_seq, lag=10, distance_matrix="grantham")      

# SOCN_Grant1  SOCN_Grant_2  SOCN_Grant_3  SOCN_Grant_4  SOCN_Grant_5 ...
# 399.125    402.153    387.820    393.111    409.096 ...  
```

Calculate Quasi Sequence Order (QSO):
```python
#using default parameters: lag=30, weight=0.1, distance_matrix="schneider-wrede"
qso = protpy.quasi_sequence_order(protein_seq)

# QSO_SW1   QSO_SW2   QSO_SW3   QSO_SW4   QSO_SW5 ...
# 0.005692  0.002643  0.004947  0.002846  0.005625 ...  

#using custom parameters: lag=10, weight=0.2, distance_matrix="grantham"
qso = protpy.quasi_sequence_order(protein_seq, lag=10, weight=0.2, distance_matrix="grantham")

# QSO_Grant1   QSO_Grant2   QSO_Grant3   QSO_Grant4   QSO_Grant5 ...
# 0.123287  0.079967  0.04332  0.039983  0.013332 ...  
```

</details>

Documentation
-------------
The documentation for `protpy` is hosted on **ReadTheDocs** and is available [here][readthedocs]. It is automatically built and published from the repository and includes full API reference, usage examples, and a changelog.

Directories
-----------
* `/tests` - unit and integration tests for `protpy` package.
* `/protpy` - source code and all required external data files for package.
* `/docs` - `protpy` documentation.

Tests
-----
To run all tests, from the main `protpy` folder run:
```
python3 -m unittest discover tests -v
-v: verbose output flag
```

Contact
-------
If you have any questions or comments, please contact amckenna41@qub.ac.uk or raise an issue on the [Issues][Issues] tab.

References
----------
\[1\]: Mckenna, A., & Dubey, S. (2022). Machine learning based predictive model for the analysis of sequence activity relationships using protein spectra and protein descriptors. Journal of Biomedical Informatics, 128(104016), 104016. https://doi.org/10.1016/j.jbi.2022.104016 <br>
\[2\]: Shuichi Kawashima, Minoru Kanehisa, AAindex: Amino Acid index database, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Page 374, https://doi.org/10.1093/nar/28.1.374 <br>
\[3\]: Dong, J., Yao, ZJ., Zhang, L. et al. PyBioMed: a python library for various molecular representations of chemicals, proteins and DNAs and their interactions. J Cheminform 10, 16 (2018). https://doi.org/10.1186/s13321-018-0270-2 <br>
\[4\]: Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
fold class predictions. Nucleic Acids Res, 22, 3616-3619. <br>
\[5\]: Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
subcellular localization prediction. Bioinformatics, 17, 721-728. <br>
\[6\]: Broto P, Moreau G, Vandicke C: Molecular structures: perception,
autocorrelation descriptor and SAR studies. Eur J Med Chem 1984, 19: 71–78. <br>
\[7\]: Ong, S.A., Lin, H.H., Chen, Y.Z. et al. Efficacy of different protein
descriptors in predicting protein functional families. BMC Bioinformatics
8, 300 (2007). https://doi.org/10.1186/1471-2105-8-300 <br>
\[8\]: Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim. Prediction
of protein folding class using global description of amino acid sequence.
Proc.Natl. Acad.Sci.USA, 1995, 92, 8700-8704. <br>
\[9\]: Juwen Shen, Jian Zhang, Xiaomin Luo, Weiliang Zhu, Kunqian Yu, Kaixian Chen,
Yixue Li, Huanliang Jiang. Predicting proten-protein interactions based only
on sequences inforamtion. PNAS. 2007 (104) 4337-4341. <br>
\[10\]: Kuo-Chen Chou. Prediction of Protein Subcellar Locations by Incorporating
Quasi-Sequence-Order Effect. Biochemical and Biophysical Research
Communications 2000, 278, 477-483. <br>
\[11\]: Kuo-Chen Chou. Prediction of Protein Cellular Attributes Using
Pseudo-Amino Acid Composition. PROTEINS: Structure, Function, and
Genetics, 2001, 43: 246-255. <br>
\[12\]: Kuo-Chen Chou. Using amphiphilic pseudo amino acid composition to predict enzyme
subfamily classes. Bioinformatics, 2005,21,10-19.

Support
-------
[<img src="https://img.shields.io/github/stars/amckenna41/protpy?color=green&label=star%20it%20on%20GitHub" width="132" height="20" alt="Star it on GitHub">](https://github.com/amckenna41/protpy)


<a href="https://www.buymeacoffee.com/amckenna41" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>

[Back to top](#TOP)

<!-- Links -->
[python]: https://www.python.org/downloads/release/python-360/
[protpy]: https://github.com/amckenna41/protpy
[aaindex]: https://github.com/amckenna41/aaindex
[varname]: https://pypi.org/project/varname/
[biopython]: https://biopython.org/
[numpy]: https://numpy.org/
[pandas]: https://pandas.pydata.org/
[PyPi]: https://pypi.org/project/protpy/
[demo]: https://colab.research.google.com/drive/12E3ayovpZOf6Gv-8ILwkN1zReKslEksB?usp=sharing
[article]: https://medium.com/@ajmckenna69/protpy-05f02f821baa
[Issues]: https://github.com/amckenna41/protpy/issues
[readthedocs]: https://protpy.readthedocs.io/en/latest