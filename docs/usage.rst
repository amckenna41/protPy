Usage
=====

.. _installation:

Installation
------------

Install the latest version of **protpy** via ``pip``:

.. code-block:: console

   pip3 install protpy --upgrade

Alternatively, clone the repository and install from source:

.. code-block:: console

   git clone https://github.com/amckenna41/protpy.git
   cd protpy
   pip install -e .

Importing protpy
----------------

.. code-block:: python

   import protpy as protpy

Loading a protein sequence from a FASTA file
--------------------------------------------

.. code-block:: python

   from Bio import SeqIO

   with open("protein.fasta") as f:
       protein_seq = str(next(SeqIO.parse(f, "fasta")).seq)

Composition Descriptors
-----------------------

Calculate **Amino Acid Composition (AAComp)**:

.. code-block:: python

   amino_acid_comp = protpy.amino_acid_composition(protein_seq)
   # A      C      D      E      F  ...
   # 6.693  3.108  5.817  3.347  6.614 ...

Calculate **Dipeptide Composition (DPComp)**:

.. code-block:: python

   dipeptide_comp = protpy.dipeptide_composition(protein_seq)
   # AA    AC    AD    AE    AF  ...
   # 0.72  0.16  0.48  0.4   0.24 ...

Calculate **Tripeptide Composition (TPComp)**:

.. code-block:: python

   tripeptide_comp = protpy.tripeptide_composition(protein_seq)
   # AAA  AAC  AAD  AAE  AAF ...
   # 1    0    0    2    0 ...

Calculate **Pseudo Amino Acid Composition (PAAComp)**:

.. code-block:: python

   #using default parameters: lamda=30, weight=0.05, properties=[]
   pseudo_comp = protpy.pseudo_amino_acid_composition(protein_seq)
   # PAAC_1  PAAC_2  PAAC_3  PAAC_4  PAAC_5 ...
   # 0.127   0.059   0.111   0.064   0.126 ...

Calculate **Amphiphilic Pseudo Amino Acid Composition (APAAComp)**:

.. code-block:: python

   #using default parameters: lamda=30, weight=0.5
   amphiphilic_comp = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq)
   # APAAC_1  APAAC_2  APAAC_3  APAAC_4  APAAC_5 ...
   # 6.624    3.076    5.757    3.312    6.546 ...

Calculate **Shannon Entropy**:

.. code-block:: python

   se = protpy.shannon_entropy(protein_seq)
   # ShannonEntropy
   # 4.163

Calculate **GRAVY (Grand Average of Hydropathy)**:

.. code-block:: python

   gravy = protpy.gravy(protein_seq)
   # GRAVY
   # -0.045

   #positive value = overall hydrophobic; negative value = overall hydrophilic

Calculate **Aromaticity**:

.. code-block:: python

   aromaticity = protpy.aromaticity(protein_seq)
   # Aromaticity
   # 0.118

Calculate **Instability Index**:

.. code-block:: python

   instability = protpy.instability_index(protein_seq)
   # InstabilityIndex
   # 31.836

   #score < 40 = stable protein; score >= 40 = unstable

Calculate **Isoelectric Point**:

.. code-block:: python

   pi = protpy.isoelectric_point(protein_seq)
   # IsoelectricPoint
   # 5.412

Calculate **Molecular Weight**:

.. code-block:: python

   mw = protpy.molecular_weight(protein_seq)
   # MolecularWeight
   # 139122.355

Calculate **Charge Distribution**:

.. code-block:: python

   #using default pH: ph=7.4
   charge = protpy.charge_distribution(protein_seq)
   # PositiveCharge  NegativeCharge  NetCharge
   # 99.526          114.956         -15.43

   #using custom pH
   charge = protpy.charge_distribution(protein_seq, ph=6.0)

Calculate **Hydrophobic/Polar/Charged Composition (HPC)**:

.. code-block:: python

   hpc = protpy.hydrophobic_polar_charged_composition(protein_seq)
   # Hydrophobic  Polar   Charged
   # 44.542       32.669  18.247

Calculate **Secondary Structure Propensity (SSP)**:

.. code-block:: python

   ssp = protpy.secondary_structure_propensity(protein_seq)
   # Helix  Sheet  Coil
   # 0.983  1.05   1.043

Calculate **k-mer Composition**:

.. code-block:: python

   #using default k=2 (dipeptide frequencies by Chou-Fasman class)
   kmer = protpy.kmer_composition(protein_seq)
   # AA     AC     AD  ...
   # 0.797  0.159  ... ...

   #using custom k
   kmer = protpy.kmer_composition(protein_seq, k=3)

Calculate **Reduced Alphabet Composition**:

.. code-block:: python

   #using default alphabet_size=6
   reduced = protpy.reduced_alphabet_composition(protein_seq)
   # Group_1  Group_2  Group_3  Group_4  Group_5  Group_6
   # 25.339   34.741   9.163    9.084    10.837   10.837

   #supported sizes: 2, 3, 4, 6
   reduced = protpy.reduced_alphabet_composition(protein_seq, alphabet_size=4)

Calculate **Motif Composition**:

.. code-block:: python

   #using default built-in motifs (N-glycosylation, RGD, KDEL, CxxC, CAAX, PKA, dileucine, PEST)
   motif = protpy.motif_composition(protein_seq)
   # NxST_glycosylation  RGD_integrin  KDEL_retention  CxxC_zinc_finger  ...
   # 23                  0             0               2                 ...

   #using custom motif list
   motif = protpy.motif_composition(protein_seq, motifs=[r'RGD', r'NxS'])

Calculate **Amino Acid Pair Composition**:

.. code-block:: python

   pair = protpy.amino_acid_pair_composition(protein_seq)
   # 400-column DataFrame with class-annotated column names
   # AA_Hydrophobic-Hydrophobic  AA_Hydrophobic-Polar  ...
   # 0.797                       0.159                ...

Calculate **Aliphatic Index**:

.. code-block:: python

   ai = protpy.aliphatic_index(protein_seq)
   # AliphaticIndex
   # 82.725

Calculate **Extinction Coefficient**:

.. code-block:: python

   ec = protpy.extinction_coefficient(protein_seq)
   # ExtCoeff_Reduced  ExtCoeff_Oxidized
   # 140960            143335

Calculate **Boman Index**:

.. code-block:: python

   bi = protpy.boman_index(protein_seq)
   # BomanIndex
   # 0.119

Calculate **Aggregation Propensity**:

.. code-block:: python

   ap = protpy.aggregation_propensity(protein_seq)
   # AggregProneRegions  AggregProneFraction
   # 58                  11.793

Calculate **Hydrophobic Moment**:

.. code-block:: python

   #using default parameters: window=11, angle=100 (alpha-helix)
   hm = protpy.hydrophobic_moment(protein_seq)
   # HydrophobicMoment_Mean  HydrophobicMoment_Max
   # 0.272                   0.813

Autocorrelation Descriptors
----------------------------

Calculate **Moreaubroto Autocorrelation (MBAuto)**:

.. code-block:: python

   #using default parameters: lag=30, properties=[...], normalize=True
   moreaubroto_auto = protpy.moreaubroto_autocorrelation(protein_seq)
   # MBAuto_CIDH920105_1  MBAuto_CIDH920105_2  MBAuto_CIDH920105_3 ...
   # -0.052               -0.104               -0.156 ...

Calculate **Moran Autocorrelation (MAuto)**:

.. code-block:: python

   moran_auto = protpy.moran_autocorrelation(protein_seq)
   # MAuto_CIDH920105_1  MAuto_CIDH920105_2  MAuto_CIDH920105_3 ...
   # -0.07786            -0.07879            -0.07906 ...

Calculate **Geary Autocorrelation (GAuto)**:

.. code-block:: python

   geary_auto = protpy.geary_autocorrelation(protein_seq)
   # GAuto_CIDH920105_1  GAuto_CIDH920105_2  GAuto_CIDH920105_3 ...
   # 1.057               1.077               1.04 ...

Conjoint Triad Descriptor
--------------------------

Calculate **Conjoint Triad (CTriad)**:

.. code-block:: python

   conjoint_triad = protpy.conjoint_triad(protein_seq)
   # 111  112  113  114  115 ...
   # 7    17   11   3    6 ...

CTD Descriptors
---------------

Calculate **CTD — Composition, Transition, Distribution**:

.. code-block:: python

   #using default parameters: property="hydrophobicity", all_ctd=True
   ctd = protpy.ctd(protein_seq)
   # hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  hydrophobicity_CTD_C_03 ...
   # 0.279                    0.386                    0.335 ...

Sequence Order Descriptors
---------------------------

Calculate a single **Sequence Order Coupling Number (SOCN)**:

.. code-block:: python

   #using default parameters: d=1, distance_matrix="schneider-wrede"
   socn = protpy.sequence_order_coupling_number_(protein_seq)
   # 401.387

Calculate all SOCNs across a lag:

.. code-block:: python

   #using default parameters: lag=30, distance_matrix="schneider-wrede"
   socn_all = protpy.sequence_order_coupling_number(protein_seq)
   # SOCN_SW1  SOCN_SW2  SOCN_SW3  SOCN_SW4  SOCN_SW5 ...
   # 401.387   409.243   376.946   393.042   396.196 ...

   #using custom parameters: lag=10, distance_matrix="grantham"
   socn_all = protpy.sequence_order_coupling_number(protein_seq, lag=10, distance_matrix="grantham")
   # SOCN_Grant1  SOCN_Grant2  SOCN_Grant3 ...
   # 399.125      402.153      387.820 ...

Calculate **Quasi Sequence Order (QSO)**:

.. code-block:: python

   #using default parameters: lag=30, weight=0.1, distance_matrix="schneider-wrede"
   qso = protpy.quasi_sequence_order(protein_seq)
   # QSO_SW1   QSO_SW2   QSO_SW3   QSO_SW4   QSO_SW5 ...
   # 0.005692  0.002643  0.004947  0.002846  0.005625 ...

   #using custom parameters: lag=10, weight=0.2, distance_matrix="grantham"
   qso = protpy.quasi_sequence_order(protein_seq, lag=10, weight=0.2, distance_matrix="grantham")
   # QSO_Grant1  QSO_Grant2  QSO_Grant3 ...
   # 0.123287    0.079967    0.04332 ...

.. note::
   A demo of the software is available `here <https://colab.research.google.com/drive/12E3ayovpZOf6Gv-8ILwkN1zReKslEksB?usp=sharing>`_.