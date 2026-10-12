Welcome to protpy's documentation 🧬!
======================================

.. image:: https://raw.githubusercontent.com/amckenna41/protpy/main/images/protpy.png
   :alt: protpy
   :width: 600px
   :align: center

|

**protpy** is a Python package for generating a variety of physicochemical, biochemical and structural
descriptors for proteins. All descriptors are calculated using sequence-derived or physicochemical
features of the constituent amino acids. These descriptors have been widely studied and applied across
Bioinformatics applications including protein engineering, SAR (sequence-activity-relationships),
predicting protein structure and function, subcellular localization, protein-protein interactions,
and drug-target interactions.

The descriptors available in **protpy** are:

Composition Descriptors (22)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

Autocorrelation Descriptors (3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Moreaubroto Autocorrelation (MBAuto)**
* **Moran Autocorrelation (MAuto)**
* **Geary Autocorrelation (GAuto)**

Conjoint Triad (1)
~~~~~~~~~~~~~~~~~~~
* **Conjoint Triad (CTriad)**

CTD Descriptors (4)
~~~~~~~~~~~~~~~~~~~~~
* **CTD Composition**
* **CTD Transition**
* **CTD Distribution**
* **CTD Combined**

Sequence Order Descriptors (5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Sequence Order Coupling Number — single (SOCN)**
* **Sequence Order Coupling Number — series**
* **Sequence Order Coupling Number — all matrices**
* **Quasi Sequence Order (QSO)**
* **Quasi Sequence Order — all matrices**

v1.5.0 Composition Descriptors (7)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **cksaap** (CKSAAP)
* **cksaagp** (CKSAAGP)
* **dde** (DDE)
* **asdc** (ASDC)
* **distance_pair** (DP)
* **natural_vector** (NV)
* **atomic_composition** (ATC)

v1.5.0 Scales Descriptors (11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **z_scales_5** (Z5)
* **vhse** (VHSE)
* **atchley_factors** (AF)
* **t_scales** (T)
* **st_scales** (ST)
* **protfp** (ProtFP)
* **fasgai** (FASGAI)
* **blosum_indices** (BLOSUM)
* **ms_whim** (MSWHIM)
* **auto_cross_covariance** (ACC)
* **physicochemical_distance_transformation** (PDT)

v1.5.0 Physicochemical Descriptors (11)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **charge_fractions** (FCR)
* **ncpr** (NCPR)
* **kappa** (Kappa)
* **sequence_charge_decoration** (SCD)
* **charge_density** (ChargeDensity)
* **hydrophobicity_profile** (HP)
* **hydrophobic_ratio** (HydrophobicRatio)
* **half_life** (HalfLife)
* **ivywrel_index** (IVYWREL)
* **cvp_bias** (CvP)
* **flexibility** (Flex)

v1.5.0 Complexity Descriptors (3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **wootton_federhen_complexity** (WF)
* **linguistic_complexity** (LC)
* **lempel_ziv_complexity** (LZ)

v1.5.0 Signal Descriptors (2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **fft_descriptors** (FFT)
* **dct_descriptors** (DCT)

v1.5.0 Structural Descriptors (10)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **aggrescan** (AGG)
* **uversky_charge_hydropathy** (UV)
* **fold_index** (FI)
* **disorder_propensity** (DP)
* **transmembrane_features** (TM)
* **signal_peptide_features** (SP)
* **functional_site_counts** (FS)
* **plaac** (PLAAC)
* **antigenicity_surface** (Antigenicity)
* **cysteine_features** (Cys)

The v1.5.0 descriptors are added to ``all_descriptors`` with ``descriptor_set="extended"`` or ``"all"``;
see :doc:`usage` and :doc:`descriptors`.

.. note::

   A demo of the software is available `on Colab <https://colab.research.google.com/drive/12E3ayovpZOf6Gv-8ILwkN1zReKslEksB?usp=sharing>`_
   and a Medium article covering the background of ``protpy`` is available on `Medium <https://medium.com/@ajmckenna69/protpy-05f02f821baa>`_.

License
=======
**protpy** is distributed under the MIT license.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   api
   descriptors
   contributing