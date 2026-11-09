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

v1.5.0 / v1.6.0 Descriptors (22)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **ctd** (CTD) — CTD
* **natural_vector** (NV) — Composition
* **atomic_composition** (ATC) — Composition
* **cksaagp** (CKSAAGP) — Composition
* **charge_fractions** (FCR) — Physicochemical
* **kappa** (Kappa) — Physicochemical
* **sequence_charge_decoration** (SCD) — Physicochemical
* **hydrophobicity_profile** (HP) — Physicochemical
* **thermostability_indices** (Thermo) — Physicochemical
* **wootton_federhen_complexity** (WF) — Complexity
* **lempel_ziv_complexity** (LZ) — Complexity
* **spectral_descriptors** (Spectral) — Spectral
* **aggrescan** (AGG) — Structural
* **charge_hydropathy** (CH) — Structural
* **plaac** (PLAAC) — Structural
* **cksaap** (CKSAAP) — Composition
* **asdc** (ASDC) — Composition
* **distance_pair** (DP) — Composition
* **auto_cross_covariance** (ACC) — Scales
* **dde** (DDE) — Composition
* **autocorrelation** (Auto) — Autocorrelation
* **aa_scale_descriptors** (Scale) — Scales

v1.6.0 merged most of the v1.5.0 additions into these family functions and deprecated the rest (a
``DeprecationWarning`` is raised; they are removed in 2.0). The classic ``all_descriptors`` set is about
1,400 columns; ``descriptor_set="extended"`` or ``"all"`` add the descriptors above. See :doc:`usage` and
:doc:`descriptors`, including the table of deprecated functions.

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