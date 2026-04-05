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

.. note::

   A demo of the software is available `on Colab <https://colab.research.google.com/drive/12E3ayovpZOf6Gv-8ILwkN1zReKslEksB?usp=sharing>`_
   and a Medium article covering the background of ``protpy`` is available `on Medium <https://medium.com/@ajmckenna69/protpy-05f02f821baa>`_.

License
=======
**protpy** is distributed under the MIT license.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   api
   descriptors
   changelog
   contributing