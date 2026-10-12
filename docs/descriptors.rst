Descriptor Reference
====================

A complete reference for all protein descriptors available in the **protpy** package,
grouped by category.

.. contents:: Contents
   :local:
   :depth: 2

----

Setup
-----

.. code-block:: console

   pip install protpy

.. code-block:: python

   import protpy

   # Load a sequence from a FASTA file
   from Bio import SeqIO
   with open("protein.fasta") as f:
       protein_seq = str(next(SeqIO.parse(f, "fasta")).seq)

----

Composition Descriptors
-----------------------

Amino Acid Composition (AAComp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fraction of each of the 20 standard amino acids in the sequence.

.. code-block:: python

   result = protpy.amino_acid_composition(protein_seq)
   # Shape: 1 x 20
   # A      C      D      E      F  ...
   # 6.693  3.108  5.817  3.347  6.614 ...

----

Dipeptide Composition (DPComp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequency of all 400 possible dipeptide (two-residue) combinations.

.. code-block:: python

   result = protpy.dipeptide_composition(protein_seq)
   # Shape: 1 x 400
   # AA    AC    AD    AE    AF  ...
   # 0.72  0.16  0.48  0.4   0.24 ...

----

Tripeptide Composition (TPComp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequency of all 8000 possible tripeptide (three-residue) combinations.

.. code-block:: python

   result = protpy.tripeptide_composition(protein_seq)
   # Shape: 1 x 8000
   # AAA  AAC  AAD  AAE  AAF ...
   # 1    0    0    2    0 ...

----

Grand Average of Hydropathy (GRAVY)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the Kyte-Doolittle hydropathy values across all residues. A positive value
indicates overall hydrophobicity; negative indicates overall hydrophilicity.

.. code-block:: python

   result = protpy.gravy(protein_seq)
   # Shape: 1 x 1
   # GRAVY
   # -0.045

----

Aromaticity
~~~~~~~~~~~

Fraction of aromatic residues (F, W, Y, H) in the sequence.

.. code-block:: python

   result = protpy.aromaticity(protein_seq)
   # Shape: 1 x 1
   # Aromaticity
   # 0.118

----

Instability Index
~~~~~~~~~~~~~~~~~

Stability classifier based on dipeptide instability weight values (DIWV). Values
below 40 indicate a stable protein; 40 or above indicates instability.

.. code-block:: python

   result = protpy.instability_index(protein_seq)
   # Shape: 1 x 1
   # InstabilityIndex
   # 31.836

----

Isoelectric Point
~~~~~~~~~~~~~~~~~

Estimated pH at which the protein carries no net charge, calculated iteratively
using standard pKa values for ionisable residues.

.. code-block:: python

   result = protpy.isoelectric_point(protein_seq)
   # Shape: 1 x 1
   # IsoelectricPoint
   # 5.412

----

Molecular Weight
~~~~~~~~~~~~~~~~

Average molecular weight of the protein calculated from residue masses, corrected
for water lost at each peptide bond.

.. code-block:: python

   result = protpy.molecular_weight(protein_seq)
   # Shape: 1 x 1
   # MolecularWeight (Da)
   # 139122.355

----

Charge Distribution
~~~~~~~~~~~~~~~~~~~~

Positive, negative, and net charge contributions of ionisable residues at a given
pH using the Henderson-Hasselbalch equation.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``ph``
     - float
     - ``7.4``
     - pH at which to calculate charge

.. code-block:: python

   # Default pH 7.4
   result = protpy.charge_distribution(protein_seq)

   # Custom pH
   result = protpy.charge_distribution(protein_seq, ph=6.0)

   # Shape: 1 x 3
   # PositiveCharge  NegativeCharge  NetCharge
   # 99.526          114.956         -15.43

----

Hydrophobic/Polar/Charged Composition (HPC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Percentage of residues belonging to each of three physicochemical groups:
hydrophobic (A, C, F, I, L, M, V, W, Y), polar (G, N, Q, S, T), and charged
(D, E, H, K, R).

.. code-block:: python

   result = protpy.hydrophobic_polar_charged_composition(protein_seq)
   # Shape: 1 x 3
   # Hydrophobic  Polar   Charged
   # 44.542       32.669  18.247

----

Secondary Structure Propensity (SSP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Average Chou-Fasman propensity values for alpha-helix, beta-sheet, and random coil
conformations across all residues.

.. code-block:: python

   result = protpy.secondary_structure_propensity(protein_seq)
   # Shape: 1 x 3
   # Helix  Sheet  Coil
   # 0.983  1.05   1.043

----

k-mer Composition
~~~~~~~~~~~~~~~~~

Frequency of all possible k-length residue subsequences, expressed as a percentage
of total k-mers.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``k``
     - int
     - ``2``
     - Length of each k-mer

.. code-block:: python

   # Default k=2 (dipeptide frequencies)
   result = protpy.kmer_composition(protein_seq)

   # Custom k
   result = protpy.kmer_composition(protein_seq, k=3)

   # Shape: 1 x 20^k  (e.g. 1 x 400 for k=2)
   # AA     AC     AD  ...
   # 0.797  0.159  ... ...

----

Reduced Alphabet Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Amino acid composition after mapping residues to a reduced alphabet of
physicochemical groups. Supported alphabet sizes: ``2``, ``3``, ``4``, ``6``.

.. list-table::
   :header-rows: 1
   :widths: 20 10 15 55

   * - Parameter
     - Type
     - Default
     - Description
   * - ``alphabet_size``
     - int
     - ``6``
     - Number of reduced groups

.. code-block:: python

   # Default alphabet_size=6
   result = protpy.reduced_alphabet_composition(protein_seq)

   # Custom size
   result = protpy.reduced_alphabet_composition(protein_seq, alphabet_size=4)

   # Shape: 1 x alphabet_size
   # Group_1  Group_2  Group_3  Group_4  Group_5  Group_6
   # 25.339   34.741   9.163    9.084    10.837   10.837

----

Motif Composition
~~~~~~~~~~~~~~~~~

Count of occurrences (including overlapping) of biological sequence motifs matched
via regular expressions. Eight built-in motifs are used by default; a custom list
can be supplied.

**Default motifs:**

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Column
     - Pattern
     - Biological meaning
   * - ``NxST_glycosylation``
     - ``N[^P][ST]``
     - N-linked glycosylation site
   * - ``RGD_integrin``
     - ``RGD``
     - Integrin-binding RGD motif
   * - ``KDEL_retention``
     - ``KDEL``
     - ER retention signal
   * - ``CxxC_zinc_finger``
     - ``C..C``
     - Zinc-finger CxxC motif
   * - ``CAAX_prenylation``
     - ``C[A-Z]{2}[CSIM]$``
     - CAAX prenylation box
   * - ``cAMP_PKA``
     - ``[RK]{2}.[ST]``
     - cAMP/PKA phosphorylation site
   * - ``dileucine_sorting``
     - ``[DE]xxxL[LI]``
     - Dileucine lysosomal sorting signal
   * - ``PEST_region``
     - ``P.{1,10}[ED]``
     - PEST degradation signal

.. list-table::
   :header-rows: 1
   :widths: 20 10 15 55

   * - Parameter
     - Type
     - Default
     - Description
   * - ``motifs``
     - list or None
     - ``None``
     - Custom regex patterns; uses built-in 8 if ``None``

.. code-block:: python

   # Default built-in motifs
   result = protpy.motif_composition(protein_seq)

   # Custom motif list
   result = protpy.motif_composition(protein_seq, motifs=[r'RGD', r'N[^P][ST]'])

   # Shape: 1 x len(motifs)
   # NxST_glycosylation  RGD_integrin  KDEL_retention  CxxC_zinc_finger  ...
   # 23                  0             0               2                 ...

----

Amino Acid Pair Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequency of all 400 residue-pair combinations with column names annotated by the
physicochemical class of each residue (Hydrophobic, Polar, Charged, or Other).

.. code-block:: python

   result = protpy.amino_acid_pair_composition(protein_seq)
   # Shape: 1 x 400
   # AA_Hydrophobic-Hydrophobic  AA_Hydrophobic-Polar  AA_Hydrophobic-Charged  ...
   # 0.797                       0.159                 ...                     ...

----

Aliphatic Index
~~~~~~~~~~~~~~~

A measure of the relative volume occupied by aliphatic side chains (Ala, Val, Ile,
Leu). Higher values indicate greater thermostability. Formula:
AI = Ala% + 2.9×Val% + 3.9×(Ile%+Leu%).

.. code-block:: python

   result = protpy.aliphatic_index(protein_seq)
   # Shape: 1 x 1
   # AliphaticIndex
   # 82.725

----

Extinction Coefficient
~~~~~~~~~~~~~~~~~~~~~~

The molar extinction coefficient at 280 nm, calculated from the number of Trp (W),
Tyr (Y), and Cys (C) residues. Reported for both reduced (no disulfide bonds) and
oxidized (all Cys paired) states.

.. code-block:: python

   result = protpy.extinction_coefficient(protein_seq)
   # Shape: 1 x 2
   # ExtCoeff_Reduced  ExtCoeff_Oxidized
   # 140960            143335

----

Boman Index
~~~~~~~~~~~

Sum of solubility values for amino acids divided by sequence length, predicting
potential for protein–protein interactions. Positive values suggest membrane-binding
or interaction potential.

.. code-block:: python

   result = protpy.boman_index(protein_seq)
   # Shape: 1 x 1
   # BomanIndex
   # 0.119

----

Aggregation Propensity
~~~~~~~~~~~~~~~~~~~~~~

Estimates aggregation-prone regions via a sliding-window approach combining
Kyte–Doolittle hydrophobicity and charge neutrality. Returns the count of qualifying
windows and the fraction of the sequence covered.

.. list-table::
   :header-rows: 1
   :widths: 25 10 15 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``window``
     - int
     - ``5``
     - Sliding window size
   * - ``hydrophobicity_threshold``
     - float
     - ``2.0``
     - Minimum mean hydrophobicity
   * - ``charge_threshold``
     - int
     - ``0``
     - Maximum charged residues per window

.. code-block:: python

   result = protpy.aggregation_propensity(protein_seq)
   # Shape: 1 x 2
   # AggregProneRegions  AggregProneFraction
   # 58                  11.793

----

Hydrophobic Moment
~~~~~~~~~~~~~~~~~~

The mean and maximum hydrophobic moment across sliding windows, using the Eisenberg
hydrophobicity scale and a helical-wheel projection. Captures amphipathicity of
putative helix segments.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``window``
     - int
     - ``11``
     - Sliding window size
   * - ``angle``
     - int
     - ``100``
     - Residue rotation angle in degrees (100° for α-helix)

.. code-block:: python

   result = protpy.hydrophobic_moment(protein_seq)
   # Shape: 1 x 2
   # HydrophobicMoment_Mean  HydrophobicMoment_Max
   # 0.272                   0.813

----

Shannon Entropy
~~~~~~~~~~~~~~~

Information-theoretic measure of amino acid diversity in a sequence. Computed as:

.. math::

   H = -\sum_i p_i \log_2 p_i

where :math:`p_i` is the fractional frequency of each amino acid type present. A
value of 0 indicates a completely repetitive (single amino acid) sequence; the
theoretical maximum of :math:`\log_2(20) \approx 4.322` bits corresponds to a
perfectly uniform distribution across all 20 canonical amino acids. Widely used as a
low-complexity filter and diversity measure in ML feature pipelines.

.. code-block:: python

   result = protpy.shannon_entropy(protein_seq)
   # Shape: 1 x 1
   # ShannonEntropy
   # 4.163

----

Pseudo Amino Acid Composition (PAAComp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Augmented amino acid composition that incorporates sequence-order effects via
correlation factors derived from physicochemical properties. Reduces the
dimensionality problem of pure sequence-order information while retaining more
sequence information than simple AAComp.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lamda``
     - int
     - ``30``
     - Number of sequence-order correlation factors to include
   * - ``weight``
     - float
     - ``0.05``
     - Weighting factor for correlation layers
   * - ``properties``
     - list
     - ``[]``
     - AAIndex accession numbers to use (uses built-in set if empty)

.. code-block:: python

   # Default parameters
   result = protpy.pseudo_amino_acid_composition(protein_seq)

   # Custom parameters
   result = protpy.pseudo_amino_acid_composition(protein_seq, lamda=10, weight=0.1)

   # Shape: 1 x (20 + lamda)  →  1 x 50 with defaults
   # PAAC_1  PAAC_2  PAAC_3  ...
   # 0.127   0.059   0.111   ...

----

Amphiphilic Pseudo Amino Acid Composition (APAAComp)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extension of PAAComp that uses both hydrophobicity and hydrophilicity properties to
capture amphiphilic patterns (dual hydrophobic/hydrophilic character) along the
sequence.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lamda``
     - int
     - ``30``
     - Number of sequence-order correlation factors
   * - ``weight``
     - float
     - ``0.5``
     - Weighting factor for correlation layers
   * - ``properties``
     - list
     - ``[]``
     - AAIndex accession numbers (defaults to hydrophobicity + hydrophilicity)

.. code-block:: python

   # Default parameters
   result = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq)

   # Custom parameters
   result = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq, lamda=15, weight=0.3)

   # Shape: 1 x (20 + 2*lamda)  →  1 x 80 with defaults
   # APAAC_1  APAAC_2  APAAC_3  ...
   # 6.624    3.076    5.757    ...

----

Autocorrelation Descriptors
----------------------------

Autocorrelation descriptors measure the correlation between physicochemical property
values of residues separated by a lag distance along the sequence. By default, 8
AAIndex properties are used, generating ``lag × 8 = 240`` features.

**Default properties:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - AAIndex ID
     - Property
   * - ``CIDH920105``
     - Normalised average hydrophobicity
   * - ``BHAR880101``
     - Average flexibility indices
   * - ``CHAM820101``
     - Polarizability parameter
   * - ``CHAM820102``
     - Free energy of solution in water
   * - ``CHOC760101``
     - Residue accessible surface area in tripeptide
   * - ``BIGC670101``
     - Residue volume
   * - ``CHAM810101``
     - Steric parameter
   * - ``DAYM780201``
     - Relative mutability

----

Moreaubroto Autocorrelation (MBAuto)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Uses raw property values as the basis for correlation measurement.

.. list-table::
   :header-rows: 1
   :widths: 15 10 25 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum lag distance
   * - ``properties``
     - list
     - (8 defaults above)
     - AAIndex accession numbers
   * - ``normalize``
     - bool
     - ``True``
     - Normalise output values

.. code-block:: python

   # Default parameters
   result = protpy.moreaubroto_autocorrelation(protein_seq)

   # Custom parameters
   result = protpy.moreaubroto_autocorrelation(protein_seq, lag=15, properties=["CIDH920105"])

   # Shape: 1 x (lag × len(properties))  →  1 x 240 with defaults
   # MBAuto_CIDH920105_1  MBAuto_CIDH920105_2  ...
   # -0.052               -0.104               ...

----

Moran Autocorrelation (MAuto)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Uses the deviation from the mean property value, making it mean-centred and thereby
less sensitive to the absolute property scale.

.. list-table::
   :header-rows: 1
   :widths: 15 10 25 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum lag distance
   * - ``properties``
     - list
     - (8 defaults above)
     - AAIndex accession numbers
   * - ``normalize``
     - bool
     - ``True``
     - Normalise output values

.. code-block:: python

   # Default parameters
   result = protpy.moran_autocorrelation(protein_seq)

   # Custom parameters
   result = protpy.moran_autocorrelation(protein_seq, lag=15)

   # Shape: 1 x (lag × len(properties))  →  1 x 240 with defaults
   # MAuto_CIDH920105_1  MAuto_CIDH920105_2  ...
   # -0.07786            -0.07879            ...

----

Geary Autocorrelation (GAuto)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Uses squared differences between property values at each lag, making it sensitive to
local dissimilarities rather than global correlation.

.. list-table::
   :header-rows: 1
   :widths: 15 10 25 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum lag distance
   * - ``properties``
     - list
     - (8 defaults above)
     - AAIndex accession numbers
   * - ``normalize``
     - bool
     - ``True``
     - Normalise output values

.. code-block:: python

   # Default parameters
   result = protpy.geary_autocorrelation(protein_seq)

   # Custom parameters
   result = protpy.geary_autocorrelation(protein_seq, lag=10, normalize=False)

   # Shape: 1 x (lag × len(properties))  →  1 x 240 with defaults
   # GAuto_CIDH920105_1  GAuto_CIDH920105_2  ...
   # 1.057               1.077               ...

----

Conjoint Triad Descriptor
--------------------------

Conjoint Triad (CTriad)
~~~~~~~~~~~~~~~~~~~~~~~~

Encodes the sequence using a 7-class reduced amino acid alphabet and computes the
frequency of all consecutive three-residue combinations (triads). The 7 classes are:
(1) AGV, (2) ILFP, (3) YMTS, (4) HNQW, (5) RK, (6) DE, (7) C.

.. code-block:: python

   result = protpy.conjoint_triad(protein_seq)
   # Shape: 1 x 343  (7 × 7 × 7 class combinations)
   # 111  112  113  114  ...
   # 7    17   11   3    ...

----

CTD Descriptors
----------------

CTD (Composition, Transition, Distribution) descriptors characterise the distribution
of residues belonging to three physicochemical classes along the sequence. Seven
physicochemical properties are supported.

**Supported properties:**

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Property key
     - Description
     - Classes
   * - ``hydrophobicity``
     - Hydrophobicity
     - Polar / Neutral / Hydrophobic
   * - ``normalized_vdwv``
     - Normalised van der Waals volume
     - 0–2.78 / 2.95–4.0 / 4.03–8.08
   * - ``polarity``
     - Polarity
     - 4.9–6.2 / 8.0–9.2 / 10.4–13.0
   * - ``charge``
     - Charge
     - Positive / Neutral / Negative
   * - ``secondary_struct``
     - Secondary structure
     - Helix / Strand / Coil
   * - ``solvent_accessibility``
     - Solvent accessibility
     - Buried / Exposed / Intermediate
   * - ``polarizability``
     - Polarizability
     - 0–0.108 / 0.128–0.186 / 0.219–0.409

----

CTD Composition
~~~~~~~~~~~~~~~

Fraction of residues in each of the three physicochemical classes.

.. list-table::
   :header-rows: 1
   :widths: 15 10 25 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``property``
     - str
     - ``"hydrophobicity"``
     - Physicochemical property to use

.. code-block:: python

   result = protpy.ctd_composition(protein_seq)
   result = protpy.ctd_composition(protein_seq, property="charge")
   # hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  hydrophobicity_CTD_C_03
   # 0.279                    0.386                    0.335

----

CTD Transition
~~~~~~~~~~~~~~

Fraction of transitions between each pair of the three physicochemical classes.

.. list-table::
   :header-rows: 1
   :widths: 15 10 25 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``property``
     - str
     - ``"hydrophobicity"``
     - Physicochemical property to use

.. code-block:: python

   result = protpy.ctd_transition(protein_seq)
   result = protpy.ctd_transition(protein_seq, property="polarity")
   # hydrophobicity_CTD_T_12  hydrophobicity_CTD_T_13  hydrophobicity_CTD_T_23
   # 0.181                    0.161                    0.179

----

CTD Distribution
~~~~~~~~~~~~~~~~

Position of the first, 25%, 50%, 75%, and last residue of each class within the
sequence (as a percentage of total length).

.. list-table::
   :header-rows: 1
   :widths: 15 10 25 50

   * - Parameter
     - Type
     - Default
     - Description
   * - ``property``
     - str
     - ``"hydrophobicity"``
     - Physicochemical property to use

.. code-block:: python

   result = protpy.ctd_distribution(protein_seq)
   result = protpy.ctd_distribution(protein_seq, property="secondary_struct")
   # hydrophobicity_CTD_D_01_01  hydrophobicity_CTD_D_02_01  ...
   # 0.0796                      0.557                       ...

----

CTD Combined (``ctd_``)
~~~~~~~~~~~~~~~~~~~~~~~

Calculate Composition, Transition **and** Distribution for one or all supported
properties in a single call.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``property``
     - str
     - ``"hydrophobicity"``
     - Property to use when ``all_ctd=False``
   * - ``all_ctd``
     - bool
     - ``True``
     - If ``True``, compute CTD for all 7 supported properties

.. code-block:: python

   # All 7 properties (default)
   result = protpy.ctd_(protein_seq)

   # Single property
   result = protpy.ctd_(protein_seq, property="charge", all_ctd=False)

   # Shape: 1 x (3 + 3 + 15) per property  →  1 x 147 for all 7 properties
   # hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  ...
   # 0.279                    0.386                    ...

----

Sequence Order Descriptors
---------------------------

Sequence order descriptors capture the effect of residue spacing along the sequence
using physicochemical distance matrices. Two distance matrices are supported:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Matrix
     - Description
   * - ``schneider-wrede``
     - Physicochemical distance based on Schneider-Wrede scale (default)
   * - ``grantham``
     - Physicochemical distance based on Grantham's amino acid difference formula

----

Sequence Order Coupling Number — single (``sequence_order_coupling_number_``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computes the sum of squared physicochemical distances between all residue pairs
separated by a gap of ``d``. Returns a single float.

.. list-table::
   :header-rows: 1
   :widths: 20 10 25 45

   * - Parameter
     - Type
     - Default
     - Description
   * - ``d``
     - int
     - ``1``
     - Gap between residue pairs
   * - ``distance_matrix``
     - str
     - ``"schneider-wrede"``
     - Distance matrix to use

.. code-block:: python

   result = protpy.sequence_order_coupling_number_(protein_seq)
   result = protpy.sequence_order_coupling_number_(protein_seq, d=5, distance_matrix="grantham")
   # Returns: 401.387  (float)

----

Sequence Order Coupling Number — series (``sequence_order_coupling_number``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates SOCN values across all gaps from 1 to ``lag``.

.. list-table::
   :header-rows: 1
   :widths: 20 10 25 45

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum gap value
   * - ``distance_matrix``
     - str
     - ``"schneider-wrede"``
     - Distance matrix to use

.. code-block:: python

   # Default parameters
   result = protpy.sequence_order_coupling_number(protein_seq)

   # Custom lag and matrix
   result = protpy.sequence_order_coupling_number(protein_seq, lag=10, distance_matrix="grantham")

   # Shape: 1 x lag  →  1 x 30 with defaults
   # SOCN_SW1   SOCN_SW2   SOCN_SW3  ...
   # 401.387    409.243    376.946   ...

----

Sequence Order Coupling Number — all matrices (``sequence_order_coupling_number_all``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates SOCN across all lags using **both** the Schneider-Wrede and Grantham
matrices and concatenates the results.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum gap value

.. code-block:: python

   result = protpy.sequence_order_coupling_number_all(protein_seq)
   result = protpy.sequence_order_coupling_number_all(protein_seq, lag=15)
   # Shape: 1 x (2 × lag)  →  1 x 60 with defaults
   # SOCN_SW1  ...  SOCN_Grant1  ...

----

Quasi Sequence Order (``quasi_sequence_order``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extends SOCN by combining standard amino acid composition with sequence-order
coupling numbers, weighted by a factor ``weight``. Captures both residue type and
spatial distribution information.

.. list-table::
   :header-rows: 1
   :widths: 20 10 25 45

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum lag value
   * - ``weight``
     - float
     - ``0.1``
     - Weighting factor for coupling terms
   * - ``distance_matrix``
     - str
     - ``"schneider-wrede"``
     - Distance matrix to use

.. code-block:: python

   # Default parameters
   result = protpy.quasi_sequence_order(protein_seq)

   # Custom parameters
   result = protpy.quasi_sequence_order(protein_seq, lag=10, weight=0.2, distance_matrix="grantham")

   # Shape: 1 x (20 + lag)  →  1 x 50 with defaults
   # QSO_SW1    QSO_SW2    QSO_SW3   ...
   # 0.005692   0.002643   0.004947  ...

----

Quasi Sequence Order — all matrices (``quasi_sequence_order_all``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates Quasi Sequence Order using **both** distance matrices and concatenates
the results.

.. list-table::
   :header-rows: 1
   :widths: 15 10 15 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``lag``
     - int
     - ``30``
     - Maximum lag value
   * - ``weight``
     - float
     - ``0.1``
     - Weighting factor for coupling terms

.. code-block:: python

   result = protpy.quasi_sequence_order_all(protein_seq)
   result = protpy.quasi_sequence_order_all(protein_seq, lag=15, weight=0.05)
   # Shape: 1 x (2 × (20 + lag))  →  1 x 100 with defaults
   # QSO_SW1  ...  QSO_Grant1  ...

----

v1.5.0 Descriptors
------------------

Descriptors added in v1.5.0, each with its formula, output shape and primary reference.
``descriptor_set="extended"`` in ``all_descriptors`` / ``encode_dataset`` adds the descriptors marked
*extended*; ``"all"`` adds the rest as well. See ``planning/BLOCKED.md`` for plan items that
could not be implemented without unverifiable data.

cksaap (CKSAAP) — Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Composition of k-spaced amino acid pairs: CKSAAP(XY,g) = N_g(XY) / (L-g-1) x 100 for g = 0..gap.

.. code-block:: text

   Calculate Composition of k-Spaced Amino Acid Pairs (CKSAAP) of protein
   sequence. CKSAAP describes the frequency of each amino acid pair separated
   by g residues within a protein sequence, for each gap g in the range
   0..gap. With 400 possible amino acid pairs per gap this creates a
   400 * (gap+1) Dimensional vector, calculated as:

   CKSAAP(XY, g) = N_g(XY) / (L - g - 1) * 100

   where CKSAAP(XY, g) is the composition of the amino acid pair XY with g
   residues between X and Y, N_g(XY) is the number of such pairs in the sequence
   and L is the length of the sequence. A gap of 0 is the same as the
   dipeptide composition. A gap with no pairs (L - g - 1 <= 0) gives values
   of 0.

.. code-block:: python

   result = protpy.cksaap(protein_seq)
   # Shape: 1 x 400·(gap+1)  (default 1 x 2400)

* **Parameters:** ``gap``: int (default=5) — maximum spacing
* **Descriptor set:** all
* **Reference:** Chen, Kurgan & Ruan (2007). BMC Struct. Biol., 7, 25. https://doi.org/10.1186/1472-6807-7-25

----

cksaagp (CKSAAGP) — Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Composition of k-spaced amino acid group pairs over the iFeature 5-group alphabet.

.. code-block:: text

   Calculate Composition of k-Spaced Amino Acid Group Pairs (CKSAAGP) of protein
   sequence. CKSAAGP is the same as CKSAAP but the amino acids are first mapped
   to 5 groups: aliphatic (GAVLMI), aromatic (FYW), positive charge (KRH),
   negative charge (DE) and uncharged (STCPNQ). With 25 possible group pairs
   per gap this creates a 25 * (gap+1) Dimensional vector, calculated as:

   CKSAAGP(G1G2, g) = N_g(G1G2) / (L - g - 1) * 100

   where CKSAAGP(G1G2, g) is the composition of the group pair G1G2 with g
   residues between them, N_g(G1G2) is the number of such pairs in the sequence
   and L is the length of the sequence. The 25 values for each gap sum to 100.

.. code-block:: python

   result = protpy.cksaagp(protein_seq)
   # Shape: 1 x 25·(gap+1)  (default 1 x 150)

* **Parameters:** ``gap``: int (default=5) — maximum spacing
* **Descriptor set:** extended, all
* **Reference:** Chen et al. (2018). Bioinformatics, 34(14), 2499-2502 (iFeature). https://doi.org/10.1093/bioinformatics/bty140

----

dde (DDE) — Composition
~~~~~~~~~~~~~~~~~~~~~~~

Dipeptide deviation from expected mean: (Dc - Tm) / sqrt(Tv), with Tm from codon counts.

.. code-block:: text

   Calculate Dipeptide Deviation from Expected Mean (DDE) of protein sequence.
   DDE compares the observed dipeptide composition of a protein sequence with
   the composition expected from the number of codons that code for each amino
   acid. For each of the 400 dipeptides it is calculated as:

   Dc(XY) = N(XY) / (L - 1)

   Tm(XY) = (C_X / 61) * (C_Y / 61)

   Tv(XY) = Tm(XY) * (1 - Tm(XY)) / (L - 1)

   DDE(XY) = (Dc(XY) - Tm(XY)) / sqrt(Tv(XY))

   where Dc(XY) is the observed composition of dipeptide XY, N(XY) is the
   number of XY dipeptides, L is the length of the sequence, Tm(XY) and Tv(XY)
   are the theoretical mean and variance, and C_X and C_Y are the number of
   codons for amino acids X and Y out of the 61 sense codons. Unlike the
   other composition descriptors, DDE values are not percentages.

.. code-block:: python

   result = protpy.dde(protein_seq)
   # Shape: 1 x 400

* **Parameters:** sequence only
* **Descriptor set:** all
* **Reference:** Saravanan & Gautham (2015). OMICS, 19(10), 648-658. https://doi.org/10.1089/omi.2015.0095

----

asdc (ASDC) — Composition
~~~~~~~~~~~~~~~~~~~~~~~~~

Adaptive skip dipeptide composition: composition of all ordered pairs (s_i, s_j), i < j.

.. code-block:: text

   Calculate Adaptive Skip Dipeptide Composition (ASDC) of protein sequence.
   ASDC is a modified dipeptide composition that counts every ordered pair of
   amino acids (s_i, s_j), where i < j, at all gaps along the sequence, rather
   than just adjacent pairs, calculated as:

   ASDC(XY) = sum_g N_g(XY) / sum_g sum_XY N_g(XY) * 100

   where ASDC(XY) is the composition of the amino acid pair XY, N_g(XY) is the
   number of XY pairs separated by g residues, and the denominator is the total
   number of ordered pairs in the sequence, L(L-1)/2, where L is the length of
   the sequence. The 400 values sum to 100.

.. code-block:: python

   result = protpy.asdc(protein_seq)
   # Shape: 1 x 400

* **Parameters:** sequence only
* **Descriptor set:** all
* **Reference:** Wei et al. (2018). Bioinformatics, 34(23), 4007-4016 (ACPred-FL). https://doi.org/10.1093/bioinformatics/bty451

----

distance_pair (DP) — Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PseAAC of distance-pairs over a reduced cluster profile (cp13/cp14/cp19/cp20).

.. code-block:: text

   Calculate Distance Pair (DP) descriptor of protein sequence. The descriptor
   incorporates amino acid distance pairs and a reduced amino acid alphabet into
   the general pseudo amino acid composition, as used in iDNA-Prot|dis.
   The amino acids are first mapped to the clusters of a cluster profile cp(n),
   and the descriptor is the composition of the clusters followed by the
   composition of cluster pairs at each distance d = 1..rng:

   DP(c) = N(c) / L * 100

   DP(c1c2, d) = N_d(c1c2) / (L - d) * 100

   where N(c) is the number of residues in cluster c, N_d(c1c2) is the number of
   residue pairs (s_i, s_i+d) whose residues are in clusters c1 and c2, and L is
   the length of the sequence. A distance of 1 is the nearest-neighbour pair.
   The available cluster profiles are cp13, cp14, cp19 and cp20 (no reduction).
   The published cp(14) lists E twice; the extra E in the first cluster is
   dropped, as it is the only reading that gives 14 clusters.

.. code-block:: python

   result = protpy.distance_pair(protein_seq)
   # Shape: 1 x (G + G²·rng)  (default 1 x 602)

* **Parameters:** ``rng``: int (default=3) — maximum pair distance; ``alphabet``: str (default='cp14') — cluster profile
* **Descriptor set:** all
* **Reference:** Liu et al. (2014). PLoS ONE, 9(9), e106691 (iDNA-Prot|dis). https://doi.org/10.1371/journal.pone.0106691

----

natural_vector (NV) — Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Natural vector: count, mean position and normalised second central moment of each amino acid.

.. code-block:: text

   Calculate Natural Vector (NV) of protein sequence. The natural vector
   describes the number of each amino acid in a protein sequence together with
   the mean and the normalized central moment of their positions, calculated as:

   n_k = number of occurrences of amino acid k

   mu_k = sum of the positions of k / n_k

   D2_k = sum of (position of k - mu_k)^2 / (n_k * L)

   where positions start at 1, L is the length of the sequence and k is each of
   the 20 canonical amino acids. mu_k and D2_k are 0 when amino acid k is
   not in the sequence.

.. code-block:: python

   result = protpy.natural_vector(protein_seq)
   # Shape: 1 x 60

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Yu et al. (2013). J. Theor. Biol., 318, 197-204. https://doi.org/10.1016/j.jtbi.2012.11.005

----

atomic_composition (ATC) — Composition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Percentage of C, H, N, O and S atoms over the residues' free amino-acid formulas.

.. code-block:: text

   Calculate Atomic Composition (ATC) of protein sequence. ATC describes the
   percentage of carbon, hydrogen, nitrogen, oxygen and sulphur atoms in a
   protein sequence, using the molecular formulas of the free amino acids,
   calculated as:

   ATC(e) = N_e / (N_C + N_H + N_N + N_O + N_S) * 100

   where ATC(e) is the atomic composition of element e and N_e is the number of
   atoms of element e summed over all residues. By default the formulas of
   the free amino acids are summed without removing water, matching Pfeature
   . With subtract_water=True, one water molecule is removed for each of the
   L - 1 peptide bonds, giving the formula of the peptide (e.g. "GG" -> C4H8N2O3).

.. code-block:: python

   result = protpy.atomic_composition(protein_seq)
   # Shape: 1 x 5

* **Parameters:** ``subtract_water``: bool (default=False) — remove (L-1) H2O
* **Descriptor set:** extended, all
* **Reference:** Pande et al. (2023). J. Comput. Biol., 30(2), 204-222 (Pfeature). https://doi.org/10.1089/cmb.2022.0241

----

z_scales_5 (Z5) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the Sandberg 5 z-scales over the sequence (optionally per domain).

.. code-block:: text

   Calculate Sandberg 5 z-scales (Z5) of protein sequence. The z-scales are the
   principal components of a set of physicochemical properties measured for 87
   coded and non-coded amino acids. z1 mainly describes lipophilicity, z2
   steric bulk and polarizability, z3 polarity and charge, and z4 and z5
   electronegativity, heat of formation, electrophilicity and hardness.

   Each amino acid in the sequence is replaced by its 5 z-scales and the mean of
   each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.z_scales_5(protein_seq)
   # Shape: 1 x 5·domains  (default 1 x 5)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Sandberg et al. (1998). J. Med. Chem., 41(14), 2481-2491. https://doi.org/10.1021/jm9700575

----

vhse (VHSE) — Scales
~~~~~~~~~~~~~~~~~~~~

Mean of the VHSE (hydrophobic, steric, electronic) scales over the sequence (optionally per domain).

.. code-block:: text

   Calculate Vectors of Hydrophobic, Steric and Electronic properties (VHSE) of
   protein sequence. The VHSE scales are principal components score vectors derived
   from physicochemical variables of the 20 coded amino acids. VHSE1-2 describe
   hydrophobic properties, VHSE3-4 steric properties and VHSE5-8 electronic
   properties.

   Each amino acid in the sequence is replaced by its 8 VHSE scales and the mean of
   each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.vhse(protein_seq)
   # Shape: 1 x 8·domains  (default 1 x 8)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Mei et al. (2005). Biopolymers, 80(6), 775-786. https://doi.org/10.1002/bip.20296

----

atchley_factors (AF) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the Atchley factors I-V over the sequence (optionally per domain).

.. code-block:: text

   Calculate Atchley Factors (AF) of protein sequence. The Atchley factors are 5
   factors obtained from a factor analysis of 54 amino acid attributes. Factor
   I mainly reflects polarity, accessibility and hydrophobicity, factor II
   propensity for secondary structure, factor III molecular size, factor IV
   relative amino acid composition and factor V electrostatic charge.

   Each amino acid in the sequence is replaced by its 5 Atchley factors and the
   mean of each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.atchley_factors(protein_seq)
   # Shape: 1 x 5·domains  (default 1 x 5)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Atchley et al. (2005). PNAS, 102(18), 6395-6400. https://doi.org/10.1073/pnas.0408677102

----

t_scales (T) — Scales
~~~~~~~~~~~~~~~~~~~~~

Mean of the topological T-scales over the sequence (optionally per domain).

.. code-block:: text

   Calculate Topological Scales (T-scales) of protein sequence. The T-scales are
   principal components of topological descriptors of amino acids, describing their
   structural and topological characteristics.

   Each amino acid in the sequence is replaced by its 5 T-scales and the mean of
   each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.t_scales(protein_seq)
   # Shape: 1 x 5·domains  (default 1 x 5)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Tian, Zhou & Li (2007). J. Mol. Struct., 830, 106-115. https://doi.org/10.1016/j.molstruc.2006.07.004

----

st_scales (ST) — Scales
~~~~~~~~~~~~~~~~~~~~~~~

Mean of the structural topology ST-scales over the sequence (optionally per domain).

.. code-block:: text

   Calculate Structural Topology Scales (ST-scales) of protein sequence. The ST-
   scales are principal components of 3D structural and topological properties of
   amino acids.

   Each amino acid in the sequence is replaced by its 8 ST-scales and the mean of
   each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.st_scales(protein_seq)
   # Shape: 1 x 8·domains  (default 1 x 8)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Yang et al. (2010). Amino Acids, 38(3), 805-816. https://doi.org/10.1007/s00726-009-0287-y

----

protfp (ProtFP) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the ProtFP (PCA, 8 components) over the sequence (optionally per domain).

.. code-block:: text

   Calculate ProtFP (ProtFP) of protein sequence. ProtFP is a set of 8 principal
   components calculated from a selection of physicochemical properties in the
   AAindex database.

   Each amino acid in the sequence is replaced by its 8 ProtFP components and the
   mean of each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.protfp(protein_seq)
   # Shape: 1 x 8·domains  (default 1 x 8)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** van Westen et al. (2013). J. Cheminform., 5, 41. https://doi.org/10.1186/1758-2946-5-41

----

fasgai (FASGAI) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the FASGAI factor scales over the sequence (optionally per domain).

.. code-block:: text

   Calculate Factor Analysis Scales of Generalized Amino Acid Information (FASGAI)
   of protein sequence. The FASGAI scales are 6 factors obtained from a factor
   analysis of amino acid properties, describing hydrophobicity, alpha and turn
   propensities, bulky properties, compositional characteristics, local flexibility
   and electronic properties.

   Each amino acid in the sequence is replaced by its 6 FASGAI factors and the mean
   of each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.fasgai(protein_seq)
   # Shape: 1 x 6·domains  (default 1 x 6)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Liang & Li (2007). QSAR Comb. Sci., 26(6), 754-763. https://doi.org/10.1002/qsar.200630145

----

blosum_indices (BLOSUM) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the BLOSUM62-derived indices over the sequence (optionally per domain).

.. code-block:: text

   Calculate BLOSUM Indices (BLOSUM) of protein sequence. The BLOSUM indices are 10
   factors derived from the BLOSUM62 substitution matrix by principal component
   analysis and VARIMAX rotation, so they describe the evolutionary similarity of
   amino acids.

   Each amino acid in the sequence is replaced by its 10 BLOSUM indices and the
   mean of each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.blosum_indices(protein_seq)
   # Shape: 1 x 10·domains  (default 1 x 10)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Georgiev (2009). J. Comput. Biol., 16(5), 703-723. https://doi.org/10.1089/cmb.2008.0173

----

ms_whim (MSWHIM) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the MS-WHIM scores over the sequence (optionally per domain).

.. code-block:: text

   Calculate Molecular Surface Weighted Holistic Invariant Molecular Scores (MS-
   WHIM) of protein sequence. The MS-WHIM scores are 3 principal components of 3D
   descriptors calculated from the molecular electrostatic potential on the surface
   of each amino acid.

   Each amino acid in the sequence is replaced by its 3 MS-WHIM scores and the mean
   of each is taken over the sequence, calculated as:

   S_j = (1/L) * sum_i z_j(s_i)

   where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
   j for the amino acid at position i and L is the length of the sequence. With
   domains > 1 the sequence is split into that many near-equal segments and the
   means of each segment are concatenated.

.. code-block:: python

   result = protpy.ms_whim(protein_seq)
   # Shape: 1 x 3·domains  (default 1 x 3)

* **Parameters:** ``domains``: int (default=1) — number of near-equal segments
* **Descriptor set:** extended, all
* **Reference:** Zaliani & Gancia (1999). J. Chem. Inf. Comput. Sci., 39(3), 525-533. https://doi.org/10.1021/ci980211b

----

auto_cross_covariance (ACC) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Auto/cross-covariance of an amino-acid scale: sum z_j(i) z_k(i+l) / (L-l).

.. code-block:: text

   Calculate Auto Cross Covariance (ACC) of protein sequence. ACC transforms the
   amino acid scale values of a sequence into a fixed-length vector that
   describes the covariance between the scale dimensions of residues separated
   by l positions, for each lag l in the range 1..lag. It is calculated as:

   ACC_jk(l) = sum_{i=1}^{L-l} z_j(s_i) * z_k(s_i+l) / (L - l)

   where ACC_jk(l) is the covariance between dimensions j and k at lag l, z_j(s_i)
   is the value of dimension j for the amino acid at position i and L is the length
   of the sequence. When j = k this is the auto covariance and when j != k the cross
   covariance. The scale can be one of the factor-scale families, which are used as
   published, or one or more AAindex accession numbers, whose values are first
   standardised over the 20 amino acids.

.. code-block:: python

   result = protpy.auto_cross_covariance(protein_seq)
   # Shape: 1 x d²·lag  (default 1 x 125)

* **Parameters:** ``scale``: str | list[str] (default='z_scales_5') — factor scale or AAindex code(s); ``lag``: int (default=5) — maximum lag; ``mode``: str (default='acc') — 'acc', 'ac' or 'cc'
* **Descriptor set:** all
* **Reference:** Wold et al. (1993). Anal. Chim. Acta, 277(2), 239-253. https://doi.org/10.1016/0003-2670(93)80437-P

----

physicochemical_distance_transformation (PDT) — Scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Physicochemical distance transformation: mean squared difference of standardised AAindex properties at lag l.

.. code-block:: text

   Calculate Physicochemical Distance Transformation (PDT) of protein sequence.
   PDT describes how much a physicochemical property changes between residues
   separated by l positions, for each lag l in the range 1..lag. For each
   AAindex property it is calculated as:

   PDT_j(l) = sum_{i=1}^{L-l} (P_j(s_i) - P_j(s_i+l))^2 / (L - l)

   where PDT_j(l) is the descriptor value for property j at lag l, P_j(s_i) is the
   value of property j for the amino acid at position i, standardised over the 20
   amino acids, and L is the length of the sequence. By default the same 8 AAindex
   properties as the autocorrelation descriptors are used.

.. code-block:: python

   result = protpy.physicochemical_distance_transformation(protein_seq)
   # Shape: 1 x len(properties)·lag  (default 1 x 24)

* **Parameters:** ``lag``: int (default=3) — maximum lag; ``properties``: list[str] | None — AAindex codes (default: 8 autocorrelation properties)
* **Descriptor set:** extended, all
* **Reference:** Liu et al. (2012). PLoS ONE, 7(9), e46633. https://doi.org/10.1371/journal.pone.0046633

----

charge_fractions (FCR) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fraction of charged residues (FCR), positive/negative fractions and charge asymmetry sigma = NCPR²/FCR.

.. code-block:: text

   Calculate the Charge Fractions (FCR) of protein sequence. The fraction of
   charged residues (FCR), the fractions of positive and negative residues and
   the overall charge asymmetry (sigma) describe the charge content of a protein
   sequence, and are key parameters of the conformations of intrinsically
   disordered proteins. They are calculated as:

   FCR = (n+ + n-) / L

   FPlus = n+ / L, FMinus = n- / L

   Sigma = NCPR^2 / FCR

   where n+ is the number of positive residues (K, R), n- the number of negative
   residues (D, E), L is the length of the sequence and NCPR = (n+ - n-) / L. His
   is treated as neutral and sigma is 0 when there are no charged residues.

.. code-block:: python

   result = protpy.charge_fractions(protein_seq)
   # Shape: 1 x 4

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Das & Pappu (2013). PNAS, 110(33), 13392-13397. https://doi.org/10.1073/pnas.1304749110

----

ncpr (NCPR) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Net charge per residue: (n+ - n-) / L with K/R positive and D/E negative.

.. code-block:: text

   Calculate the Net Charge Per Residue (NCPR) of protein sequence. NCPR is the
   net charge of a protein sequence divided by its length, calculated as:

   NCPR = (n+ - n-) / L

   where n+ is the number of positive residues (K, R), n- the number of negative
   residues (D, E) and L is the length of the sequence.

.. code-block:: python

   result = protpy.ncpr(protein_seq)
   # Shape: 1 x 1

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Mao et al. (2010). PNAS, 107(18), 8183-8188. https://doi.org/10.1073/pnas.0911107107

----

kappa (Kappa) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Das-Pappu kappa: normalised blob charge-asymmetry deviation measuring segregation of opposite charges (0 mixed, 1 segregated).

.. code-block:: text

   Calculate the Kappa (Kappa) charge patterning parameter of protein sequence.
   Kappa describes how well oppositely charged residues are mixed or segregated
   along a protein sequence, which strongly affects the conformations of
   intrinsically disordered proteins. The sequence is split into overlapping
   blobs of g residues and kappa is calculated as:

   sigma_b = (f+ - f-)^2 / (f+ + f-)

   delta_g = mean over all blobs of (sigma_b - sigma)^2

   kappa_g = delta_g / delta_max,g

   where sigma_b is the charge asymmetry of a blob, f+ and f- are its fractions of
   positive (K, R) and negative (D, E) residues, sigma = NCPR^2 / FCR is the charge
   asymmetry of the whole sequence and delta_max,g is the largest delta over the
   most segregated arrangements of the same composition. The final kappa is the
   mean of kappa_g over the blob sizes (5 and 6 by default). Kappa ranges from 0
   (well mixed) to 1 (fully segregated). It is reported as 0 when the sequence has
   no charged residues or only one type of charge, and blobs longer than the
   sequence are reduced to its length.

.. code-block:: python

   result = protpy.kappa(protein_seq)
   # Shape: 1 x (1 + len(blob_sizes))  (default 1 x 3)

* **Parameters:** ``blob_sizes``: tuple[int] (default=(5, 6)) — blob sizes averaged
* **Descriptor set:** extended, all
* **Reference:** Das & Pappu (2013). PNAS, 110(33), 13392-13397. https://doi.org/10.1073/pnas.1304749110

----

sequence_charge_decoration (SCD) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequence charge decoration: (1/L) sum_{m>n} q_m q_n sqrt(m-n).

.. code-block:: text

   Calculate the Sequence Charge Decoration (SCD) of protein sequence. SCD
   describes the patterning of charges along a protein sequence, weighting each
   pair of charges by the square root of their separation. It is calculated as:

   SCD = (1/L) * sum_{m=2}^{L} sum_{n=1}^{m-1} q_m * q_n * (m - n)^(1/2)

   where q_m and q_n are the charges of the residues at positions m and n (+1 for
   K and R, -1 for D and E, 0 otherwise) and L is the length of the sequence. More
   negative values indicate more segregated opposite charges and uncharged
   sequences have an SCD of 0.

.. code-block:: python

   result = protpy.sequence_charge_decoration(protein_seq)
   # Shape: 1 x 1

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Sawle & Ghosh (2015). J. Chem. Phys., 143(8), 085101. https://doi.org/10.1063/1.4929391

----

charge_density (ChargeDensity) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Net charge at the given pH divided by molecular weight (charges per Da).

.. code-block:: text

   Calculate the Charge Density of protein sequence. Charge density is the net
   charge of a protein sequence at a given pH divided by its molecular weight,
   calculated as:

   ChargeDensity = NetCharge(pH) / MW

   where NetCharge(pH) is calculated using the charge_distribution function and MW
   using the molecular_weight function.

.. code-block:: python

   result = protpy.charge_density(protein_seq)
   # Shape: 1 x 1

* **Parameters:** ``ph``: float (default=7.4) — pH
* **Descriptor set:** extended, all
* **Reference:** Müller et al. (2017). Bioinformatics, 33(17), 2753-2755 (modlAMP). https://doi.org/10.1093/bioinformatics/btx285

----

hydrophobicity_profile (HP) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary statistics of a sliding-window hydrophobicity profile (ProtScale-style edge weighting).

.. code-block:: text

   Calculate the Hydrophobicity Profile (HP) of protein sequence. Each amino acid
   in the sequence is replaced by its value on a hydrophobicity scale and the
   values are averaged over a sliding window, giving a profile along the sequence
   . The profile is summarised by its mean, maximum, minimum and standard
   deviation, and optionally the fraction of windows above a threshold, e.g. 1.6
   for transmembrane segments with the Kyte-Doolittle scale. The window value
   is calculated as:

   HP(i) = sum_j w_j * H(s_i+j) / sum_j w_j

   where H is the hydrophobicity scale and w_j are the window weights. With
   edge=1.0 all weights are 1 (a plain mean); with edge < 1 the weights fall
   linearly from 1 at the centre of the window to edge at its ends, as in ExPASy
   ProtScale. Windows longer than the sequence use the whole sequence.

.. code-block:: python

   result = protpy.hydrophobicity_profile(protein_seq)
   # Shape: 1 x 4 (1 x 5 with threshold)

* **Parameters:** ``scale``: str (default='kyte_doolittle') — single-value scale or AAindex code; ``window``: int (default=9) — window length; ``edge``: float (default=1.0) — relative end weight; ``threshold``: float | None (default=None) — adds fraction of windows above threshold
* **Descriptor set:** extended, all
* **Reference:** Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132. https://doi.org/10.1016/0022-2836(82)90515-0

----

hydrophobic_ratio (HydrophobicRatio) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

APD hydrophobic ratio: percentage of L, I, V, F, M, A, C, W residues.

.. code-block:: text

   Calculate the Hydrophobic Ratio of protein sequence. The hydrophobic ratio is
   the percentage of hydrophobic residues in a protein sequence, as defined by
   the Antimicrobial Peptide Database (APD), calculated as:

   HydrophobicRatio = N(LIVFMACW) / L * 100

   where N(LIVFMACW) is the number of L, I, V, F, M, A, C and W residues and L
   is the length of the sequence.

.. code-block:: python

   result = protpy.hydrophobic_ratio(protein_seq)
   # Shape: 1 x 1

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Wang (2020). Protein Sci., 29(1), 8-18 (APD); Wang, Li & Wang (2016). Nucleic Acids Res., 44, D1087 (APD3). https://doi.org/10.1002/pro.3702

----

half_life (HalfLife) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N-end rule half-life (hours) in mammalian, yeast and E. coli systems from the N-terminal residue.

.. code-block:: text

   Calculate the estimated Half Life of protein sequence. The half-life is a
   prediction of the time it takes for half of the amount of a protein in a cell
   to disappear after its synthesis. It is estimated from the N-terminal residue
   using the N-end rule, for 3 model systems: mammalian reticulocytes
   (in vitro), yeast and E. coli (in vivo), using the values tabulated by ExPASy
   ProtParam. Values given as lower bounds in ProtParam (">20 hours" for yeast,
   ">10 hours" for E. coli) are reported as the bound. The E. coli value for an
   N-terminal proline is unknown in ProtParam and returned as NaN. The estimate is
   not applicable to N-terminally modified proteins.

.. code-block:: python

   result = protpy.half_life(protein_seq)
   # Shape: 1 x 3

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Bachmair, Finley & Varshavsky (1986). Science, 234(4773), 179-186. https://doi.org/10.1126/science.3018930

----

ivywrel_index (IVYWREL) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fraction of I, V, Y, W, R, E, L residues (thermophilic adaptation index).

.. code-block:: text

   Calculate the IVYWREL Index of protein sequence. The IVYWREL index is the
   fraction of I, V, Y, W, R, E and L residues in a protein sequence, which
   correlates with the optimal growth temperature of the organism. It is
   calculated as:

   IVYWREL = N(IVYWREL) / L

   where N(IVYWREL) is the number of I, V, Y, W, R, E and L residues and L is the
   length of the sequence.

.. code-block:: python

   result = protpy.ivywrel_index(protein_seq)
   # Shape: 1 x 1

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Zeldovich, Berezovsky & Shakhnovich (2007). PLoS Comput. Biol., 3(1), e5. https://doi.org/10.1371/journal.pcbi.0030005

----

cvp_bias (CvP) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Charged versus polar bias: %(D+E+K+R) - %(N+Q+S+T).

.. code-block:: text

   Calculate the Charged versus Polar Bias (CvP) of protein sequence. The CvP bias
   is the difference between the percentages of charged and polar residues in a
   protein sequence, an indicator of hyperthermostability, calculated as:

   CvP = %(D + E + K + R) - %(N + Q + S + T)

   where %(X) is the percentage of residues of type X in the sequence.

.. code-block:: python

   result = protpy.cvp_bias(protein_seq)
   # Shape: 1 x 1

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Suhre & Claverie (2003). J. Biol. Chem., 278(19), 17198-17202. https://doi.org/10.1074/jbc.M301327200

----

flexibility (Flex) — Physicochemical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary statistics of the Vihinen average-B-value flexibility profile (window 9, end weight 0.25).

.. code-block:: text

   Calculate the Flexibility (Flex) profile of protein sequence. Each amino acid
   in the sequence is replaced by its normalised flexibility parameter (B-value)
   and the values are averaged over a sliding window of 9 residues whose end
   weights fall linearly to 0.25 (weights 0.25, 0.4375, 0.625, 0.8125, 1, ...,
   0.25), as in Vihinen et al.. The profile is summarised by its mean,
   maximum, minimum and standard deviation. The average B-values of Vihinen et
   al. are used; their parameter sets for residues with 0, 1 or 2 rigid
   neighbours are not. Compared with Biopython's ProteinAnalysis.flexibility,
   which uses the same values and weights, each window is centred on its middle
   residue and the last window is included.

.. code-block:: python

   result = protpy.flexibility(protein_seq)
   # Shape: 1 x 4

* **Parameters:** ``window``: int (default=9) — window length
* **Descriptor set:** extended, all
* **Reference:** Vihinen, Torkkila & Riikonen (1994). Proteins, 19(2), 141-149. https://doi.org/10.1002/prot.340190207

----

wootton_federhen_complexity (WF) — Complexity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wootton-Federhen local complexity (K1, K2) over sliding windows and the SEG-trigger low-complexity fraction.

.. code-block:: text

   Calculate the Wootton-Federhen Complexity (WF) of protein sequence. The
   Wootton-Federhen measures describe the local compositional complexity of a
   protein sequence over a sliding window, and are the basis of the SEG program
   for finding low-complexity regions. For each window of W residues
   they are calculated as:

   K1 = (1/W) * log_20(W! / prod_i n_i!)

   K2 = -sum_i (n_i / W) * log2(n_i / W)

   where K1 is the complexity state, K2 is the Shannon entropy of the window in
   bits and n_i is the number of residues of amino acid i in the window. The
   low-complexity fraction is the fraction of residues covered by at least one
   window with K2 <= k2_trigger. This is the trigger stage of SEG only; SEG's
   extension and merging of low-complexity segments are not performed. A window
   longer than the sequence uses the whole sequence.

.. code-block:: python

   result = protpy.wootton_federhen_complexity(protein_seq)
   # Shape: 1 x 4

* **Parameters:** ``window``: int (default=12) — window length; ``k2_trigger``: float (default=2.2) — K2 trigger (bits)
* **Descriptor set:** extended, all
* **Reference:** Wootton & Federhen (1993). Comput. Chem., 17(2), 149-163. https://doi.org/10.1016/0097-8485(93)85006-X

----

linguistic_complexity (LC) — Complexity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linguistic complexity: distinct subwords observed / maximum possible (ratio of sums over word lengths).

.. code-block:: text

   Calculate the Linguistic Complexity (LC) of protein sequence. Linguistic
   complexity describes how repetitive a sequence is, as the number of distinct
   subwords (k-mers) it contains relative to the maximum possible number for a
   sequence of its length. It is calculated as:

   LC = sum_{k=1}^{K} U_k / sum_{k=1}^{K} min(20^k, L - k + 1)

   where U_k is the number of distinct k-mers in the sequence, L is the length of
   the sequence and K is the longest subword length (by default L, i.e. all
   subwords). Repetitive sequences have a low LC and random sequences an LC
   close to 1.

.. code-block:: python

   result = protpy.linguistic_complexity(protein_seq)
   # Shape: 1 x 1

* **Parameters:** ``max_k``: int | None (default=None) — longest word length (default L)
* **Descriptor set:** extended, all
* **Reference:** Troyanskaya et al. (2002). Bioinformatics, 18(5), 679-688. https://doi.org/10.1093/bioinformatics/18.5.679

----

lempel_ziv_complexity (LZ) — Complexity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lempel-Ziv (LZ76) complexity c(n) and its normalised value c(n)·log20(n)/n.

.. code-block:: text

   Calculate the Lempel-Ziv Complexity (LZ) of protein sequence. The Lempel-Ziv
   (LZ76) complexity is the number of new patterns, c(n), found when a sequence is
   read from left to right, each pattern being the shortest substring that cannot
   be copied from the preceding part of the sequence. As c(n) grows with the
   length of the sequence, a normalised value is also calculated as:

   LZ_Normalized = c(n) * log_20(n) / n

   where n is the length of the sequence and 20 is the size of the amino acid
   alphabet. The normalised value is 0 for a single residue.

.. code-block:: python

   result = protpy.lempel_ziv_complexity(protein_seq)
   # Shape: 1 x 2

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Lempel & Ziv (1976). IEEE Trans. Inf. Theory, 22(1), 75-81. https://doi.org/10.1109/TIT.1976.1055501

----

fft_descriptors (FFT) — Signal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mean power of the FFT spectrum of a mean-centred property signal (default EIIP) in equal-width frequency bins.

.. code-block:: text

   Calculate Fast Fourier Transform (FFT) descriptors of protein sequence. Each
   amino acid in the sequence is replaced by its value on a property scale, by
   default the electron-ion interaction potential (EIIP) used in the resonant
   recognition model, and the power spectrum of the mean-centred signal is
   calculated as:

   P(f) = abs(FFT(x)(f))^2 / L

   where x is the property signal, f is the frequency in cycles per residue
   (0 to 0.5) and L is the length of the sequence. To give the same number of
   features for any sequence length, the frequency range 0-0.5 is split into
   n_bins equal-width bins and the mean power in each bin is returned; bins
   containing no frequencies are 0.

.. code-block:: python

   result = protpy.fft_descriptors(protein_seq)
   # Shape: 1 x n_bins  (default 1 x 10)

* **Parameters:** ``scale``: str | dict (default='eiip') — single-value scale, AAindex code or dict; ``n_bins``: int (default=10) — number of frequency bins
* **Descriptor set:** extended, all
* **Reference:** Cosic (1994). IEEE Trans. Biomed. Eng., 41(12), 1101-1114. https://doi.org/10.1109/10.335859

----

dct_descriptors (DCT) — Signal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First n orthonormal DCT-II coefficients of a mean-centred property signal (default Kyte-Doolittle).

.. code-block:: text

   Calculate Discrete Cosine Transform (DCT) descriptors of protein sequence.
   Each amino acid in the sequence is replaced by its value on a property scale,
   by default the Kyte-Doolittle hydropathy, and the first n_coeffs coefficients
   of the orthonormal type-II DCT of the mean-centred signal are returned. The
   DCT compresses a sequence-length signal into a fixed number of coefficients
   , calculated as:

   X_k = s_k * sum_{n=0}^{L-1} x_n * cos(pi * k * (2n + 1) / (2L))

   where X_k is coefficient k, x_n is the property value at position n, L is the
   length of the sequence, and s_0 = sqrt(1/L) and s_k = sqrt(2/L) for k > 0. X_0
   is 0 because the signal is mean-centred.

.. code-block:: python

   result = protpy.dct_descriptors(protein_seq)
   # Shape: 1 x n_coeffs  (default 1 x 10)

* **Parameters:** ``scale``: str | dict (default='kyte_doolittle') — single-value scale, AAindex code or dict; ``n_coeffs``: int (default=10) — number of coefficients (<= L)
* **Descriptor set:** extended, all
* **Reference:** Wang et al. (2017). PLoS ONE, 12(9), e0185587. https://doi.org/10.1371/journal.pone.0185587

----

aggrescan (AGG) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AGGRESCAN aggregation hot-spot descriptors from the a3v propensity profile (hot spots: >=5 residues above -0.02, no Pro).

.. code-block:: text

   Calculate AGGRESCAN (AGG) aggregation descriptors of protein sequence.
   AGGRESCAN predicts aggregation "hot spots", short stretches that promote
   protein aggregation, from experimental aggregation propensities (a3v) of the
   20 amino acids. The a3v values are averaged over a sliding window, giving
   the aggregation profile (a4v), with each window value assigned to its central
   residue and the off-centre N- and C-terminal residues given the first and last
   window values. A hot spot is a run of 5 or more residues with an a4v above the
   hot-spot threshold (HST) of -0.02 and containing no proline. Areas are
   calculated by trapezoidal integration of the profile, with HST as the zero axis.

   The 12 output values are: the mean a3v of the sequence (a3vSA); the number of
   hot spots (nHS) and the number per 100 residues (NnHS); the area of the profile
   above HST (AAT) and per residue (AATr); the total area of the hot spots (THSA)
   and per residue (THSAr); the total area of the profile (TA); the sum of a4v
   (a4vSS) and per 100 residues (Na4vSS); the mean a4v in the hot spots (a4vAHS,
   0 if there are none); and the largest single hot-spot area (HSA_max).

   AGGRESCAN chooses a window of 5, 7, 9 or 11 residues from the sequence length
   and gives the two extreme residues specific values; neither rule is given in
   the published supplementary material, so the window is a parameter and the
   extreme residues are treated like the other off-centre residues.

.. code-block:: python

   result = protpy.aggrescan(protein_seq)
   # Shape: 1 x 12

* **Parameters:** ``window``: int (default=5) — odd window length (AGGRESCAN uses 5, 7, 9, 11)
* **Descriptor set:** extended, all
* **Reference:** Conchillo-Solé et al. (2007). BMC Bioinformatics, 8, 65. https://doi.org/10.1186/1471-2105-8-65

----

uversky_charge_hydropathy (UV) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Uversky mean hydropathy <H>, mean net charge <R> and boundary distance HQ = 2.785<H> - <R> - 1.151.

.. code-block:: text

   Calculate the Uversky Charge-Hydropathy (UV) descriptors of protein sequence.
   Uversky et al. showed that natively unfolded proteins can be separated from
   folded proteins by their mean net charge and mean hydropathy. These are
   calculated as:

   <H> = mean of (KD + 4.5) / 9

   <R> = abs(n+ - n-) / L

   HQ = 2.785 <H> - <R> - 1.151

   where KD is the Kyte-Doolittle hydropathy of each residue, n+ and n- are the
   numbers of positive (K, R) and negative (D, E) residues and L is the length
   of the sequence. HQ is the distance from the boundary <R> = 2.785 <H> - 1.151
   (as restated in); the protein is predicted to be unfolded when HQ < 0.

.. code-block:: python

   result = protpy.uversky_charge_hydropathy(protein_seq)
   # Shape: 1 x 4

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Uversky, Gillespie & Fink (2000). Proteins, 41(3), 415-427. https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7

----

fold_index (FI) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FoldIndex profile summary: windowed 2.785<H> - abs(<R>) - 1.151 and disordered-residue statistics.

.. code-block:: text

   Calculate the FoldIndex (FI) of protein sequence. FoldIndex applies the
   Uversky charge-hydropathy boundary over a sliding window to find disordered
   regions of a protein sequence. For each window it is calculated as:

   FI = 2.785 <H> - abs(<R>) - 1.151

   where <H> is the mean Kyte-Doolittle hydropathy of the window rescaled to 0-1
   and <R> is the mean net charge of the window (K, R positive; D, E negative).
   Residues covered by at least one window with FI < 0 are counted as disordered.
   A window longer than the sequence uses the whole sequence, so FI then equals
   UV_HQ from uversky_charge_hydropathy.

.. code-block:: python

   result = protpy.fold_index(protein_seq)
   # Shape: 1 x 4

* **Parameters:** ``window``: int (default=51) — window length
* **Descriptor set:** extended, all
* **Reference:** Prilusky et al. (2005). Bioinformatics, 21(16), 3435-3438. https://doi.org/10.1093/bioinformatics/bti537

----

disorder_propensity (DP) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fractions of disorder- and order-promoting residues and mean TOP-IDP value.

.. code-block:: text

   Calculate the Disorder Propensity (DP) of protein sequence. Intrinsically
   disordered proteins are enriched in disorder-promoting residues (A, R, G, Q,
   S, P, E, K) and depleted in order-promoting residues (W, C, F, I, Y, V, L, N)
   . The descriptor gives the fraction of each, together with the mean value of
   the TOP-IDP scale, a scale of amino acid propensity for intrinsic disorder,
   calculated as:

   DP_TOP_IDP_mean = sum_i TOP-IDP(s_i) / L

   where TOP-IDP(s_i) is the TOP-IDP value of the amino acid at position i and L
   is the length of the sequence. Higher values indicate a higher propensity for
   disorder.

.. code-block:: python

   result = protpy.disorder_propensity(protein_seq)
   # Shape: 1 x 3

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Campen et al. (2008). Protein Pept. Lett., 15(9), 956-963. https://doi.org/10.2174/092986608785849164

----

transmembrane_features (TM) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kyte-Doolittle transmembrane-segment heuristics (segment count, max hydropathy, TM fraction, positive-inside bias).

.. code-block:: text

   Calculate Transmembrane (TM) segment features of protein sequence. Putative
   transmembrane segments are found from the Kyte-Doolittle hydropathy profile:
   windows of 19 residues with a mean hydropathy above 1.6 indicate a likely
   transmembrane segment, and overlapping windows are merged into segments.
   The positive-inside rule states that the cytoplasmic side of a transmembrane
   segment is enriched in positive residues. The 4 output values are: the
   number of segments; the maximum window hydropathy; the fraction of residues
   inside a segment; and the positive-inside bias, the number of K and R residues
   in the 15 residues before the first segment minus the number in the 15 residues
   after it (0 if there is no segment). This is a hydropathy heuristic and not a
   topology predictor.

.. code-block:: python

   result = protpy.transmembrane_features(protein_seq)
   # Shape: 1 x 4

* **Parameters:** ``window``: int (default=19) — hydropathy window; ``kd_threshold``: float (default=1.6) — segment threshold
* **Descriptor set:** extended, all
* **Reference:** Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132. https://doi.org/10.1016/0022-2836(82)90515-0

----

signal_peptide_features (SP) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N-terminal signal-peptide heuristics: n-region net charge and maximum 8-residue hydropathy.

.. code-block:: text

   Calculate Signal Peptide (SP) features of protein sequence. Signal peptides
   have a positively charged N-terminal n-region followed by a hydrophobic
   h-region. The 2 output values are: the net charge of the first 5 residues
   (K and R +1, D and E -1); and the maximum mean Kyte-Doolittle hydropathy over
   8-residue windows within the first n_term residues. The von Heijne
   cleavage-site weight matrix is not included.

.. code-block:: python

   result = protpy.signal_peptide_features(protein_seq)
   # Shape: 1 x 2

* **Parameters:** ``n_term``: int (default=30) — N-terminal region length
* **Descriptor set:** extended, all
* **Reference:** von Heijne (1983). Eur. J. Biochem., 133(1), 17-21. https://doi.org/10.1111/j.1432-1033.1983.tb07424.x

----

functional_site_counts (FS) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Counts and per-100-residue densities of PROSITE N-glycosylation, ER-retention and PTS1 motifs and tryptic sites.

.. code-block:: text

   Calculate Functional Site (FS) counts of protein sequence. The number of each
   functional sequence motif, and its density per 100 residues, is counted for:
   N-glycosylation sites (PROSITE PS00001, N-{P}-[ST]-{P}), the C-terminal
   endoplasmic reticulum retention signal (PS00014, [KRHQSA]-[DENQ]-E-L at the
   C-terminus), the C-terminal peroxisomal targeting signal 1 (PS00342,
   [STAGCN]-[RKH]-[LIVMAFY] at the C-terminus) and tryptic cleavage sites
   (after K or R, not before P). Overlapping matches are counted. These motifs
   complement the default motifs of motif_composition.

.. code-block:: python

   result = protpy.functional_site_counts(protein_seq)
   # Shape: 1 x 8

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Sigrist et al. (2013). Nucleic Acids Res., 41(D1), D344-D347 (PROSITE). https://doi.org/10.1093/nar/gks1067

----

plaac (PLAAC) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~

PLAAC prion-like domain scores: max window LLR, HMM core score, PrD length and presence.

.. code-block:: text

   Calculate the Prion-Like Amino Acid Composition (PLAAC) descriptors of protein
   sequence. PLAAC identifies prion-like domains, whose composition resembles
   that of yeast prion domains, using the amino acid frequencies of prion domains
   and of the S. cerevisiae proteome. Each residue has a log-likelihood ratio:

   LLR(s_i) = ln(f_prion(s_i) / f_background(s_i))

   The 4 output values are: the maximum sum of LLR over any window of core_length
   residues (PLAAC_LLR); the maximum sum of LLR over a window of core_length
   residues lying entirely inside a prion-like region of the two-state hidden
   Markov model Viterbi parse (PLAAC_COREscore); the length of the prion-like
   region containing that core (PLAAC_PrD_length); and whether the sequence has a
   prion-like region of at least core_length residues (PLAAC_has_PrD). As in
   PLAAC, PLAAC_LLR and PLAAC_COREscore are NaN when they are undefined (a
   sequence shorter than core_length, or no prion-like region of that length).
   The parameters are ported from the PLAAC source code (MIT License), using the
   S. cerevisiae background (alpha = 1).

.. code-block:: python

   result = protpy.plaac(protein_seq)
   # Shape: 1 x 4

* **Parameters:** ``core_length``: int (default=60) — minimal core length
* **Descriptor set:** extended, all
* **Reference:** Lancaster et al. (2014). Bioinformatics, 30(17), 2501-2502. https://doi.org/10.1093/bioinformatics/btu310

----

antigenicity_surface (Antigenicity) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Profile statistics of the Parker hydrophilicity and Janin surface scales over 7-residue windows.

.. code-block:: text

   Calculate the Antigenicity and Surface (Antigenicity) profiles of protein
   sequence. Antigenic sites tend to be hydrophilic and exposed on the protein
   surface. Each amino acid in the sequence is replaced by its value on the
   Parker hydrophilicity scale, derived from HPLC peptide retention times,
   and on the Janin scale of free energy of transfer from the inside to the
   outside of a globular protein. Each is averaged over a sliding window of 7
   residues, as used by Parker et al., and the profile is summarised by its mean,
   maximum, minimum and standard deviation.

.. code-block:: python

   result = protpy.antigenicity_surface(protein_seq)
   # Shape: 1 x 8

* **Parameters:** ``window``: int (default=7) — window length
* **Descriptor set:** extended, all
* **Reference:** Parker, Guo & Hodges (1986). Biochemistry, 25(19), 5425-5432. https://doi.org/10.1021/bi00367a013

----

cysteine_features (Cys) — Structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cysteine count, fraction, parity, maximum disulfides, CxxC motifs and mean spacing.

.. code-block:: text

   Calculate Cysteine (Cys) features of protein sequence. Cysteines form
   disulfide bonds and functional motifs such as the CxxC motif of thiol
   oxidoreductases. The 6 output values are: the number of cysteines; their
   fraction; whether the number is even; the maximum number of disulfide bonds
   (the number of cysteines // 2); the number of overlapping C-x-x-C motifs; and
   the mean gap between consecutive cysteines (0 if there are fewer than 2).

.. code-block:: python

   result = protpy.cysteine_features(protein_seq)
   # Shape: 1 x 6

* **Parameters:** sequence only
* **Descriptor set:** extended, all
* **Reference:** Fomenko & Gladyshev (2003). Biochemistry, 42(38), 11214-11225. https://doi.org/10.1021/bi034459s

----

Descriptor Summary
------------------

Speed ratings reflect typical computation time for a single average-length protein
(~500 residues) on a standard CPU:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Rating
     - Meaning
   * - ✅ Fast
     - < 1 ms — simple residue counting or scalar formula
   * - 🟡 Moderate
     - 1–50 ms — sliding window or O(n²) pass
   * - 🔴 Slow
     - > 50 ms — large feature space, iterative convergence, or many property lookups

.. note::

   Autocorrelation, PseAAC, and APAAComp scale with both sequence length **and**
   ``lag`` — reduce ``lag`` or the number of properties to speed them up. Tripeptide
   composition always produces 8000 columns regardless of sequence length.

.. list-table::
   :header-rows: 1
   :widths: 28 35 12 15 7 8

   * - Descriptor
     - Function
     - Output shape
     - Category
     - Speed
     - Complexity
   * - Amino Acid Composition
     - ``amino_acid_composition(seq)``
     - 1 × 20
     - Composition
     - ✅ Fast
     - O(n)
   * - Dipeptide Composition
     - ``dipeptide_composition(seq)``
     - 1 × 400
     - Composition
     - ✅ Fast
     - O(n)
   * - Tripeptide Composition
     - ``tripeptide_composition(seq)``
     - 1 × 8000
     - Composition
     - 🟡 Moderate
     - O(n)
   * - GRAVY
     - ``gravy(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Aromaticity
     - ``aromaticity(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Instability Index
     - ``instability_index(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Isoelectric Point
     - ``isoelectric_point(seq)``
     - 1 × 1
     - Composition
     - 🟡 Moderate
     - O(n · iter)
   * - Molecular Weight
     - ``molecular_weight(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Charge Distribution
     - ``charge_distribution(seq, ph=7.4)``
     - 1 × 3
     - Composition
     - ✅ Fast
     - O(n)
   * - Hydrophobic/Polar/Charged
     - ``hydrophobic_polar_charged_composition(seq)``
     - 1 × 3
     - Composition
     - ✅ Fast
     - O(n)
   * - Secondary Structure Propensity
     - ``secondary_structure_propensity(seq)``
     - 1 × 3
     - Composition
     - ✅ Fast
     - O(n)
   * - k-mer Composition
     - ``kmer_composition(seq, k=2)``
     - 1 × 20^k
     - Composition
     - 🟡 Moderate
     - O(n · 20^k)
   * - Reduced Alphabet Composition
     - ``reduced_alphabet_composition(seq, alphabet_size=6)``
     - 1 × alphabet_size
     - Composition
     - ✅ Fast
     - O(n)
   * - Motif Composition
     - ``motif_composition(seq, motifs=None)``
     - 1 × len(motifs)
     - Composition
     - 🟡 Moderate
     - O(n · m)
   * - Amino Acid Pair Composition
     - ``amino_acid_pair_composition(seq)``
     - 1 × 400
     - Composition
     - ✅ Fast
     - O(n)
   * - Aliphatic Index
     - ``aliphatic_index(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Extinction Coefficient
     - ``extinction_coefficient(seq)``
     - 1 × 2
     - Composition
     - ✅ Fast
     - O(n)
   * - Boman Index
     - ``boman_index(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Aggregation Propensity
     - ``aggregation_propensity(seq, window=5)``
     - 1 × 2
     - Composition
     - 🟡 Moderate
     - O(n · win)
   * - Hydrophobic Moment
     - ``hydrophobic_moment(seq, window=11, angle=100)``
     - 1 × 2
     - Composition
     - 🟡 Moderate
     - O(n · win)
   * - Shannon Entropy
     - ``shannon_entropy(seq)``
     - 1 × 1
     - Composition
     - ✅ Fast
     - O(n)
   * - Pseudo AAComp
     - ``pseudo_amino_acid_composition(seq, lamda=30, weight=0.05)``
     - 1 × (20 + lamda)
     - Composition
     - 🔴 Slow
     - O(n · lamda · props)
   * - Amphiphilic Pseudo AAComp
     - ``amphiphilic_pseudo_amino_acid_composition(seq, lamda=30, weight=0.5)``
     - 1 × (20 + 2×lamda)
     - Composition
     - 🔴 Slow
     - O(n · lamda · props)
   * - Moreaubroto Autocorrelation
     - ``moreaubroto_autocorrelation(seq, lag=30)``
     - 1 × (lag × props)
     - Autocorrelation
     - 🔴 Slow
     - O(n · lag · props)
   * - Moran Autocorrelation
     - ``moran_autocorrelation(seq, lag=30)``
     - 1 × (lag × props)
     - Autocorrelation
     - 🔴 Slow
     - O(n · lag · props)
   * - Geary Autocorrelation
     - ``geary_autocorrelation(seq, lag=30)``
     - 1 × (lag × props)
     - Autocorrelation
     - 🔴 Slow
     - O(n · lag · props)
   * - Conjoint Triad
     - ``conjoint_triad(seq)``
     - 1 × 343
     - Conjoint Triad
     - ✅ Fast
     - O(n)
   * - CTD Composition
     - ``ctd_composition(seq, property="hydrophobicity")``
     - 1 × 3
     - CTD
     - ✅ Fast
     - O(n)
   * - CTD Transition
     - ``ctd_transition(seq, property="hydrophobicity")``
     - 1 × 3
     - CTD
     - ✅ Fast
     - O(n)
   * - CTD Distribution
     - ``ctd_distribution(seq, property="hydrophobicity")``
     - 1 × 15
     - CTD
     - ✅ Fast
     - O(n)
   * - CTD Combined
     - ``ctd_(seq, property="hydrophobicity", all_ctd=True)``
     - 1 × 147
     - CTD
     - 🟡 Moderate
     - O(n · props)
   * - SOCN (single)
     - ``sequence_order_coupling_number_(seq, d=1)``
     - float
     - Sequence Order
     - ✅ Fast
     - O(n)
   * - SOCN (series)
     - ``sequence_order_coupling_number(seq, lag=30)``
     - 1 × lag
     - Sequence Order
     - 🟡 Moderate
     - O(n · lag)
   * - SOCN (all matrices)
     - ``sequence_order_coupling_number_all(seq, lag=30)``
     - 1 × (2 × lag)
     - Sequence Order
     - 🟡 Moderate
     - O(n · lag)
   * - Quasi Sequence Order
     - ``quasi_sequence_order(seq, lag=30, weight=0.1)``
     - 1 × (20 + lag)
     - Sequence Order
     - 🟡 Moderate
     - O(n · lag)
   * - Quasi Sequence Order (all)
     - ``quasi_sequence_order_all(seq, lag=30, weight=0.1)``
     - 1 × (2 × (20 + lag))
     - Sequence Order
     - 🟡 Moderate
     - O(n · lag)

----

References
----------

The descriptors implemented in **protpy** are based on the following published methods.

Composition
~~~~~~~~~~~

- Amino acid, dipeptide, and tripeptide composition: Nakashima, H., Nishikawa, K., &
  Ooi, T. (1986). The folding type of a protein is relevant to the amino acid
  composition. *Journal of Biochemistry*, 99(1), 153–162.
- GRAVY: Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the
  hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132.
- Aromaticity: Lobry, J. R., & Gautier, C. (1994). Hydrophobicity, expressivity and
  aromaticity are the major trends of amino-acid usage in 999 *Escherichia coli*
  chromosome-encoded genes. *Nucleic Acids Research*, 22(15), 3174–3180.
- Instability index: Guruprasad, K., Reddy, B. V. B., & Pandit, M. W. (1990).
  Correlation between stability of a protein and its dipeptide composition.
  *Protein Engineering*, 4(2), 155–161.
- Isoelectric point: Bjellqvist, B., et al. (1994). The focusing positions of
  polypeptides in immobilized pH gradients can be predicted from their amino acid
  sequences. *Electrophoresis*, 14(1), 1023–1031.
- Molecular weight, isoelectric point, charge: Gasteiger, E., et al. (2005). Protein
  identification and analysis tools on the ExPASy server. In *The Proteomics Protocols
  Handbook*, Humana Press, 571–607.
- Secondary structure propensity: Chou, P. Y., & Fasman, G. D. (1974). Prediction of
  protein conformation. *Biochemistry*, 13(2), 222–245.
- Aliphatic index: Ikai, A. J. (1980). Thermostability and aliphatic index of globular
  proteins. *Journal of Biochemistry*, 88(6), 1895–1898.
- Extinction coefficient: Pace, C. N., et al. (1995). How to measure and predict the
  molar absorption coefficient of a protein. *Protein Science*, 4(11), 2411–2423.
- Boman index: Boman, H. G. (2003). Antibacterial peptides: basic facts and emerging
  concepts. *Journal of Internal Medicine*, 254(3), 197–215.
- Hydrophobic moment: Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1982). The
  helical hydrophobic moment: a measure of the amphiphilicity of a helix. *Nature*,
  299, 371–374.
- Shannon entropy: Shannon, C. E. (1948). A mathematical theory of communication.
  *Bell System Technical Journal*, 27(3), 379–423.
- Pseudo amino acid composition (PseAAC): Chou, K.-C. (2001). Prediction of protein
  cellular attributes using pseudo-amino acid composition. *Proteins: Structure,
  Function, and Bioinformatics*, 43(3), 246–255.
- Amphiphilic PseAAC (APseAAC): Chou, K.-C. (2005). Using amphiphilic pseudo amino
  acid composition to predict enzyme subfamily classes. *Bioinformatics*, 21(1), 10–19.

Autocorrelation
~~~~~~~~~~~~~~~

- Moreau-Broto autocorrelation: Moreau, G., & Broto, P. (1980). The autocorrelation of
  a topological structure: A new molecular descriptor. *Nouveau Journal de Chimie*,
  4, 359–360.
- Moran autocorrelation: Moran, P. A. P. (1950). Notes on continuous stochastic
  phenomena. *Biometrika*, 37(1–2), 17–23.
- Geary autocorrelation: Geary, R. C. (1954). The contiguity ratio and statistical
  mapping. *The Incorporated Statistician*, 5(3), 115–145.
- AAIndex properties: Kawashima, S., & Kanehisa, M. (2000). AAindex: amino acid index
  database. *Nucleic Acids Research*, 28(1), 374.

Conjoint Triad
~~~~~~~~~~~~~~

- Liu, B., et al. (2008). Prediction of protein-protein interactions based on the naive
  Bayes classifier with amino acid composition features. *Biochemical and Biophysical
  Research Communications*, 368(2), 462–468.

CTD
~~~

- Dubchak, I., et al. (1995). Prediction of protein folding class using global
  description of amino acid sequence. *PNAS*, 92(19), 8700–8704.
- Dubchak, I., et al. (1999). Recognition of a protein fold in the context of the SCOP
  classification. *Proteins*, 35(4), 401–407.

Sequence Order
~~~~~~~~~~~~~~

- Grantham, R. (1974). Amino acid difference formula to help explain protein evolution.
  *Science*, 185(4154), 862–864.
- Schneider, G., & Wrede, P. (1994). The rational design of amino acid sequences by
  artificial neural networks and simulated molecular evolution. *Biophysical Journal*,
  66(2), 335–344.
- Chou, K.-C. (2000). Prediction of protein subcellular locations by incorporating
  quasi-sequence-order effect. *Biochemical and Biophysical Research Communications*,
  278(2), 477–483.
