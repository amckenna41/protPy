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

Loading Sequences from Files
----------------------------

Use ``read_sequences_from_file()`` to parse sequences from a file. The format
is detected automatically from the extension (``.fasta`` / ``.fa``, ``.txt``,
``.csv``, ``.tsv``).

**FASTA** — returns every record as a ``list[str]`` (a trailing stop codon ``*`` is removed),
or with ``id_col`` a DataFrame of the header identifiers and sequences:

.. code-block:: python

   seqs = protpy.read_sequences_from_file("proteins.fasta")
   # ['MKTAYIAKQRQISFVKSHFS...', 'ACDEFGHIKLM...', ...]

   df = protpy.read_sequences_from_file("proteins.fasta", id_col="id")
   #    id                    sequence
   # 0  sp|P59594|SPIKE_SARS  MFIFLLFLTLTSGSDLDRCT...
   X = protpy.encode_dataset(df)   # the identifiers are kept in the first column

**Plain TXT** — one sequence per non-blank line:

.. code-block:: python

   seqs = protpy.read_sequences_from_file("sequences.txt")

**CSV / TSV** — specify which column holds the sequences:

.. code-block:: python

   # select by column name
   seqs = protpy.read_sequences_from_file("data.csv", seq_col="sequence")

   # select by zero-based column index
   seqs = protpy.read_sequences_from_file("data.tsv", seq_col=1)

**Multi-column files (e.g. thermostability datasets)** — pass ``id_col`` to
retain the identifier column and drop everything else. The function returns a
``pd.DataFrame`` with exactly two columns instead of a plain list:

.. code-block:: python

   # thermostability.txt has columns: name, sequence, T50, is_train, m
   df = protpy.read_sequences_from_file(
       "thermostability.txt",
       seq_col="sequence",
       id_col="name",
   )
   # returns a DataFrame with only ['name', 'sequence'] — all other columns dropped
   #        name                           sequence
   # 0  00000000  MTIKEMPQPKTFGELKNLPLLNTD...
   # 1  11111111  KETSPIPQPKTFGPLGNLPLIDKD...
   # ...

For headerless files set ``has_header=False``; columns can then only be
selected by integer index:

.. code-block:: python

   seqs = protpy.read_sequences_from_file("seqs.csv", has_header=False)

All returned sequences are uppercased and blank entries are excluded.

**Encoding all sequences from a file** — use ``encode_dataset()`` to compute
all descriptors for every sequence in one call.  The first column of the
result is always a unique sequence identifier:

.. code-block:: python

   # Named identifiers from a multi-column file (e.g. thermostability.txt)
   dataset = protpy.read_sequences_from_file(
       "thermostability.txt", seq_col="sequence", id_col="name"
   )
   result = protpy.encode_dataset(dataset, seq_col="sequence")
   # result is (n_sequences × 10 758):
   #        name     A      C      D  ...  QSO_Grant50
   # 0  00000000  6.693  3.108  5.817 ...  0.034943
   # 1  11111111  8.057  ...

   # Plain list — identifiers are auto-generated as seq1, seq2, …
   seqs = ["MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD", "ACDEFGHIKLM..."]
   result = protpy.encode_dataset(seqs)
   # result first column is "sequence_id" with values ["seq1", "seq2", …]

Each completed row (with its identifier) is flushed to the output CSV
immediately, so progress is preserved if the process is interrupted.  Pass
``output_file=None`` to suppress file output.

Composition Descriptors
-----------------------

Calculate **Amino Acid Composition**:

.. code-block:: python

   aac = protpy.amino_acid_composition(protein_seq)
   # A      C      D      E      F ...
   # 6.693  3.108  5.817  3.347  6.614 ...
   # shape: 1 x 20

Calculate **Dipeptide Composition**:

.. code-block:: python

   dpc = protpy.dipeptide_composition(protein_seq)
   # AA     AC     AD     AE     AF ...
   # 0.797  0.159  0.478  0.399  0.239 ...
   # shape: 1 x 400

   #dipeptide deviation from expected mean (alias: protpy.dde)
   dde = protpy.dipeptide_composition(protein_seq, normalize="dde")
   # DDE_AA  DDE_AC  DDE_AD  DDE_AE  DDE_AF ...
   # 1.989   -0.424  2.014   1.405   0.185 ...
   # shape: 1 x 400

Calculate **Tripeptide Composition**:

.. code-block:: python

   #alias of kmer_composition(k=3); 8000 columns, not in any all_descriptors set
   tpc = protpy.tripeptide_composition(protein_seq)
   # AAA    AAC  AAD  AAE    AAF ...
   # 0.080  0    0    0.160  0 ...
   # shape: 1 x 8000

Calculate **k-mer / Reduced Alphabet Composition**:

.. code-block:: python

   #using default parameters: k=2
   kmer = protpy.kmer_composition(protein_seq)
   # AA     AC     AD     AE     AF ...
   # 0.797  0.159  0.478  0.399  0.239 ...
   # shape: 1 x 400

   #hydrophobic / polar / charged groups (formerly hydrophobic_polar_charged_composition)
   hpc = protpy.kmer_composition(protein_seq, k=1, alphabet=3)
   # ReducedAlphabet_1  ReducedAlphabet_2  ReducedAlphabet_3
   # 44.542             32.669             22.789
   # shape: 1 x 3

   #grouped dipeptide composition over the iFeature 5-group alphabet (GDPC)
   gdpc = protpy.kmer_composition(protein_seq, k=2, alphabet="ifeature5")
   # GDPC_aliphatic.aliphatic  GDPC_aliphatic.aromatic  GDPC_aliphatic.positive  GDPC_aliphatic.negative  GDPC_aliphatic.uncharged ...
   # 12.839                    4.625                    3.429                    2.951                    12.121 ...
   # shape: 1 x 25

   #a custom grouping: a list of residue strings covering the 20 amino acids
   custom = protpy.kmer_composition(protein_seq, k=1, alphabet=["AGV", "ILFP", "YMTS", "HNQW", "RK", "DE", "C"])
   # RA_custom_AGV  RA_custom_ILFP  RA_custom_YMTS  RA_custom_HNQW  RA_custom_RK ...
   # 20.239         25.259          21.434          12.908          7.888 ...
   # shape: 1 x 7

Calculate **GRAVY**:

.. code-block:: python

   gravy = protpy.gravy(protein_seq)
   # GRAVY
   # -0.045
   # shape: 1 x 1

Calculate **Aromaticity**:

.. code-block:: python

   aromaticity = protpy.aromaticity(protein_seq)
   # Aromaticity
   # 0.118
   # shape: 1 x 1

Calculate **Instability Index**:

.. code-block:: python

   #score < 40 = stable protein; DIWV table as in Biopython / ExPASy
   instability = protpy.instability_index(protein_seq)
   # InstabilityIndex
   # 32.419
   # shape: 1 x 1

Calculate **Isoelectric Point**:

.. code-block:: python

   #using default parameters: pka_scale="bjellqvist" (as ExPASy / Biopython)
   pi = protpy.isoelectric_point(protein_seq)
   # IsoelectricPoint
   # 5.562
   # shape: 1 x 1

Calculate **Molecular Weight**:

.. code-block:: python

   mw = protpy.molecular_weight(protein_seq)
   # MolecularWeight
   # 139125.120
   # shape: 1 x 1

Calculate **Charge Distribution**:

.. code-block:: python

   #using default parameters: ph=7.4, pka_scale="bjellqvist"
   charge = protpy.charge_distribution(protein_seq)
   # PositiveCharge  NegativeCharge  NetCharge
   # 99.398          114.920         -15.522
   # shape: 1 x 3

   #charge per Da (formerly charge_density)
   charge = protpy.charge_distribution(protein_seq, normalize="mass")
   # PositiveCharge_perDa  NegativeCharge_perDa  NetCharge_perDa
   # 0.000714              0.000826              -0.000112
   # shape: 1 x 3

Calculate **Motif Composition**:

.. code-block:: python

   #using the 8 default motifs
   motif = protpy.motif_composition(protein_seq)
   # NxS/T_glycosylation  RGD_integrin  KDEL_retention  CxxC_zinc_finger  CAAX_prenylation ...
   # 23                   0             0               2                 0 ...
   # shape: 1 x 8

   #built-in PROSITE motifs by name, with counts per 100 residues
   motif = protpy.motif_composition(protein_seq, motifs=["N_glycosylation", "PTS1", "tryptic_site"], per_100=True)
   # N_glycosylation  PTS1  tryptic_site  N_glycosylation_per100  PTS1_per100 ...
   # 23               0     93            1.833                   0 ...
   # shape: 1 x 6

Calculate **Aliphatic Index**:

.. code-block:: python

   ai = protpy.aliphatic_index(protein_seq)
   # AliphaticIndex
   # 82.725
   # shape: 1 x 1

Calculate **Extinction Coefficient**:

.. code-block:: python

   ec = protpy.extinction_coefficient(protein_seq)
   # ExtCoeff_Reduced  ExtCoeff_Oxidized
   # 140960            143335
   # shape: 1 x 2

Calculate **Boman Index**:

.. code-block:: python

   bi = protpy.boman_index(protein_seq)
   # BomanIndex
   # 1.198
   # shape: 1 x 1

Calculate **Hydrophobic Moment**:

.. code-block:: python

   #using default parameters: window=11, angle=100 (alpha-helix)
   hm = protpy.hydrophobic_moment(protein_seq)
   # HydrophobicMoment_Mean  HydrophobicMoment_Max
   # 0.272                   0.813
   # shape: 1 x 2

Calculate **Pseudo Amino Acid Composition**:

.. code-block:: python

   #using default parameters: lamda=30, weight=0.05
   paac = protpy.pseudo_amino_acid_composition(protein_seq)
   # PAAC_1  PAAC_2  PAAC_3  PAAC_4  PAAC_5 ...
   # 1.871   0.869   1.626   0.935   1.849 ...
   # shape: 1 x 50

Calculate **Amphiphilic Pseudo Amino Acid Composition**:

.. code-block:: python

   #using default parameters: lamda=30, weight=0.5
   apaac = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq)
   # APAAC_1  APAAC_2  APAAC_3  APAAC_4  APAAC_5 ...
   # 6.049    2.809    5.257    3.025    5.977 ...
   # shape: 1 x 80

Calculate **CKSAAP / CKSAAGP / ASDC / Distance Pair**:

.. code-block:: python

   #using default parameters: gap=5
   cksaap = protpy.cksaap(protein_seq)
   # CKSAAP_AA_gap0  CKSAAP_AC_gap0  CKSAAP_AD_gap0  CKSAAP_AE_gap0  CKSAAP_AF_gap0 ...
   # 0.797           0.159           0.478           0.399           0.239 ...
   # shape: 1 x 2400

   #grouped pairs (alias: protpy.cksaagp)
   cksaagp = protpy.cksaap(protein_seq, alphabet="ifeature5")
   # CKSAAGP_aliphatic.aliphatic_gap0  CKSAAGP_aliphatic.aromatic_gap0  CKSAAGP_aliphatic.positive_gap0  CKSAAGP_aliphatic.negative_gap0  CKSAAGP_aliphatic.uncharged_gap0 ...
   # 12.839                            4.625                            3.429                            2.951                            12.121 ...
   # shape: 1 x 150

   asdc = protpy.asdc(protein_seq)
   # ASDC_AA  ASDC_AC  ASDC_AD  ASDC_AE  ASDC_AF ...
   # 0.443    0.218    0.329    0.252    0.329 ...
   # shape: 1 x 400

   dp = protpy.distance_pair(protein_seq)
   # DP_IMV  DP_L   DP_F   DP_WY  DP_G ...
   # 15.060  7.888  6.614  5.179  6.295 ...
   # shape: 1 x 602

Calculate **Natural Vector / Atomic Composition**:

.. code-block:: python

   nv = protpy.natural_vector(protein_seq, normalize=True)
   # NV_n_A  NV_n_C  NV_n_D  NV_n_E  NV_n_F ...
   # 0.067   0.031   0.058   0.033   0.066 ...
   # shape: 1 x 60

   atc = protpy.atomic_composition(protein_seq)
   # ATC_C   ATC_H   ATC_N  ATC_O   ATC_S
   # 27.011  52.281  6.952  13.501  0.255
   # shape: 1 x 5

Autocorrelation Descriptors
---------------------------

Calculate **Autocorrelation**:

.. code-block:: python

   #Moreau-Broto, Moran or Geary; the three names below are aliases
   auto = protpy.autocorrelation(protein_seq, method="moran", lag=10)
   # MAuto_CIDH920105_1  MAuto_CIDH920105_2  MAuto_CIDH920105_3  MAuto_CIDH920105_4  MAuto_CIDH920105_5 ...
   # -0.058              -0.078              -0.042              -0.021              -0.054 ...
   # shape: 1 x 80

   #using default parameters: lag=30, the 8 default AAindex properties, normalize=True
   mb = protpy.moreaubroto_autocorrelation(protein_seq)
   # MBAuto_CIDH920105_1  MBAuto_CIDH920105_2  MBAuto_CIDH920105_3  MBAuto_CIDH920105_4  MBAuto_CIDH920105_5 ...
   # -0.052               -0.072               -0.036               -0.015               -0.048 ...
   # shape: 1 x 240

   moran = protpy.moran_autocorrelation(protein_seq)
   # MAuto_CIDH920105_1  MAuto_CIDH920105_2  MAuto_CIDH920105_3  MAuto_CIDH920105_4  MAuto_CIDH920105_5 ...
   # -0.058              -0.078              -0.042              -0.021              -0.054 ...
   # shape: 1 x 240

   #Geary is the autocorrelation in the classic all_descriptors set
   geary = protpy.geary_autocorrelation(protein_seq)
   # GAuto_CIDH920105_1  GAuto_CIDH920105_2  GAuto_CIDH920105_3  GAuto_CIDH920105_4  GAuto_CIDH920105_5 ...
   # 1.057               1.077               1.040               1.020               1.052 ...
   # shape: 1 x 240

Conjoint Triad Descriptors
--------------------------

Calculate **Conjoint Triad**:

.. code-block:: python

   #normalised (f - min) / max as in Shen et al. 2007
   ctriad = protpy.conjoint_triad(protein_seq)
   # conj_triad_111  conj_triad_112  conj_triad_113  conj_triad_114  conj_triad_115 ...
   # 0.409           0.773           0.500           0.136           0.273 ...
   # shape: 1 x 343

   #raw triad counts
   ctriad = protpy.conjoint_triad(protein_seq, normalize=False)
   # conj_triad_111  conj_triad_112  conj_triad_113  conj_triad_114  conj_triad_115 ...
   # 9               17              11              3               6 ...
   # shape: 1 x 343

CTD Descriptors
---------------

Calculate **CTD**:

.. code-block:: python

   #using default parameters: all 7 properties, parts="CTD" (ctd_ is an alias)
   ctd = protpy.ctd(protein_seq)
   # CTD_C_01_hydrophobicity  CTD_C_02_hydrophobicity  CTD_C_03_hydrophobicity  CTD_C_01_normalized_vdwv  CTD_C_02_normalized_vdwv ...
   # 0.279                    0.386                    0.335                    0.420                     0.355 ...
   # shape: 1 x 147

   #transition only, for one property (same as ctd_transition)
   ctd_t = protpy.ctd(protein_seq, property="charge", all_ctd=False, parts="T")
   # CTD_T_12_charge  CTD_T_13_charge  CTD_T_23_charge
   # 0.131            0.019            0.147
   # shape: 1 x 3

Sequence Order Descriptors
--------------------------

Calculate **Sequence Order Coupling Number**:

.. code-block:: python

   #using default parameters: lag=30, distance_matrix="schneider-wrede"
   socn = protpy.sequence_order_coupling_number(protein_seq)
   # SOCN_SW1  SOCN_SW2  SOCN_SW3  SOCN_SW4  SOCN_SW5 ...
   # 401.387   409.243   376.946   393.042   396.196 ...
   # shape: 1 x 30

   #both matrices, divided by the number of pairs (removes the length dependence)
   socn = protpy.sequence_order_coupling_number(protein_seq, lag=10, distance_matrix="both", normalize=True)
   # SOCN_SW1  SOCN_SW2  SOCN_SW3  SOCN_SW4  SOCN_SW5 ...
   # 0.320     0.327     0.301     0.314     0.317 ...
   # shape: 1 x 20

Calculate **Quasi Sequence Order**:

.. code-block:: python

   #using default parameters: lag=30, weight=0.1, distance_matrix="schneider-wrede"; the 50 values sum to 1
   qso = protpy.quasi_sequence_order(protein_seq)
   # QSO_SW1   QSO_SW2   QSO_SW3   QSO_SW4   QSO_SW5 ...
   # 0.000057  0.000026  0.000049  0.000028  0.000056 ...
   # shape: 1 x 50

   #both matrices (as in all_descriptors)
   qso = protpy.quasi_sequence_order(protein_seq, distance_matrix="both")
   # QSO_SW1   QSO_SW2   QSO_SW3   QSO_SW4   QSO_SW5 ...
   # 0.000057  0.000026  0.000049  0.000028  0.000056 ...
   # shape: 1 x 100

Physicochemical, Complexity, Spectral and Structural Descriptors
----------------------------------------------------------------

Calculate **Charge**:

.. code-block:: python

   fcr = protpy.charge_fractions(protein_seq)
   # Charge_FCR  Charge_FPlus  Charge_FMinus  Charge_NCPR  Charge_Sigma
   # 0.171       0.079         0.092          -0.013       0.000953
   # shape: 1 x 5

   kappa = protpy.kappa(protein_seq)
   # Kappa
   # 0.153
   # shape: 1 x 1

   scd = protpy.sequence_charge_decoration(protein_seq)
   # Charge_SCD
   # 1.392
   # shape: 1 x 1

Calculate **Hydrophobicity Profile**:

.. code-block:: python

   hp = protpy.hydrophobicity_profile(protein_seq)
   # HP_kyte_doolittle_w9_mean  HP_kyte_doolittle_w9_max  HP_kyte_doolittle_w9_min  HP_kyte_doolittle_w9_std
   # -0.053                     3.189                     -2.822                    0.876
   # shape: 1 x 4

   #hydropathy transmembrane segments (formerly transmembrane_features)
   tm = protpy.hydrophobicity_profile(protein_seq, window=19, threshold=1.6, segments=True)
   # HP_kyte_doolittle_w19_mean  HP_kyte_doolittle_w19_max  HP_kyte_doolittle_w19_min  HP_kyte_doolittle_w19_std  HP_kyte_doolittle_w19_frac_above ...
   # -0.058                      2.811                      -1.911                     0.631                      0.015 ...
   # shape: 1 x 8

   #Vihinen flexibility profile (formerly flexibility)
   flex = protpy.hydrophobicity_profile(protein_seq, scale="vihinen_flexibility", window=9, edge=0.25)
   # HP_vihinen_flexibility_w9_mean  HP_vihinen_flexibility_w9_max  HP_vihinen_flexibility_w9_min  HP_vihinen_flexibility_w9_std
   # 0.995                           1.060                          0.931                          0.019
   # shape: 1 x 4

Calculate **Thermostability Indices**:

.. code-block:: python

   thermo = protpy.thermostability_indices(protein_seq)
   # IVYWREL  CvP_Bias
   # 0.330    -9.323
   # shape: 1 x 2

Calculate **Complexity**:

.. code-block:: python

   h = protpy.shannon_entropy(protein_seq)
   # ShannonEntropy
   # 4.163
   # shape: 1 x 1

   wf = protpy.wootton_federhen_complexity(protein_seq)
   # WF_K2_mean  WF_K2_min  WF_LowComplexity_frac
   # 3.001       2.055      0.030
   # shape: 1 x 3

   lz = protpy.lempel_ziv_complexity(protein_seq)
   # LZ_Normalized
   # 0.860
   # shape: 1 x 1

Calculate **Spectral Descriptors**:

.. code-block:: python

   fft = protpy.spectral_descriptors(protein_seq)
   # FFT_eiip_bin1  FFT_eiip_bin2  FFT_eiip_bin3  FFT_eiip_bin4  FFT_eiip_bin5 ...
   # 0.001361       0.001641       0.001838       0.001610       0.001683 ...
   # shape: 1 x 10

   dct = protpy.spectral_descriptors(protein_seq, transform="dct")
   # DCT_kyte_doolittle_1  DCT_kyte_doolittle_2  DCT_kyte_doolittle_3  DCT_kyte_doolittle_4  DCT_kyte_doolittle_5 ...
   # -1.812                1.051                 2.011                 1.613                 -3.117 ...
   # shape: 1 x 10

Calculate **Factor Scales**:

.. code-block:: python

   #sequence mean of a factor scale
   scales = protpy.aa_scale_descriptors(protein_seq, "atchley_factors")
   # AF_1    AF_2   AF_3    AF_4   AF_5
   # -0.095  0.080  -0.141  0.337  -0.122
   # shape: 1 x 5

   #auto/cross covariance of the 5 z-scales
   acc = protpy.auto_cross_covariance(protein_seq)
   # ACC_1_1_lag1  ACC_1_2_lag1  ACC_1_3_lag1  ACC_1_4_lag1  ACC_1_5_lag1 ...
   # -0.494        -0.181        0.021         -0.101        -0.119 ...
   # shape: 1 x 125

Calculate **Structural**:

.. code-block:: python

   agg = protpy.aggrescan(protein_seq)
   # AGG_a3vSA  AGG_NnHS  AGG_AATr  AGG_THSAr  AGG_a4vAHS
   # 0.072      4.701     0.230     0.155      0.409
   # shape: 1 x 5

   #Uversky charge-hydropathy score of the whole sequence
   ch = protpy.charge_hydropathy(protein_seq)
   # CH_HQ
   # 0.215
   # shape: 1 x 1

   #FoldIndex profile
   fi = protpy.charge_hydropathy(protein_seq, window=51)
   # CH_w51_mean  CH_w51_max  CH_w51_min  CH_w51_std
   # 0.186        0.582       -0.190      0.132
   # shape: 1 x 4

   plaac = protpy.plaac(protein_seq, fill_value=0.0)
   # PLAAC_LLR  PLAAC_COREscore  PLAAC_PrD_length  PLAAC_has_PrD
   # -1.578     0                0                 0
   # shape: 1 x 4

Encoding All Descriptors at Once
---------------------------------

Calculate a **set of descriptors** in a single call using ``all_descriptors()``. With the
default ``descriptor_set="classic"`` this returns a ``DataFrame`` of shape ``1 × 1,402``: amino
acid and dipeptide composition, the ProtParam scalars (GRAVY, aromaticity, instability index,
pI, molecular weight, aliphatic index, extinction coefficient), Boman index, charge
distribution, hydrophobic moment, Shannon entropy, motifs, PseAAC, APseAAC, Geary
autocorrelation, conjoint triad, CTD (7 properties) and quasi-sequence-order with both
distance matrices. (Before v1.6.0 the classic set had 9,957 columns, 8,000 of them tripeptide
composition, plus Moreau-Broto and Moran autocorrelations that duplicate Geary.)

.. code-block:: python

   df = protpy.all_descriptors(protein_seq)              # 1 x 1402
   df = protpy.all_descriptors(protein_seq, decimals=3)  # values are not rounded unless asked

**Progress bar** — pass ``show_progress=True`` to display a live ``tqdm`` progress bar
that names the descriptor currently being computed and tracks overall progress:

.. code-block:: python

   df = protpy.all_descriptors(protein_seq, show_progress=True)
   # [Composition] amino_acid_composition  5%|█                 |  01/20 [00:00<00:00]

Custom per-descriptor parameters can be passed via the ``params`` argument — a dict keyed by
registered descriptor name, whose values are dicts of keyword arguments forwarded to that
descriptor. Descriptors absent from ``params`` use their defaults (the classic set computes
``quasi_sequence_order`` with ``distance_matrix="both"`` and keeps only the
``ExtCoeff_Oxidized`` column of ``extinction_coefficient``).

.. code-block:: python

   df = protpy.all_descriptors(protein_seq, params={
       'charge_distribution':            {'ph': 5.5},
       'pseudo_amino_acid_composition':  {'lamda': 10, 'weight': 0.1},
       'geary_autocorrelation':          {'lag': 15, 'normalize': False},
       'quasi_sequence_order':           {'lag': 20, 'weight': 0.2, 'distance_matrix': 'both'},
   })

Every registered descriptor that takes parameters can be configured this way; see
``protpy.get_descriptor_info(name)["parameters"]``. Passing an unrecognised descriptor name or a
descriptor that takes no extra parameters raises a ``ValueError``; a non-dict ``params`` raises a
``TypeError``.

**Lags** — the autocorrelation, PseAAC/APseAAC and sequence-order descriptors need a lag (or
``lamda``) smaller than the sequence length. Since v1.6.0 calling one of them directly with a
larger lag raises a ``ValueError`` instead of silently replacing it by 30. Inside
``all_descriptors`` / ``encode_dataset`` a sequence that is too short (30 residues or fewer with
the default lags) is computed with lag = length - 1 and the higher-lag columns are NaN, with a
``UserWarning``, so short peptides and long proteins get the same columns; pass smaller lags
through ``params`` to avoid the NaN columns.

**Descriptor sets** — ``descriptor_set`` chooses how much is computed; the sets are nested.
``"extended"`` adds the low-dimensional descriptors (natural vector, atomic composition,
CKSAAGP, charge fractions, kappa, SCD, hydrophobicity profile, thermostability indices,
Wootton-Federhen and Lempel-Ziv complexity, spectral descriptors, AGGRESCAN,
charge-hydropathy, PLAAC); ``"all"`` also adds the high-dimensional ones (CKSAAP, ASDC,
distance pair, auto/cross covariance, Moreau-Broto and Moran autocorrelation, sequence order
coupling numbers). ``include`` and ``exclude`` add or remove individual registered
descriptors by name; tripeptide composition and DDE are in no set and must be included:

.. code-block:: python

   df = protpy.all_descriptors(protein_seq, descriptor_set="extended")
   df = protpy.all_descriptors(protein_seq, descriptor_set="all")
   df = protpy.all_descriptors(protein_seq,
                               include=["kappa", "tripeptide_composition"],
                               exclude=["ctd"],
                               params={"kappa": {"blob_sizes": (5,)}})

   # the same arguments work for whole datasets
   df = protpy.encode_dataset(sequences, descriptor_set="extended", output_file=None)

A few extended descriptors are ``NaN`` where their method leaves them undefined:
``plaac`` LLR/COREscore for sequences shorter than the core length or without a
prion-like domain (pass ``params={"plaac": {"fill_value": 0.0}}`` to replace them), and
``kappa`` for a sequence without charged residues.

**Deprecated functions** — v1.6.0 merged many descriptors into family functions
(``kmer_composition(alphabet=)``, ``cksaap(alphabet=)``, ``autocorrelation(method=)``,
``ctd(parts=)``, ``hydrophobicity_profile(segments=)``, ``spectral_descriptors(transform=)``,
``charge_hydropathy(window=)``, ``thermostability_indices()``, ``aa_scale_descriptors(scale=)``)
and deprecated the rest. Deprecated functions still work, raise a ``DeprecationWarning`` and are
removed in 2.0; :doc:`descriptors` lists their replacements.

Scaling Features for Machine Learning
-------------------------------------

protpy returns each descriptor on its natural scale: percentages, counts, indices, pH
values, Da and so on. Many classifiers expect features on a common scale, for example
[0, 1]. Rather than rescaling inside protpy, fit a scaler on your **training** data and
reuse it for validation, test and new sequences. This scales every column (whatever
descriptors you computed) to the range actually seen in your data, and it is the
standard way to avoid leaking information from the test set:

.. code-block:: python

   from sklearn.preprocessing import MinMaxScaler

   X_train = protpy.encode_dataset(train_seqs, output_file=None, verbose=False).set_index("sequence_id")
   X_test = protpy.encode_dataset(test_seqs, output_file=None, verbose=False).set_index("sequence_id")

   scaler = MinMaxScaler(clip=True).fit(X_train)   # fit on training data only
   X_train_scaled = scaler.transform(X_train)      # every column in [0, 1]
   X_test_scaled = scaler.transform(X_test)        # clip=True keeps unseen extremes in [0, 1]

Notes:

* ``StandardScaler`` (zero mean, unit variance) is a common alternative; tree-based models
  (random forests, gradient boosting) do not need scaling at all.
* Wrap the scaler and model in a scikit-learn ``Pipeline`` when cross-validating, so the
  scaler is re-fitted on each training fold.
* ``descriptor_set="extended"`` can produce a few ``NaN`` values (see above); impute them
  first, e.g. with ``sklearn.impute.SimpleImputer``.
* Scaling with each descriptor's theoretical minimum and maximum is generally a poor
  substitute: real proteins sit far from those extremes, so values end up compressed into
  a narrow band (e.g. instability-index bounds are roughly [-159, 583] while typical values
  are 20-60).

.. note::
   A demo of the software is available `here <https://colab.research.google.com/drive/12E3ayovpZOf6Gv-8ILwkN1zReKslEksB?usp=sharing>`_.