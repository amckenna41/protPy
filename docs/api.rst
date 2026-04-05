API Reference
=============

All public functions are importable directly from the ``protpy`` namespace after installation.

.. contents:: Modules
   :local:
   :depth: 1

Composition
-----------

.. automodule:: protpy.composition
   :members:
      amino_acid_composition,
      dipeptide_composition,
      tripeptide_composition,
      pseudo_amino_acid_composition,
      amphiphilic_pseudo_amino_acid_composition,
      gravy,
      aromaticity,
      instability_index,
      isoelectric_point,
      molecular_weight,
      charge_distribution,
      hydrophobic_polar_charged_composition,
      secondary_structure_propensity,
      kmer_composition,
      reduced_alphabet_composition,
      motif_composition,
      amino_acid_pair_composition,
      aliphatic_index,
      extinction_coefficient,
      boman_index,
      aggregation_propensity,
      hydrophobic_moment,
      shannon_entropy
   :undoc-members:
   :show-inheritance:

Autocorrelation
---------------

.. automodule:: protpy.autocorrelation
   :members:
      moreaubroto_autocorrelation,
      moran_autocorrelation,
      geary_autocorrelation
   :undoc-members:
   :show-inheritance:

Conjoint Triad
--------------

.. automodule:: protpy.conjoint_triad
   :members: conjoint_triad
   :undoc-members:
   :show-inheritance:

CTD
---

.. automodule:: protpy.ctd
   :members:
      ctd,
      ctd_composition,
      ctd_transition,
      ctd_distribution
   :undoc-members:
   :show-inheritance:

Sequence Order
--------------

.. automodule:: protpy.sequence_order
   :members:
      sequence_order_coupling_number_,
      sequence_order_coupling_number,
      sequence_order_coupling_number_all,
      quasi_sequence_order,
      quasi_sequence_order_all
   :undoc-members:
   :show-inheritance:

.. note::
    A demo of the software and API is available `here <https://colab.research.google.com/drive/1btfEx23bgWdkUPiwdwlDqKkmUp1S-_7U?usp=sharing/>`_.
