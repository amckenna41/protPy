# protpy Unit Tests

All tests in the project were ran using Python's unittest testing framework (https://docs.python.org/3/library/unittest.html).

**Run all unittests from main project directory:**
```
python3 -m unittest discover tests
```

You can add the flag *-b* to suppress some of the verbose output when running the unittests.

Unit tests
----------
* `test_protpy.py` - unit tests for overall `protpy` package.
* `test_autocorrellation.py` - unit tests for autocorrelation module and functionality.
* `test_composition.py` - unit tests for composition module and functionality, including pseudo composition and amphiphillic composition.
* `test_conjoint_triad.py` - unit tests for conjoint triad module and functionality.
* `test_ctd.py` - unit tests for CTD module and functionality.
* `test_sequence_order.py` - unit tests for sequence order module and functionality.

Test Files
----------
* `test_data/test_fasta1.fasta` - Spike glycoprotein (P59594 · SPIKE_SARS) [[1]](#references).
* `test_data/test_fasta2.fasta` - Nucleoprotein (P59595 · NCAP_SARS) [[2]](#references).
* `test_data/test_fasta3.fasta` - Envelope small membrane protein (P59637 · VEMP_SARS) [[3]](#references).
* `test_data/test_fasta4.fasta` - Nucleoprotein (P0DTC9 · NCAP_SARS2) [[4]](#references).

\[1\]: https://www.uniprot.org/uniprotkb/P59594/entry
\[2\]: https://www.uniprot.org/uniprotkb/P59595/entry
\[3\]: https://www.uniprot.org/uniprotkb/P59637/entry
\[4\]: https://www.uniprot.org/uniprotkb/P0DTC9/entry