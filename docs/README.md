# Documentation for protpy

[![Documentation Status](https://readthedocs.org/projects/protpy/badge/?version=latest)](https://protpy.readthedocs.io/en/latest/?badge=latest)

The documentation for `protpy` is hosted on **ReadTheDocs** and is available [here](https://protpy.readthedocs.io/en/latest). It is built and published automatically from the repo.

Files
-----
* `index.rst` - landing page.
* `usage.rst` - installation and usage examples.
* `descriptors.rst` - every descriptor: formulas, parameters, output columns, deprecations (v1.6.0) and references; kept in sync with `DESCRIPTORS.md`.
* `api.rst` - API reference generated from the module docstrings.
* `contributing.rst` - contribution guide.
* `references.bib` - bibliography; every DOI in the descriptor registry must appear here (checked by `tests/test_helpers.py`).
* `conf.py` - Sphinx configuration (version 1.6.0).

Build the documentation locally
-------------------------------
Install the requirements from the project root:
```bash
pip install -e . -r docs/requirements.txt
```

Build the HTML in the /docs folder, treating warnings as errors:
```bash
cd docs && make html SPHINXOPTS="-W"
```

Open the generated page:
```bash
open _build/html/index.html
```

The full tutorial for using the ReadTheDocs platform is available [here](https://docs.readthedocs.io/en/stable/tutorial/index.html).
