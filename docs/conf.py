# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'protpy'
copyright = '2026, AJ McKenna'
author = 'AJ McKenna'
release = '1.5.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

# autodoc settings: show both class and __init__ docstrings, preserve param order
autodoc_member_order = 'bysource'
autodoc_default_options = {'members': True}
napoleon_numpy_docstring = True
napoleon_google_docstring = False

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


#protpy docstrings write parameters as ":name: type" under a "Parameters" heading; napoleon
#expects numpy's "name : type", so rewrite those lines before napoleon (priority 500) runs
import re as _re
_FIELD = _re.compile(r"^(\s*):(\w+):\s*(.*)$")

def _house_style_to_numpy(app, what, name, obj, options, lines):
    section = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped in ("Parameters", "Returns", "Raises") and i + 1 < len(lines) and set(lines[i + 1].strip()) <= {"=", "-"}:
            section = stripped
        match = _FIELD.match(line)
        if section in ("Parameters", "Returns") and match:
            lines[i] = f"{match.group(1)}{match.group(2)} : {match.group(3)}"

def setup(app):
    app.connect("autodoc-process-docstring", _house_style_to_numpy, priority=400)
