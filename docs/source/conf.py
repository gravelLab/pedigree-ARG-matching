# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../../alignment'))
sys.path.insert(0, os.path.abspath('../../graph'))
sys.path.insert(0, os.path.abspath('../../scripts'))
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pedigree-ARG-matching'
copyright = '2024, Andrii Serdiuk and Simon Gravel'
author = 'Andrii Serdiuk, Simon Gravel'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.autosectionlabel',
]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,  # Excludes private members by default
    'special-members': '__init__',
}

napoleon_google_docstring = True
napoleon_numpy_docstring = False

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
