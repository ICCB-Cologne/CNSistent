import os
import sys
sys.path.insert(0, os.path.abspath('../'))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    "nbsphinx",
]

html_theme = 'sphinx_rtd_theme'