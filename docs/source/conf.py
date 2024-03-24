# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../../genereporter/'))

project = 'genereporter'
copyright = '2024, Samantha Bening'
author = 'Samantha Bening'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 
              'sphinx.ext.autosummary',
              "nbsphinx",
              "nbsphinx_link",]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#7C4DFF",
        "color-brand-content": "#7C4DFF",
        "color-problematic": "#7C4DFF",
        "color-api-background": "#efeff4",
        "color-api-background-hover": "#d4d4d9",

    },
    "dark_css_variables": {
        "color-brand-primary": "#7C4DFF",
        "color-brand-content": "#7C4DFF",
        "color-problematic": "#7C4DFF",
        "color-api-background": "#1e2124",
        "color-api-background-hover": "#45454a",

    },

}
html_static_path = ['_static']
