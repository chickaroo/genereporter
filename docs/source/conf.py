# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../../genereporter/'))

project = 'genereporter'
copyright = '2025, Samantha Bening'
author = 'Samantha Bening'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon',
              'sphinx.ext.autodoc', 
              'sphinx.ext.autosummary',
              "nbsphinx",
              "nbsphinx_link",]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_theme_options = {
    "source_repository": "https://github.com/chickaroo/genereporter",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#7C4DFF",
        "color-brand-content": "#7C4DFF",
        "color-problematic": "#7C4DFF",
        "color-api-background": "#efeff4",
        "color-api-background-hover": "#d4d4d9",
        "color-link": "#7C4DFF",
        "color-link--hover": "#7C4DFF",
        "color-link--visited": "#7C4DFF",
        # Add custom purple color for light mode
        "color-purple-custom": "#2F1D61", # darker purple for genereporter text

    },
    "dark_css_variables": {
        "color-brand-primary": "#7C4DFF",
        "color-brand-content": "#7C4DFF",
        "color-problematic": "#7C4DFF",
        "color-api-background": "#1e2124",
        "color-api-background-hover": "#45454a",
        "color-link": "#7C4DFF",
        "color-link--hover": "#7C4DFF",
        "color-link--visited": "#7C4DFF",
        # Add custom purple color for dark mode
        "color-purple-custom": "#8069BF",  # Lighter for dark mode

    },

}
html_static_path = ['_static']
html_css_files = ['custom.css']

rst_prolog = """
.. role:: purple
.. role:: purple-bold
"""
