.. genereporter documentation master file, created by
   sphinx-quickstart on Tue Mar 12 18:56:32 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to genereporter's documentation!
========================================

**genereporter** is a Python library Sami developed as part (or most) or her Bachelor's thesis. 
genereporter focuses on generating an in-depth, gene-level view of a specific gene of interest. 
Using mostly scRNA-seq data (for now), genereporter provides analysis methods that can facilitate wet-lab
researchers accessability of typically computational tasks. genereporter has three different methods of analysis: 
cell type specific, gene set enrichment analysis, and sample (patient) specific. Through both the python library and the
automized Jupyter Notebook pipelines, genereporter is built to be user-friendly and easy to use for both 
bioinformaticians and wet-lab researchers. 
Check out the :doc:`usage` section for further information.

.. note:: This is a work in progress. The documentation and the library may still contain bugs and errors. 

.. toctree::
   :maxdepth: 1
   
   usage
   cellpipeline
   grnpipeline
   samplepipeline
   Cell_Example
   GRN_Example
   Sample_Example




