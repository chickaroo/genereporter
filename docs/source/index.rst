.. genereporter documentation master file, created by
   sphinx-quickstart on Tue Mar 12 18:56:32 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to :purple:`genereporter's` documentation!
========================================

:purple-bold:`genereporter` is a Python library Sami developed as part of her Bachelor's thesis. 

:purple-bold:`genereporter` can help researchers of any background quickly visualize and understand a gene of interest's expression levels in their specific dataset.

Using scRNA-sequencing data (for now), :purple-bold:`genereporter` generates automatic analysis reports and visualizations, granting wet-lab
researchers access to output usually generated manually by computational researchers. 

:purple-bold:`genereporter` has three different modes of analysis: 

* Cell type specific
* Gene regulatory network (GRN) specific
* Sample (patient) specific

:purple-bold:`genereporter` is available to implement as a no-code front-end interface for wet-lab researchers, or as a Python library for bioinformaticians to integrate into their pipelines.

Check out the :doc:`usage` section for further information.

.. note:: This is a work in progress. Feel free to open an issue on GitHub with any comments or questions. 

.. toctree::
   :maxdepth: 1
   
   usage
   cellpipeline
   grnpipeline
   samplepipeline
   Cell_Example
   GRN_Example
   Sample_Example
   scenicpipeline




