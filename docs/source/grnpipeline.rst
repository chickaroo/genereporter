GRN Pipeline API
--------------------------------

GRNPipeline is a pipeline for gene regulatory network analysis and visualization, specifically for the gene of interest. 
It is designed to be used with the output of the SCENIC package, specifically with the adjacencies and regulons files. More gene sets are also implemented, specifically from the Reactome database.
This pipeline assumes the SCENIC analysis has already been run, and the output files are provided as input to the pipeline.
For more details on running the SCENIC pipelines, see :doc:`scenicpipeline`.
See the :doc:`GRN_Example` for a detailed example of how to use it.

.. autoclass:: grn_pipeline.GRNPipeline
    :members:
    :special-members: __init__
    :no-undoc-members:
    :no-inherited-members:
