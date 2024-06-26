<h1 align="center">genereporter</h1>
<p align="center">
  VEO-IBD Gene-Level Report Generator
  <br>
  <br>
  <img src="data2/geneformer_logo.jpeg" alt="Demo image" width="300" height="300">
</p>

<!-- start elevator-pitch -->

genereporter automatically generates an analysis and visualization of the gene of interest. Across three different areas&mdash;cell type specific, gene set enrichment analysis, and sample (patient) specific&mdash;genereporter can help researchers quickly access and understand scRNA-seq data easily, which usually is a heavily computational task. However, genereporter is also built to be customizable for users that wish to dive into the code. 

The <a href="https://genereporter.readthedocs.io/en/latest/index.html">documentation</a> includes sample automatic pipelines as well as the API. 


<!-- end elevator-pitch -->

## Quickstart

<!-- start quickstart -->

Installing genereporter is straight forward&mdash;just clone the repository and install the package using pip. 

First, clone the github repository.
<pre>
  <code> $ git clone https://github.com/chickaroo/genereporter </code>
</pre>

Next, you might want to make a conda environment to install the package in.

<pre>
  <code> 
    $ conda create -n genereporter python=3.11
    $ conda activate genereporter
  </code>
</pre>

Then, navigate to the directory and install the package using pip.

<pre>
  <code>
    $ cd genereporter
    $ pip install .
  </code>
</pre>

This will install the package and all of its dependencies using pip.

<!-- end quickstart -->
