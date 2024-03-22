<h1 align="center">genereporter</h1>
<p align="center">
  VEO-IBD Gene-Level Report Generator
  <br>
  <img src="data/geneformer_logo.jpeg" alt="Demo image" width="300" height="300">
</p>
<!-- <a href="https://pradyunsg.me/furo/">
  <img align="center" src="https://github.com/pradyunsg/furo/raw/main/docs/_static/demo.png" alt="Demo image">
</a> -->

## Elevator pitch

<!-- start elevator-pitch -->

- **Automatic Gene Report** --- genereporter automatic generates an analysis and visualization of the gene of interest. Across three different areas--cell type specific, gene set enrichment analysis, and sample (patient) specific--genereporter can help researchers quickly access and understand the scRNA-seq data and information easily, which usually is a heavily computational task. However, genereporter is built to be customizable and adaptable to users that wish to dive into the code. 
 
- **Sami Worked Very Hard On This!!** 


<!-- end elevator-pitch -->

## Quickstart

<!-- start quickstart -->

Installing genereporter is straight forward--just clone the repository and install the package using pip. 

First, clone the github repository.

<code> $ git clone https://github.com/chickaroo/genereporter </code>

Next, you might want to make a conda environment to install the package in.

<code>$ conda create -n genereporter python=3.11 <br> $ conda activate genereporter</code>

Then, navigate to the directory and install the package using pip.

<code>$ cd genereporter <br> $ pip install . </code>

This will install the package and all of its dependencies using pip. The standard dependencies are
only geared towards the core functionality of the package. If you want to generate any of the GRN data
using <a href="https://scenic.aertslab.org/">SCENIC</a> for example, you will need to install additional dependencies according to the pySCENIC <a href="https://pyscenic.readthedocs.io/en/latest/installation.html">documentation</a>.

You can now test if the package is installed correctly by running the following command:

<code>$ python <br> &gt;&gt;&gt;import genereporter </code>

If you don't get any errors, the package is installed correctly.

<!-- end quickstart -->



