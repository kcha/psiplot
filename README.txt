Using PSI Plotter script
========================

Dependencies
------------

* preprocess_sample_colors.R - a helper script for re-organizing the PSI input
data

Input File Requirements
------------------------

Two input files are required:

1) The PSI input data - this contains the unmodified PSI data from Manuel's pipeline.
2) Tissue grouping files (e.g. Tissues.Hsa.txt or Tissues.Mmu.txt)


Execution
---------

From command-line:

~~~~
> Rscript PSI_Plotter.R --args PSI_Input.tab Tissues.Mmu.txt
~~~~

PSI_Plotter.R will output a pdf file named after the input file (e.g.
PSI_Input.PSI_plots.pdf). The pdf will be saved in the same directory as the
input file. 
