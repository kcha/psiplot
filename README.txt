Using PSI Plotter (aka "Nuno plots") script
===========================================

Dependencies
------------

* preprocess_sample_colors.R - a helper script for re-organizing the PSI input
data

Input File Requirements
------------------------

Two input files are required:

1) The PSI input data - this contains the unmodified PSI data from Manuel's pipeline.
2) Tissue grouping files (e.g. Tissues.Hsa.txt or Tissues.Mmu.txt)
    - A copy of both human and mouse files can be found in the test_data
      sub-directory.
    - More details about the formatting of this file is described in the
      document preprocess_sample_colors.pdf.

      Briefly, this file is made up of 4 tab-deliminted columns:
      Order    SampleName    GroupName    RColorCode

      * Order takes an integer and defines ordering of the sample on the
      PSI Plot from right (order=1) to left. 
      * SampleName is the name of the sample. This MUST match the column name in
      the PSI input data.
      * GroupName is used to group similar samples together, such as neural
      tissues.
      * RColorCode is the index of the vector from colors().


Execution
---------

From command-line:

~~~~
> Rscript PSI_Plotter.R --args \
    test_data/INCLUSION_LEVELS-ALL3m-Mmu89-SELECTED.test.tab \
    test_data/Tissues.Mmu.txt
~~~~

PSI_Plotter.R will output a pdf file named after the input file and
will be saved in the same directory as the input file (e.g.
test_data/INCLUSION_LEVELS-ALL3m-Mmu89-SELECTED.test.PSI_plots.pdf)
