-----------------------------------------------------------------------------------
  Supplementary Material to main paper: an example, source code, simulation files
-----------------------------------------------------------------------------------


1. `anExample.rmd`: An R notebook that contains two examples for our sequential sorting procedure. One with number of variables p=5, and another with p=10,000. 


2. `anExample.pdf`: The nice summary document for the R notebook.


3. `source/scorelingam/source.cpp`: The C++ source code for the sorting procedure, making use of linear algebra modules. This source is compiled when the .Rmd notebook cells are knit (compiled) into the PDF.
	- Formal R Package on Github, with Installation Instructions: https://gabriel-ruiz.github.io/scorelingam/.

4. `source/scorelingam/helperFunctions.r`: wrapper, data generating, and accuracy checking functions. This source is compiled by the .Rmd notebook.


5. `source/directlingam/wrapper.py`: wrapper for the `lingam` python module that implements the DirectLiNGAM procedure.


6. `simulations/section3.1`: Section 3.1 simulation files.


7. `simulations/section3.2`: Section 3.2 simulation files.


8. `simulations/section3.3`: Section 3.3 simulation files.

