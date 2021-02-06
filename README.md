# odna-segregation

Simulations and data analysis exploring oDNA segregation across eukaryotes.

Core technologies: bash, C, R, Gnuplot, Mathematica. The pipeline is rather strongly linked to bash -- implementation on other platforms is untested. Despite making up a large proportion of the code payload, Mathematica is only used for demonstrating symbolic algebra and visualisation; all analysis takes place using free platforms (and we are in the process of translating the existing Mathematica code).

All directories except `Analytics/` and `Development/` contain a `run-me.sh` file which runs the analysis. The root `./` also contains a master `run-me.sh` which runs each of these. The root also contains `plot-all.sh` which invokes individual plotting scripts in each directory to produce manuscript figures using Gnuplot (note that Figs 2 and S5, the phylogenies, are produced by Mathematica -- see `Phylogenetics/` below). Several directories contain `Data/` subdirectories containing data to be analysed.

`Analytics/`
------------
Mathematica notebooks demonstrating the algebra behind the stochastic theory.

`Development/`
--------------
Data and Gnuplot scripts for visualising observations of oDNA copy number during early plant development (Fig S4)

`Mouse/`
--------
Data, R, and Gnuplot scripts to fit and visualise models of mtDNA variance during mouse development and ageing (Figs 1C-E, S2). `model-fit-neutral.R` performs model fitting in the zero-selection case (to the HB mouse model). `model-fit-selection.R` performs model fitting with nonzero selection (to the LE mouse model).

`Partitioning/`
---------------
C code and Gnuplot scripts to simulate and visualise the genetic impact of physical partitioning of oDNA at cell divisions (Fig S3). `partitioning-simulation.c` is the simulation code. 

`Phylogenetics/`
----------------
Data, bash, R, and Mathematica scripts for visualising presence/absence of oDNA recombination genes across taxa (Figs 2, S5). Data is not raw data but is derived from previously-performed NCBI Gene searches https://www.ncbi.nlm.nih.gov/gene , blastx searches https://blast.ncbi.nlm.nih.gov/Blast.cgi?LINK_LOC=blasthome&PAGE_TYPE=BlastSearch&PROGRAM=blastx and taxonomic tree construction https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi , outlined in `readme.txt`. Various `.sh` bash scripts curate and clean the raw data; `construct-barcodes.R` assigns presence/absence barcodes to each species; `process-individual-tree.sh` invokes various small bits of code to (rather awkwardly) cast the constructed trees into a workable format. This pipeline is spaghetti code and will be refactored.

`SAM/`
------
Data, R code, and Gnuplot scripts for analysing and visualising gene expression in plant shoot apical meristems. Data is not raw but is from the University of Toronto's bio-analytic resource for plant biology (BAR) http://bar.utoronto.ca/ . `parse-sam.R` performs the (simple) analysis.

`Stochastics/`
--------------
C code and Gnuplot scripts to simulate and visualise the output of stochastic simulations to test the theory (Figs 1B, S1). `stochastic-simulations.c` simulates the full system (with physical dynamics); `stochastic-simulations-simple.c` simulates a reduced subset of processes (turnover and partitioning).
