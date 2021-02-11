# odna-segregation

Simulations and data analysis exploring oDNA segregation across eukaryotes.

Core technologies: bash, C, R, Gnuplot, Mathematica. The pipeline is rather strongly linked to bash -- implementation on other platforms is untested. Despite making up a large proportion of the code payload, Mathematica is only used for demonstrating symbolic algebra and visualisation; all analysis takes place using free platforms (and we are in the process of translating the existing Mathematica code).

All directories except `Analytics/` and `Development/` contain a `run-me.sh` file which runs the analysis. The root `./` also contains a master `run-me.sh` which runs each of these. The root also contains `plot-all.sh` which invokes individual plotting scripts in each directory to produce manuscript figures using Gnuplot (note that Figs 2 and S5, the phylogenies, are produced by Mathematica -- see `Phylogenetics/` below). Several directories contain `Data/` subdirectories containing data to be analysed.

`Analytics/`
------------
Mathematica notebooks demonstrating the algebra behind the stochastic theory.

`Development/`
--------------
Data and Gnuplot scripts for visualising observations of oDNA copy number during early mouse and plant development (Fig S4). Data from Refs. (1-9)

`Mouse/`
--------
Data, R, and Gnuplot scripts to fit and visualise models of mtDNA variance during mouse development and ageing (Figs 1C-E, S2). `model-fit-neutral-hb.R` performs model fitting in the zero-selection case (to the HB mouse model). `model-fit-selection-le.R` performs model fitting with nonzero selection (to the LE mouse model). `model-fit-elife.R` performs model fitting with some different zero selection mechanisms (to the NZB-BALB/C mouse model). Data from Refs. (1-3, 10)

`Partitioning/`
---------------
C code and Gnuplot scripts to simulate and visualise the genetic impact of physical partitioning of oDNA at cell divisions (Fig S3). `partitioning-simulation.c` is the simulation code. 

`Phylogenetics/`
----------------
Data, bash, R, and Mathematica scripts for visualising presence/absence of oDNA recombination genes across taxa (Figs 2, S5). Data is derived from previously-performed NCBI Gene searches https://www.ncbi.nlm.nih.gov/gene , blastx searches https://blast.ncbi.nlm.nih.gov/Blast.cgi?LINK_LOC=blasthome&PAGE_TYPE=BlastSearch&PROGRAM=blastx and taxonomic tree construction https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi , outlined in `readme.txt`. Various `.sh` bash scripts curate and clean the data; `construct-barcodes.R` assigns presence/absence barcodes to each species; `process-individual-tree.sh` invokes various small bits of code to (rather awkwardly) cast the constructed trees into a workable format. This pipeline is spaghetti code and will be refactored.

`SAM/`
------
Data, R code, and Gnuplot scripts for analysing and visualising gene expression in plant shoot apical meristems. Data is not raw but is from the University of Toronto's bio-analytic resource for plant biology (BAR) http://bar.utoronto.ca/ , compiling multiple sources especially Ref. 11 (BAR itself is described in Ref. 12). `parse-sam.R` performs the (simple) analysis.

`Stochastics/`
--------------
C code and Gnuplot scripts to simulate and visualise the output of stochastic simulations to test the theory (Figs 1B, S1). `stochastic-simulations.c` simulates the full system (with physical dynamics); `stochastic-simulations-simple.c` simulates a reduced subset of processes (turnover and partitioning).

References
----------
  1. Timothy Wai, Daniella Teoli, and Eric A Shoubridge. The mitochondrial DNA genetic bottleneck results from replication of a subpopulation of genomes. Nature genetics, 40(12):1484, 2008.
  2. Liqin Cao, Hiroshi Shitara, Takuro Horii, Yasumitsu Nagao, Hiroshi Imai, Kuniya Abe, Takahiko Hara, Jun-Ichi Hayashi, and Hiromichi Yonekawa. The mitochondrial bottleneck occurs without reduction of mtDNA content in female mouse germ cells. Nature genetics, 39(3):386, 2007.
  3. Lynsey M Cree, David C Samuels, Susana Chuva de Sousa Lopes, Harsha Karur Rajasimha, Passorn Wonnapinij, Jeffrey R Mann, Hans-Henrik M Dahl, and Patrick F Chinnery. A reduction of mitochondrial DNA molecules during embryogenesis explains the rapid segregation of genotypes. Nature genetics, 40(2):249, 957 2008.
  4. Tobias Preuten, Emilia Cincu, Jörg Fuchs, Reimo Zoschke, Karsten Liere, and Thomas Börner. Fewer genes than organelles: extremely low and variable gene copy numbers in mitochondria of somatic plant cells. The Plant Journal, 64(6):948–959, 2010.
  5. Hideki Takanashi, Takayuki Ohnishi, Mirai Mogi, Takashi Okamoto, Shin-ichi Arimura, and Nobuhiro Tsutsumi. Studies of mitochondrial morphology and DNA amount in the rice egg cell. Current genetics, 56(1):33–41, 980 2010.
  6. H Kuroiwa, T Ohta, and T Kuroiwa. Studies on the development and three-dimensional reconstruction of giant mitochondria and their nuclei in egg cells of Pelargonium zonale ait. Protoplasma, 192(3-4):235–244, 1996.
  7. Gayle K Lamppa and Arnold J Bendich. Changes in mitochondrial DNA levels during development of pea (Pisum sativum l.). Planta, 162(5):463–468, 1984.
  8. Dan-Yang Wang, Quan Zhang, Yang Liu, Zhi-Fu Lin, Shao-Xiang Zhang, Meng-Xiang Sun, et al. The levels of male gametic mitochondrial dna are highly regulated in angiosperms with regard to mitochondrial inheritance. The Plant Cell, 22(7):2402–2416, 2010.
  9. Long Gao, Xue Guo, Xue-Qiong Liu, Li Zhang, Jilei Huang, Li Tan, Zhen Lin, Shingo Nagawa, and Dan-Yang Wang. Changes in mitochondrial DNA levels during early embryogenesis in Torenia fournieri and Arabidopsis thaliana. The Plant Journal, 95(5):785–795, 2018.
  10. Joerg P Burgstaller, Thomas Kolbe, Vitezslav Havlicek, Stephanie Hembach, Joanna Poulton, Jaroslav Piálek, Ralf Steinborn, Thomas Rülicke, Gottfried Brem, Nick S Jones, et al. Large-scale genetic analysis reveals mammalian mtDNA heteroplasmy dynamics and variance increase through lifetimes and generations. Nature communications, 9(1):1–12, 2018.
  11. Markus Schmid, Timothy S Davison, Stefan R Henz, Utz J Pape, Monika Demar, Martin Vingron, Bernhard Schölkopf, Detlef Weigel, and Jan U Lohmann. A gene expression map of Arabidopsis thaliana development. Nature genetics, 37(5):501–506, 2005.
  12. Debbie Winter, Ben Vinegar, Hardeep Nahal, Ron Ammar, Greg V Wilson, and Nicholas J Provart. An electronic fluorescent pictograph browser for exploring and analyzing large-scale biological data sets. PloS ONE, 2(8):e718, 2007.
  
