# the data here includes hit lists from blastx searches and preconstructed text trees from NCBI Common Taxonomy Tool.
# see readme.txt for details on reconstructing these searches and trees

# construct "barcodes" of gene presence/absence in different organisms
R CMD BATCH construct-barcodes.R

# apply code to process individual trees to each tree we want to visualise
# "depths" -- to pull all leaves down to the same level -- are chosen manually
chmod +x phylo-individual.sh
./phylo-individual.sh commontree.txt 7
./phylo-individual.sh commontree-metazoa.txt 7
./phylo-individual.sh msh1-blastx-tree.txt 8
./phylo-individual.sh mgm101-blastx-tree.txt 9
./phylo-individual.sh mhr1-blastx-tree.txt 8
./phylo-individual.sh commontree-genera.txt 7
