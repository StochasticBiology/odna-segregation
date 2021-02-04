# apply code to process individual trees to each tree we want to visualise
# "depths" -- to pull all leaves down to the same level -- are chosen manually

R CMD BATCH construct-barcodes.R
chmod +x phylo-individual.sh
./phylo-individual.sh commontree.txt 7
./phylo-individual.sh commontree-metazoa.txt 7
./phylo-individual.sh msh1-blastx-tree.txt 8
./phylo-individual.sh mgm101-blastx-tree.txt 9
./phylo-individual.sh mhr1-blastx-tree.txt 8
./phylo-individual.sh commontree-genera.txt 7
