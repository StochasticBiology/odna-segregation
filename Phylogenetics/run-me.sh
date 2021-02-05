# the data here includes search results, hit lists from blastx searches and preconstructed text trees from NCBI Common Taxonomy Tool.
# see readme.txt for details on reconstructing these searches and trees

# get names of hits from Gene search
chmod +x extract-names-from-data.sh
./extract-names-from-data.sh

# process blastx results
chmod +x process-blast-results.sh
./process-blast-results.sh

# to reconstruct the pipeline from scratch, use the outputs of those scripts to construct Common Taxonomy Trees as required below. these are included, pre-made, in this repo for convenience.
cp Data/*tree*txt .

# construct "barcodes" of gene presence/absence in different organisms
R CMD BATCH construct-barcodes.R

# apply code to process individual trees to each tree we want to visualise
# "depths" -- to pull all leaves down to the same level -- are chosen manually
chmod +x process-individual-tree.sh
chmod +x tree-to-old-format.sh
chmod +x pad-tree.sh
./process-individual-tree.sh commontree.txt 7
#./process-individual-tree.sh commontree-metazoa.txt 7
./process-individual-tree.sh msh1-blastx-tree.txt 8
./process-individual-tree.sh mgm101-blastx-tree.txt 9
./process-individual-tree.sh mhr1-blastx-tree.txt 8
#./process-individual-tree.sh commontree-genera.txt 7

# after completion, use Mathematica notebook for visualisation (to be refactored)
