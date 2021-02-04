### having downloaded summaries of gene searches in eukaryotes

grep -v "Drosophila\|Rattus\|Mus\|Homo\|[Vv]irus\|[Pp]hage" Data/msh1-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > msh1-species-set.txt
grep -v "[Vv]irus" Data/mgm101-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > mgm101-species-set.txt
grep -v "[Vv]irus" Data/reca-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > reca-species-set.txt
grep -v "Haloprofundus\|Homo\|Danio" Data/mhr1-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > mhr1-species-set.txt

# remove some awkward records
cat msh1-species-set.txt mgm101-species-set.txt reca-species-set.txt mhr1-species-set.txt | sort | uniq | grep -v "Exserohilum" | grep -v "Keratoisidinae" | grep -v "Pyrus x" > all-species-set.txt

# these species have synonyms
#grep -v "costaricaensis\|nomius\|lacticoffeatus\|fumosorosea\|haematococca\|Candida\|gattii\|japonicum\|fasciculatum\|pallidum\|lamblia\|lucimarinus\|Escherichia" all-species-set.txt > all-species-set-clean.txt

echo "Please construct taxonomic tree using all-species-set.txt then run XXX"
# tr -d "\n" < phyliptree.phy | sed 's/:4//g' | sed 's/ /_/g' | sed "s/'//g" > strip-tree.phy

# blast queries
blastx -db nr -query Data/arabidopsis-msh1.fasta -out msh1-blast-1.txt -outfmt "6 ssciname qseqid sseqid evalue" -evalue 1e-50 -max_target_seqs 1000 -remote
blastx -db nr -query Data/coral-msh1.fasta -out msh1-blast-2.txt -outfmt 6 -evalue 1e-50 -max_target_seqs 5000 -remote
blastx -db nr -query Data/yeast-mgm101.fasta -out mgm101-blast.txt -outfmt 6 -evalue 1 -max_target_seqs 5000 -remote

# issues
#Organism name 'Exserohilum turcica' not found
#Organism name 'Keratoisidinae sp.' not found
#Organism name 'Pyrus x' not found


#Aspergillus costaricensis (Aspergillus costaricaensis)
#Aspergillus nomiae (Aspergillus nomius)
#Aspergillus niger (Aspergillus lacticoffeatus)
#Cordyceps fumosorosea (Isaria fumosorosea)
#[Nectria] haematococca (Nectria haematococca)
#[Candida]
#Cryptococcus gattii VGI (Cryptococcus gattii)
#Corallium japonicum (Paracorallium japonicum)
#Cavenderia fasciculata (Dictyostelium fasciculatum)
#Heterostelium pallidum (Polysphondylium pallidum)
#Giardia intestinalis (Giardia lamblia)
#Ostreococcus sp. 'lucimarinus' (Ostreococcus lucimarinus)
#Escherichia coli

