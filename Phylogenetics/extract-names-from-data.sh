### having downloaded summaries of gene searches in eukaryotes

grep -v "Drosophila\|Rattus\|Mus\|Homo\|[Vv]irus\|[Pp]hage" Data/msh1-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > msh1-species-set.txt
grep -v "[Vv]irus" Data/mgm101-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > mgm101-species-set.txt
grep -v "[Vv]irus" Data/reca-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > reca-species-set.txt
grep -v "Haloprofundus\|Homo\|Danio" Data/mhr1-results.txt | awk 'BEGIN{FS="\t";n=0;}{if(n++ > 0){split($2, arr, " "); print arr[1], arr[2];}}' | sort | uniq > mhr1-species-set.txt

# issues
#Organism name 'Exserohilum turcica' not found
#Organism name 'Keratoisidinae sp.' not found
#Organism name 'Pyrus x' not found
# ... so remove these awkward records
cat msh1-species-set.txt mgm101-species-set.txt reca-species-set.txt mhr1-species-set.txt | sort | uniq | grep -v "Exserohilum" | grep -v "Keratoisidinae" | grep -v "Pyrus x" > all-species-set.txt

# these species have synonyms
grep -v "costaricaensis\|nomius\|lacticoffeatus\|fumosorosea\|haematococca\|Candida\|gattii\|japonicum\|fasciculatum\|pallidum\|lamblia\|lucimarinus\|Escherichia" all-species-set.txt > all-species-set-clean.txt


# synonyms from taxonomy tool
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

