grep -o "\[[^]]*\]" Data/arabidopsis-msh1-blast.xml | sort | uniq | sed 's/\[//g' | sed 's/\]//g' > tmp
grep -o "\[[^]]*\]" Data/dendronepthya-msh1-blast.xml | sort | uniq | sed 's/\[//g' | sed 's/\]//g' >> tmp
grep -o "\[[^]]*\]" Data/heliospora-msh1-blast.xml | sort | uniq | sed 's/\[//g' | sed 's/\]//g' >> tmp
grep -o "\[[^]]*\]" Data/saccharomyces-msh1-blast.xml | sort | uniq | sed 's/\[//g' | sed 's/\]//g' >> tmp
sort tmp | uniq > msh1-blast-species.txt
grep -o "\[[^]]*\]" Data/saccharomyces-mgm101-blast.xml | sort | uniq | sed 's/\[//g' | sed 's/\]//g' > mgm101-blast-species.txt
grep -o "\[[^]]*\]" Data/saccharomyces-mhr1-blast.xml | sort | uniq | sed 's/\[//g' | sed 's/\]//g' > mhr1-blast-species.txt

awk '{print $1;}' msh1-blast-species.txt | sort | uniq | grep -v "\[" > msh1-blast-genera.txt
awk '{print $1;}' mgm101-blast-species.txt | sort | uniq | grep -v "\[" > mgm101-blast-genera.txt
awk '{print $1;}' mhr1-blast-species.txt | sort | uniq | grep -v "\[" > mhr1-blast-genera.txt
awk '{print $1;}' reca-species-set.txt | sort | uniq | grep -v "\[" > reca-blast-genera.txt
cat msh1-blast-genera.txt mgm101-blast-genera.txt mhr1-blast-genera.txt reca-blast-genera.txt | sort | uniq > all-genera.txt
