# code to construct "barcodes" of gene presence/absence from the superset of all hits from our blast searches

all_names = read.table("all-species-set-clean.txt", sep="\t", stringsAsFactors = FALSE)
msh1_names = read.table("msh1-species-set.txt", sep="\t", stringsAsFactors = FALSE)
mgm101_names = read.table("mgm101-species-set.txt", sep="\t", stringsAsFactors = FALSE)
mhr1_names = read.table("mhr1-species-set.txt", sep="\t", stringsAsFactors = FALSE)
reca_names = read.table("reca-species-set.txt", sep="\t", stringsAsFactors = FALSE)
presence = rep(0, length(all_names))
all_set = cbind(cbind(cbind(cbind(all_names, presence), presence), presence), presence)
colnames(all_set) = c("Name", "msh1", "mgm101", "mhr1", "reca")
for(i in 1:length(all_set[,1]))
{
  if(length(grep(all_names[i,1], msh1_names) != 0)) { all_set$msh1[i] = 1; }
  if(length(grep(all_names[i,1], mgm101_names) != 0)) { all_set$mgm101[i] = 1; }
  if(length(grep(all_names[i,1], mhr1_names) != 0)) { all_set$mhr1[i] = 1; }
  if(length(grep(all_names[i,1], reca_names) != 0)) { all_set$reca[i] = 1; }
}
write.table(all_set, "barcodes.txt")