library(XML)
library(ggplot2)

# set up labels for organisms, tissues, and genes
plants = c("arabidopsis", "barley", "maize", "medicago", "potato", "rice", "soybean")
sam_tissues = c("Shoot Apex", "Immature Inflorescence", "SAM", "Vegetative Bud", "S4 [Shoot Apex]", "SAM", "SAM")
control_tissues = c("Lea", "Lea", "Lea", "Lea", "Lea", "Lea", "Lea")
control_tissues = c("Roo", "Roo", "Roo", "Roo", "Roo", "Roo", "Roo")
genes = c("msh1", "reca2", "reca3", "why2", "odb1", "osb1", "osb4", "actin")

# number of genes we're testing
ncasegenes = 7

# control gene for general SAM elevation
controlgene = 8

df = data.frame()
dfout = data.frame()
counter = 0

# loop through genes we're testing
for(j in 1:ncasegenes)
{
  # loop through species
  for(i in 1:length(plants))
  {
    # try to read in the corresponding record
    l = -1
    fnamecase = paste("Data/", plants[i], "-", genes[j], ".html", sep="")
    fnamectl = paste("Data/", plants[i], "-", genes[controlgene], ".html", sep="")
    l = readHTMLTable(fnamecase)
    lctl = readHTMLTable(fnamectl)

    # if we read something
    if(length(l) != 0 && length(lctl) != 0)
    {
      # extract the appropriate parts of the data structure
      l = l[[1]]
      lctl = lctl[[1]]

      # first we look at the current test gene of interest
      # the inner grepl finds all l$Tissue labels matching the regex in our list of tissues (i.e. get SAM-related tissues)
      # we pull the corresponding mean expression values (col 3) and paste them together
      samexp = as.numeric(paste(l[grepl(sam_tissues[i],l$Tissue,fixed=TRUE),3]))

      # same for SD (col 4)
      samsds = as.numeric(paste(l[grepl(sam_tissues[i],l$Tissue,fixed=TRUE),4]))

      # same now for the opposite grep -- i.e. non-SAM-related tissues
      otherexp = as.numeric(paste(l[!grepl(sam_tissues[i],l$Tissue,fixed=TRUE),3]))
      othersds = as.numeric(paste(l[!grepl(sam_tissues[i],l$Tissue,fixed=TRUE),4]))

      # now the same for our control gene
      samexpctl = as.numeric(paste(lctl[grepl(sam_tissues[i],lctl$Tissue,fixed=TRUE),3]))
      samsdsctl = as.numeric(paste(lctl[grepl(sam_tissues[i],lctl$Tissue,fixed=TRUE),4]))
      otherexpctl = as.numeric(paste(lctl[!grepl(sam_tissues[i],lctl$Tissue,fixed=TRUE),3]))
      othersdsctl = as.numeric(paste(lctl[!grepl(sam_tissues[i],lctl$Tissue,fixed=TRUE),4]))

      # compute uncertainties
      samsd = sqrt(sum(samsds**2))/length(samexp)
      othersd = sqrt(sum(othersds**2))/length(otherexp)
      samsdctl = sqrt(sum(samsdsctl**2))/length(samexpctl)
      othersdctl = sqrt(sum(othersdsctl**2))/length(otherexpctl)
      # compute magnitude of overexpression in SAM relative to control gene
      ratio1 = mean(samexp)/mean(otherexp)
      ratio2 = mean(samexpctl)/mean(otherexpctl)
      ratio = (mean(samexp)/mean(otherexp))/(mean(samexpctl)/mean(otherexpctl))

      # compute uncertainties
      ratio1sd = sqrt(ratio1*((samsd/mean(samexp))**2 + (othersd/mean(otherexp))**2))
      ratio2sd = sqrt(ratio2*((samsdctl/mean(samexpctl))**2 + (othersdctl/mean(otherexpctl))**2)) 

      ratiosd = sqrt(ratio*((ratio1sd/ratio1)**2 + (ratio2sd/ratio2)**2))
     
      tmp = data.frame(fnamecase, ratio, ratiosd)
      df = rbind(df, tmp)
      dfout = rbind(dfout, data.frame(counter, plants[i], genes[j], ratio, ratiosd))
      counter = counter + 1
    }
  }
}

ggplot(df) +
   geom_bar( aes(x=fnamecase, y=ratio), stat="identity", fill="black", alpha=0.7) +
   geom_errorbar( aes(x=fnamecase, ymin=ratio-ratiosd, ymax=ratio+ratiosd), width=0.4, colour="black", alpha=0.9, size=1.3) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

write.table(dfout, "geneexp-summary.txt", sep = " ", col.names=FALSE, row.names=FALSE, quote=FALSE)