The [gene]-results.txt files were generated from NCBI's Gene database. The gene was searched for online and the results "sent to" a tabular text file. An example search is https://www.ncbi.nlm.nih.gov/gene/?term=msh1 .

Following hit list generation, the data are cleaned using extract-names-from-data.sh (this removes some known false positives) and NCBI's Common Taxonomy Tree tool is used to construct text tree taxonomies for the remaining species.

The blastx hit lists in Data/ were generated with the following searches. All are blastx of the given sequence against the nr database, against eukaryotes unless otherwise stated. Searches were performed online and the results "Downloaded" to XML format.

Arabidopsis Msh1 (NM_113339.4), limit E < 1e-50 [to avoid Msh2 etc]: 668 hits, viridiplantae and protists; doesn't hit corals

Dendronephthya putteri Msh1 (NC_036022.1 (6348..9287)) vs metazoa, limit E < 1e-50: 2593 hits, soft corals, anthozoans, sea pens, blue corals

Heliopora coerulea msh1 (NC_020375.1 (2498..5473)) limit E < 1e-20 threshold: 2551 hits; blue corals etc

Saccharomyces cerevisiae S288C msh1 (NC_001140.6 (349574..352453)), limit E < 1e-50: 1927 hits, fungi
  
Saccharomyces cerevisiae S288C mgm101 (NC_001142.9 (700882..701691)), limit E < 1: 1278 hits, lots of fungi, quercus, naegleria, tieghemostelium (slime mold), orbicella (coral), actinia (sea anemone), stony corals, placozoans, anemones, etc

Saccharomyces cerevisiae S288C mhr1 NC_001136.10 (1055212..1055892), limit E < 1: 146 hits, fungi

These results are processed with process-blast-results.sh .
