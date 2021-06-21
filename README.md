# quokka
QUality Of Known Klusters Assessed


## To do

### Lineage QC

- check SNPs for each lineage, are there some sequences that have missing sites at what otherwise be defining SNPs
- are there any obvious sequences that look like they belong in different lineages
- check ambiguity content
- check monophyly within a tree (figtree2, usher and nextclade), or within updown-topranking gofasta
- check decision tree for each lineage, see if any conflicting decisions or any that dont seem to fit
- get 95% consensus of a lineage, check if outliers and remove
- Do all records on pango-designation match the gisaid metadata file? If not which ones and why (e.g. misspellings, or wrong year etc).
- Do all lineages have records?
- Are there sequences designated a lineage that scorpio disagrees with?
- Are there sequences that are in the metadata that scorpio calls that we should
- Do all sequences in a certain lineage have the set of mutations that they should have? Are there ones that really stick out?
- Primer sites/ problematic sites check
- Input fasta file, output ambiguity content

### Lineage discovery

- diversity scanning
- polecat -> fast spreading, locations
- geographic clustering
- How to flag new lineages?
