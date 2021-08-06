Code for performing phylogeographic rarefaction/sensitivity analyses in Wilkinson et al. 'A year of genomic surveillance reveals how the SARS-CoV-2 pandemic unfolded in Africa'

rarefaction.py performs the rarefaction analysis subsampling genomes from either Africa (internal) or non-African countries (external)
NOTE: For this analysis we used the 2021-05-03_timetree.tre tree and metadata.tsv data

time_sensitivity.py performs the sensitivity analysis looking at how the the proportion of external introductions changes depending on whether sampling was temporally uniform or non-uniform through time
NOTE: For this analyses we used the newtimetree.nwk tree and metadata.tsv data

DateTimeUtils.py is a helper script that converts python DateTime objects to decimal years and back