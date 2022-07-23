# Find the best partition scheme
iqtree  -p pomo_partition.nex -m MF+MERGE -redo -nt AUTO

# Run phylogeny with ultrafast bootstraps for best scheme
iqtree  -p pomo_partition.nex.best_scheme.nex -nt AUTO -B 1000 -redo

