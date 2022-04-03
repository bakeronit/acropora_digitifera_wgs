# See https://github.com/ComparativeGenomicsToolkit/cactus

singularity run -B $(pwd):/data cactus_v2.0.5.sif halSnps /data/corals.hal adig amil,aten --tsv /data/snps.tsv


# Our ancestral allele calculations assume that all ancestors align at any position. 
# This allows us to check this assumption 
#
singularity run -B $(pwd):/data cactus_v2.0.5.sif halAlignmentDepth /data/corals.hal adig --targetGenomes amil,aten --outWiggle corals.wig
