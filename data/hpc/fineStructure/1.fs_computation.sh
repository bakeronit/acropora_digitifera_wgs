#!/usr/bin/bash

#use hpc mode
# there several stages
# stage 1, estimate cp parameters
# we used 20% samples (-s1indfrac 0.2) to perform chromopainter parameter estimate
# post s1: estimating the parameters mu and Ne from output of s1
# s2 chromopainter painting
# post s2: chromcombine
# s3: fs msmc inference  -s3iters will assign half to burnin and half to sampling.
# s4: fs tree inference

#./fs_4.1.1/fs adigitifera.cp -n -phasefiles phase_files/* -recombfiles recom/* -idfile samples.id \
#-hpc 1 -s1indfrac 0.3 -s3iters 2000000 -go

# commands for stage1 were prepared
#cat adigitifera/commandfiles/commandfile1.txt | parallel -j 20


#~/soft/fs_4.1.0/fs adigitifera.cp -go
#cat adigitifera/commandfiles/commandfile2.txt | parallel -j 20

#~/soft/fs_4.1.0/fs adigitifera.cp -go
#cat adigitifera/commandfiles/commandfile3.txt | parallel -j 2

#~/soft/fs_4.1.0/fs adigitifera.cp -go
#cat adigitifera/commandfiles/commandfile4.txt | parallel -j 2

#if you can install gui successfully, but I couldn't
#finegui -c adigitifera_linked.chunkcounts.out -m adigitifera_linked_mcmc.xml -t adigitifera_linked_tree.xml -m2 adigitifera/stage3/adigitifera_linked_mcmc_run1.xml -t2 adigitifera/stage4/adigitifera_linked_tree_run1.xml

#or use Rcode
#extract more info:
~/soft/fs_4.1.0/fs fs -X -Y -e X2 adigitifera_linked.chunkcounts.out adigitifera_linked_tree.xml adigitifera_linked.mapstate.csv
~/soft/fs_4.1.0/fs fs -X -Y -e X2 adigitifera_linked.chunkcounts.out adigitifera_linked_mcmc.xml adigitifera_linked.meancoincidence.csv
