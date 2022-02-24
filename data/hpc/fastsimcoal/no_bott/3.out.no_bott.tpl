//Parameters for the coalescence simulation program : fastsimcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NI
NN
NS
//Samples sizes and samples age 
58
40
52
//Growth rates
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Current migration matrix 0
0 MIGN MIGS
MIGN 0 MIGO
MIGS MIGO 0
//No migration matrix 1
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
TDIV1 2 1 1 RESIZE1 0 1
TDIV2 1 0 1 RESIZE2 0 1
//Number of independent loci [chromosome] 
20 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
DNA 2000000 3.2e-8 1.2e-8 0.33
