# cell_cycle
A collection of code used for the undergraduate project on determining the role of the cell cycle stage in response to stress (Swain lab)

timeToLastNextPeak - determines the cell cycle stage at stress and sorts the cells as pre-START and post-START

measureMultipleLineage2 - calculates temporal fraction of the cell cycle between the fluorescent peak and the bud detection using one sliding window (only along the bud lineage source)

measureFreq - calculates temporal fraction of the cell cycle between the fluorescent peak and the bud detection using two sliding window (along the bud and fluorescence lineage sources)

plotMultipleLineage - plots information from multiple lineages and multiple cells on one graph 

nextavg - calculates average duration of 1st cell cycle during stress for all the cells with the same 'last' (same time from the stress to the previous reporter peak)
