# Index construction with chromosome size balancing

`divide_and_index.sh` is the wrapper script for balanced index construction. 
It takes the reference genome and outputs a partitioned index optimised for reduced peak memory.
  
Its usage is as follows.

```
usage : ./divide_and_index.sh <reference.fa> <num_parts> <out.idx> <minimap2_exe> <minimap2_profile>

reference.fa - path to the fasta file containing the reference genome
num_parts - number of partitions in the index
out.idx - path to the file to which the index should be dumped
minimap2_exe - path to the minimap2 executable
minimap2_profile - minimap2 pre-set for indexing (map-pb or map-ont)

Example : ./divide_and_index.sh hg19.fa 4 hg19.idx minimap2 map-ont

```

Functionality of `divide_and_index.sh` is as follows.
1. Compiling `divide.c` using `gcc` to produce `divide`.
2. Call the compiled binary `divide` to split the reference genome into partitions such that the total length of chromosomes in each partition are approximately equal.
3. Calling the minimap2 binary separately on each reference partition to produce a separate index file for each partition.
4. Combining all the index files to produce a single partitioned index file.


