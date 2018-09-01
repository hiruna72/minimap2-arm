#Evaluation scripts for partitioned indexes

##`synth_eval.sh`

The bash shell script for evaluating the mapping accuracy from the synthetic reads.
Performs the mapping using minimap2 and evaluation using Paftools for 
single-idx, part-idx-merged and part-idx-no-merge.
See the comments in  the script for more information.

###`synth_plot.m`
The Matlab script that generate the plots the fraction of mapped reads against the error rate of mapped reads. 
Should be run after `synth_eval.sh`. Tested on Matlab 2017. 

##`chimeric_eval.sh`
The shell script for mapping the chromothriptic read with different options.
Performs the mapping using minimap2 for single-idx, part-idx-merged and part-idx-no-merge.
See the comments in  the script for more information.

###`chimeric_plot.m`
The Matlab script that generates the plots to visualise the mappings of a chromothriptic read against the reference.
Should be run after `chimeric_eval.sh`. Tested on Matlab 2017. 

##`real_get_stat.sh`
The shell script for mapping a set of real reads and deriving mapping statistics.
Performs the mapping using minimap2 for single-idx, part-idx-merged and part-idx-no-merge.
See the comments in  the script for more information.

##`mem_usage_no_chr_balance.sh`
The bash script that evaluates the peak memory usage for an index without chromosome balancing.
Uses the inbuilt indexing option in minimap2 

##`Mem_usage_with_chr_balance.sh`
The bash script that evaluates the peak memory usage for the index with chromosome balancing.
Uses `divide_and_index.sh` to construct the index.

