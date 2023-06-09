# cal-sse-Jaccard
It is used for calculating the Sum of Seqared Error (SSE) between the estimated Jaccard and the true Jaccard.

The `true Jaccard` are computed by the RabbitKSSD by shuffle file with a dimention reduction level of 0 (e.g. L0K8.shuf), which means all the $k$-mers are saved in the sketch, and the Jaccard index computed between these unsampled sketches are served as `true Jaccard`.

### input
1. --standard-dist: the distance file generated from the RabbitKSSD with dimention reduction level of 0, for each line, it contains five columns, splited by `\t`.
```
query_name, ref_name, common|size0|size1, Jaccard, mash_distance
```
2. --rabbit-dist: the distance file generated by RabbitKSSD with the sampled sketch.
for each line, it contains five columns, splited by `\t`.
```
query_name, ref_name, common|size0|size1, Jaccard, mash_distance
```
3. --bindash-dist: the distance file generated by the bindash.
Each line has file tab-separated fields:
```
query_name, target_genome, mutation_distance, p_value, Jaccard
```
4. --dashing2-dist: the distance file generated by the dashing2.
It contains a distance matrix and the format is as follows:
```
first line: application descirption
second line: application options
third line: all genome file list. the reference files and the query files.
next lines: containing N+1 fields, where 1 means reference genome name, N is the number of the query genomes.
So the next lines consist of <ref_name, N double Jaccards between this reference and all N querys>
```

### output
1. --output: output comparison file between estimated Jaccards and true Jaccard.  
**For the spcific usage, there are only one reference genome (the number of query genomes are not 1), the output files only contain the query_name but no ref_name.**  
Each lines contains these fields:
```
query_name, true_Jaccard, estimated_Jaccard, squared_error
```
2. for muitiple input distance files, there are multiple output results.  
For example, if there are `--rabbit-dist` and `--bindashing-dist` for the input, there will be `output_rabbit.res` and `output_bindash.res` as well.

3. The Sum of squared error between the true Jaccard and estimated Jaccard will be output to the screen.



