# Scale benchmark (sparse-partition architecture)

Threads: 10. Synthetic Zipf corpora; two toy models per grid.
A dense J x W rare mask at the largest size would need 200 GB;
the sparse partition sizes below are the complete objects.

|     J|     W|      nnz| shards|ks  | partition_s| partition_mb| doc_index_s| optimal_topic_s| micro_dev| k_pick|
|-----:|-----:|--------:|------:|:---|-----------:|------------:|-----------:|---------------:|---------:|------:|
| 1e+05| 30000|  4928396|      1|5/8 |         0.1|         51.2|         0.1|            40.0|   -0.3612|      5|
| 5e+05| 50000| 24990158|      4|5/8 |         0.6|        226.4|         0.4|           357.2|   -0.3280|      5|
| 1e+06| 50000| 42378862|      4|5/8 |         0.9|        400.3|         0.8|           739.4|   -0.2483|      5|

## Standing notes

- 0.17.0 kernel changes: the envelope initial boundary is a fixed 1024
  and the collapsed-tail term uses exact complements (the former
  full-vocabulary totals were identically 1), removing the dominant
  non-BLAS cost of optimal_topic at large W. Against the 0.15.0 numbers
  (102.6 / 809.4 / 1616.8 s) the optimal_topic rows above are 2.2x to
  2.6x faster; micro_dev and k_pick are unchanged to all printed digits.
  The remaining cost is the per-model BLAS product over the vocabulary,
  inherent to the statistic.
- Doc-kernel gemv batching: evaluated and rejected. The merge-join doc
  kernel reads contiguous phi columns with a gathered theta row; a per-
  document gemv would first gather phi(:, NR_j) into a dense buffer,
  copying the same memory the dot products read, and the kernel runs at
  ~1 s per million documents, three orders below the pipeline bottleneck
  (optimal_topic). No batching is warranted.
