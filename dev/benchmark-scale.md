# Scale benchmark (0.15.0 architecture)

Threads: 10. Synthetic Zipf corpora; two toy models per grid.
A dense J x W rare mask at the largest size would need 200 GB;
the sparse partition sizes below are the complete objects.

|     J|     W|      nnz| shards|ks  | partition_s| partition_mb| doc_index_s| optimal_topic_s| micro_dev| k_pick|
|-----:|-----:|--------:|------:|:---|-----------:|------------:|-----------:|---------------:|---------:|------:|
| 1e+05| 30000|  4928396|      1|5/8 |         0.2|         51.2|         0.1|           102.6|   -0.3612|      5|
| 5e+05| 50000| 24990158|      4|5/8 |         1.1|        226.4|         0.6|           809.4|   -0.3280|      5|
| 1e+06| 50000| 42378862|      4|5/8 |         1.7|        400.3|         1.0|          1616.8|   -0.2483|      5|
