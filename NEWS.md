# OpTop 0.6.0

### Minor Changes

* `agg_topic_stability()` can now compute smoothed tests and plot the results
accordingly.

* Support for final convertion to a `tibble` table spreaded out to all
functions,

* All eligible functions get better plots.


---

# OpTop 0.6.0

### Major Changes

* `optimal_topic()` gains the parameter `q` which allows to select the quantile
of the cumulative probability of word weights to consider as relevant.

* `optimal_topic()` now finds the optimal number of topics either by significance
levels or by forcing the algorithm to reach the global minimum. This is 
controlled by the new parameter `alpha`. 

### Minor Changes

* `optimal_topic()` drops both `threshold` and `q_type`.

* In `optimal_topic()`, `convert` now supports `tibble` structure.

---

# OpTop 0.5.0

### Major Changes

* Function `agg_topic_stability()` has been widely improved.

* All functions which return a test now gain the new argument `do_plot`. This
plot the test statistic as a function of the number of topics.

* The argument `test` has been removed from `topic_stability()` which now 
returns only the aggregate statistic. 

* The argument `compute_res` has been finally removed from `topic_stability()`.

* `topic_stability()` now returns either a data.frame or a data.table with 
the LDA specifiction associated to each statistic (i.e. column `topic`).

### Other Changes

* Improved documentation for some functions.

---

# OpTop 0.4.0

### New Functions

Since we have two more functions, I feel like this deserves a jump in 
package version. 

* `topic_match()`: detect and extract informative and uninformative components.

* `agg_topic_stability()`: implements _Test 4_ from the methodological paper 
[Lewis and Grossetti (2019)].

### Other Changes

* Improved documentation for some functions.

---

# OpTop 0.2.0

### New Functions

* `get_topic_models()`: handy function to immediately get the list of topic models
the user wants to process from a specified environment;

* `topic_stability()`: implements _Tests 2_ and _3_ from the methodological paper 
[Lewis and Grossetti (2019)].

### Other Changes

* Formal declaration of `LDA_VEM` objects as functions input.

* All the functions now have more detailed and better documentations.

* Added Continuous Integration with Travis CI and AppVeyor.

### Bug Fixes

* Choice of quantile algorithms is now fully supported.

---

# OpTop 0.1.0

First version! 

### New Functions

* `word_proportions()`: computes word proportions from a `corpus` object created 
by __quanteda__ [Benoit et al. (2018)];

* `optimal_topic()`: implements _Test 1_ from the methodological paper 
[Lewis and Grossetti (2019)].
