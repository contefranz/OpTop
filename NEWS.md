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


# OpTop 0.1.0

First version! 

### New Functions

* `word_proportions()`: computes word proportions from a `corpus` object created 
by __quanteda__ [Benoit et al. (2018)];

* `optimal_topic()`: implements _Test 1_ from the methodological paper 
[Lewis and Grossetti (2019)].
