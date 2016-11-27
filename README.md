## mf-gp-ucb 
This is a Matlab implementation of the Multi-fidelity Gaussian Process Upper Confidence
Bound algorithm for Bayesian optimisation with multi-fidelity approximations. For more
details please read our paper (below).


### Download
You can download the code from github.
```bash
$ git clone https://github.com/kirthevasank/mf-gp-ucb 
```

### Installation & Getting Started
- Simply execute the following where path-to-mf-gp-ucb is the download directory.
```matlab
>> addpath(genpath('path-to-mf-gp-ucb'))
```
- The demos directory has two examples on how to use this library.
- demo_short.m: Easy set up and uses all default configurations.
- demo_long.m: This demonstrates how to set all parameters individually in case you
  want to customise them. **(TBD)**
- To be able to use the software, you will need to construct an
[mfFunction Object](https://github.com/kirthevasank/mf-gp-ucb/blob/master/mfBO/mfFunction).
This is quite straightforward. See the examples in demos/functions/.
[demos/functions](https://github.com/kirthevasank/mf-gp-ucb/blob/master/demos/functions).


### Some notes
- We choose the GP hyper-parameters every 25 iterations via marginal likelihood
  maximisation for each GP. The chosen values are printed out.
- We report progress on the optimisation every 10 iterations. We report the cost incurred,
  the number of queries at each fidelity and the maximum value found so far.


### Citation
If you use this library in your academic work please cite our
[NIPS 2016 paper](http://www.cs.cmu.edu/~kkandasa/pubs/kandasamyNIPS16mfbo.pdf):
"Gaussian Process Bandit Optimisation with Multi-fidelity Evaluations",
Kirthevasan Kandasamy, Gautam Dasarathy, Junier Oliva, Jeff Schneider, Barnabas Poczos.

We use DiRect to optimise the acquisition function. The implementation was taken from
Dan Finkel (2004).


### License
This software is released under the MIT license. For more details, please refer
[LICENSE.txt](https://github.com/kirthevasank/mf-gp-ucb/blob/master/LICENSE.txt).

"Copyright 2015 Kirthevasan Kandasamy"

- For questions/ bug reports please email kandasamy@cs.cmu.edu

