# backShift
R Code for 'backShift', an algorithm to estimate the connectivity matrix of a directed (possibly cyclic) graph with hidden variables.
The underlying system is required to be linear and we assume that observations under different shift interventions are available. 

## Installation

### From CRAN
```r
install.packages("backShift")
```

### From Github with `devtools`
```r
devtools::install_github("christinaheinze/backShift")
```

## Code
See the [package vignettes](https://github.com/christinaheinze/backShift/tree/master/vignettes) for more details.

## References
D. Rothenhaeusler, C. Heinze, J. Peters and N. Meinshausen. "backShift: Learning causal cyclic graphs from unknown shift interventions". [arXiv](http://arxiv.org/abs/1506.02494).
