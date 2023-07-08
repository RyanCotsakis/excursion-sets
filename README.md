# excursion-sets
### Functions for Analyzing Binary and Continuous Fields in R

This repository contains a set of functions for analyzing binary and continuous fields in R. The functions can be used for a variety of applications, including image processing, data analysis, and statistical inference.

### Demos
- `demos/consistency_of_perimeter_estimate.R` recreates Figures 7 and 9 in Cotsakis et al. (2022), and shows how to estimate the perimeter of pixelated excursion sets of 2D anisotropic random fields.
- `demos/effect_of_anisotropy_on_perimeter_estimate.R` recreates Figure 6 in Cotsakis et al. (2022)

### References

- Hermine Biermé, H. and Agnès Desolneux. The effect of discretization on the mean geometry of a 2D random field. 2021. Annales Henri Lebesgue **4** 1295–1345
- Ryan Cotsakis, Elena Di Bernardino, and Thomas Opitz. On the perimeter estimation of pixelated excursion sets of 2D anisotropic random fields. 2022. [hal-03582844v2](https://hal.science/hal-03582844v2).
- Ryan Cotsakis. Computable bounds for the reach and r-convexity of subsets of Rd. 2022. [arXiv:2212.01013](https://arxiv.org/abs/2212.01013).

Please cite these references in your work if our functions are of use.

## Installation of RandomFields
If the R package [**RandomFields**](http://cran.nexr.com/web/packages/RandomFields/index.html) is not already installed on your computer, you will be prompted for an installation upon running `demos/utils.R`.
