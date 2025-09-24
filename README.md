Code for "Robust domain selection for functional data via interval-wise testing and effect size mapping"

"**Rselection_source_functions.R**"
 : source code to implement robust domain selection, specifically including functions for (i) functional M-estimator using penalized spline functions under Huber loss functions, and (ii) computation of adjusted p-value functions 


"**group_mean_functions.RData**"
 : R data containing $\mu_1(t)$ and $\mu_2(t)$ in sumulation study of Section 3

"**sim_example_t3_sigma3.RData**"
 : one set of generated trajectories from simulation study (t3 process under $\sigma_e=3$ and partial sampling) to illustrate code implementation

"**Sim_generation.R**"
 : simulation data generation code under various types of outliers, levels of error noise, sample size, and sampling structure. See Section 3 for details.


"**Domain_selection_effectsize_map.R**"
 : The first chuck implementing our domain selection algorithm to "sim_example_t3_sigma3.RData". The second chunk is to calculate the effect size and generate the heatmap.
