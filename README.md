# Bayesian-GARCH-Modeling-of-Functional-Sports-Data
Implementation and data for the Bayesian Additive model for shot put performance analysis

The folder contains:
* the cleaned and prepared data set
* the R file where the Gibbs sampler is set up and executed: "set_up_Gibbs.r"
* the R file for the posterior analysis: "trajecotry plot.r"
* the C++ code implementing the Gibbs sampler and the Adaptive Metropolis algorithms: "Functional_GARCH_AMH.cpp"

Note that a C++ compiler is required to run the code, despite execution is launched within R.
