# MATLAB_Examples_Samples_UCD
1. As a graduate student at UC Davis, I sat in on a graduate cosmology course taught by
my PhD supervisor in which I wrote matlab code and generated matlab plots as a part of
several homework assignments. Below I will give the name of the matlab code files along with a brief description of what the code does. Further details are given in comment lines in the code itself. I also include and describe any corresponding plots that were generated.

(a) hw1prob1.m -- In this homework exercise I learned how to utilize some basic matlab syntax and features relevant to the use of functions, data structures, and plotting. The data structure "par" contains all the fixed parameters. I also utilized the matlab function subroutines f1.m, f2.m, and fSin.m (all also attached) to define vectors of coordinates. Further details can be found in comment lines in the code. The corresponding plot that was generated is in the file "hw1_prob1.jpg" and is also included.

(b) hw1prob2.m and hw1prob2_loglog.m -- Similar to (a) above except a log-log plot was also
generated and the function subroutines g2.m, g3.m, and g4.m were also used. See comment lines in the code for further details. The corresponding plots that were generated are in the files:
hw1_prob2_linear.pdf, hw1.2_linear.jpg, hw1_prob2_linear.jpg, hw1_prob2_linear_axes.jpg,
and hw1_prob2_log_log.jpg.

(c) hw2prob2.m -- In this homework exercise I used matlab to generate log-log plots of the cosmological mass density, radiation density, and cosmological constant density as a function of scale factor (a) given various (constant) parameters. Also, the function subroutines rhol.m, rhom.m, and rhor.m were used.  The corresponding labeled plot is in the file hw2_prob2.jpg

(d) hw2prob4.m -- Similar to (d) except semi-log plots were generated and some slightly different
cosmological density parameters were plotted, with the use of the function subroutines omegal.m,
omegar.m, and omegam.m. The corresponding plot was generated in the file "hw2_prob4.jpg".

(e) hw3prob4.m and hw3prob5.m -- For these homework problems I generated log-log plots
similar to those discussed in (d) but with different parameterizations and normalizations for
the cosmological densities as a function of scale factor a. See all of the comment lines in these
two matlab programs for all of the details. The different parameterizations for the different 
cosmological densities are implemented via the following function subroutines: 

lomegakminus.m
lomegakplus.m
lomegal.m
lomegam.m
lomegar.m
rhol100.m
rhom100.m
rhor100.m

A sample generated plot is in the file hw3.5.jpg.

2. Under the supervision of my PhD thesis supervisor at UC Davis, I wrote a number of matlab
programs that enabled me to learn and experiment with some of the statistical and numerical
topics and tools I would need to understand and utilize for the MCMC Dark Energy work that
I would later engage in and which was a part of my PhD thesis (and which is described
in my resume/CV.) I wrote these matlab programs so that I could learn more about and better
understand (and also visualize) statistical/numerical topics and concepts that I would need when
running the MCMC simulations, including confidence levels, Fisher Matrices, Jacobian Matrices,
Figures of Merit, uncertainties and standard deviations, eigenvalues and eigenvectors, probability
density functions, the chi-squared statistic, etc.. In the course of this work, I learned more about
matlab programming and how to implement various matlab statistical/numerical tools, functions and features.

(a) I wrote the programs chi_trials_2and3_I.m, chi_trials_2and3_I.m, and chi_trails_2and3_scenario1.m in order to create a toy model that included 3 equations -- 3 parameters and 3 observables, given uncertainties for each observable,
and given values for each parameter. The idea was to be able to calculate chi-squared statistics,
standard deviations, and figures of merit, etc. for the toy models under different assumptions and scenarios (e.g., allowing two parameters to float but keeping the third fixed, or allowing all three parameters to float, etc.). Corresponding labeled plots (e.g., Chi-squared vs. random values drawn from a normal distribution) were also generated and have been included in the repository as well. I've included the following matlab programs that I wrote for your perusal:

chi_trails_2and3_scenario1.m
chi_trials_2and3_I.m
chi_trials_2and3_Ia.m
chi_trails_2_3_4.m
chi_trails_2and3.m
chitrials_2and3_plots.m

A sampling of some of the corresponding plots that were generated by these matlab programs are:

plot_as_bs_cs_CLA_chisqA.jpg
plot_chi2Bred-chi2Ared_vs_chi2Bred_abc_var.jpg
plot_chisqB-chisqA_chisqA_abc_var.jpg
plot_cs_CLB-CLA_abc_var.jpg
plot_cs_var_chisqA.jpg
plot3_CLB-CLA_CLA_CLB_abc_var.jpg
plot3_CL_B_bs_cs.jpg
plot3_chisqA_bs_cs.jpg
plot3_chisqB-chisqA_chisqA_chisqB_abc_var.jpg
plot3_chisqBred-chisqAred_chisqAred_chisqBred_abc_var.jpg
