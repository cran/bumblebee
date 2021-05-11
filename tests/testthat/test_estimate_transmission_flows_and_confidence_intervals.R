### library(bumblebee); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also:
###     https://magosil86.github.io/bumblebee/
###
### Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
### HIV transmission in the context of a universal testing and treatment
### trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.


### context("Checking bumblebee example: Estimate transmission flows and confidence intervals ")

test_that("estimates of transmission flows and confidence intervals produced by the estimate_transmission_flows_ci function are correct.", {

    ### Test may take longer than 1 min hence test will be run locally rather than on CRAN
	skip_on_cran()

	### Load libraries
	library(dplyr)


    ### Load data of HIV transmissions within and between intervention and control
    ### communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the 
    ### data ?counts_hiv_transmission_pairs, ?sampling_frequency and ?estimated_hiv_transmission_flows
 
	### The input data comprises counts of observed directed HIV transmission pairs between  
	### individuals sampled from intervention and control communities (i.e. num_linked_pairs_observed); 
	### sampling information; and the estimated HIV transmissions within and between intervention and control 
	### communities in the BCPP/Ya Tsie trial population adjusted for sampling heterogneity
	### (i.e. \code{est_linkedpairs_in_population}).

    data(list = c("counts_hiv_transmission_pairs", "sampling_frequency", "estimated_hiv_transmission_flows"), package = "bumblebee")


	### Estimate transmission flows within and between intervention and control communities
	### accounting for variable sampling among population groups. 
	
	### Note: detailed_report option set to get confidence intervals for the Goodman, Sison-Glaz and Queensbury-Hurst methods

    results_estimate_transmission_flows_and_ci_detailed <- estimate_transmission_flows_and_ci(group_in = sampling_frequency$population_group, 
 																							  individuals_sampled_in = sampling_frequency$number_sampled, 
 																							  individuals_population_in = sampling_frequency$number_population, 
 																							  linkage_counts_in = counts_hiv_transmission_pairs, 
 																							  detailed_report = TRUE, 
 																							  verbose_output = TRUE)

	### Retrieve dataset of estimated transmission flows 
	df_estimated_transmission_flows <- results_estimate_transmission_flows_and_ci_detailed$flows_dataset

	### checking p_hat values, probability of linkage between individuals randomly sampled from their respective population groups
    expect_equal(signif(df_estimated_transmission_flows$p_hat, 4), signif(c(1.269065e-05, 2.662174e-05, 3.674568e-06, 1.130636e-06), 4))
    
    ### checking prob_group_pairing_and_linked values, joint probability of linkage
    expect_equal(signif(df_estimated_transmission_flows$prob_group_pairing_and_linked, 4), signif(c(3.427258e-06, 6.141503e-06, 1.834503e-06, 5.644626e-07), 4))
    
	### checking theta_hat values, conditional probability of linkage or estimated transmission flows within and between population groups
	### adjusted for sampling heterogeneity
    expect_equal(signif(df_estimated_transmission_flows$theta_hat, 4), signif(c(0.2863750, 0.5131720, 0.1532875, 0.0471654), 4))
    
    ### checking c_hat values, probability of clustering
    expect_equal(signif(df_estimated_transmission_flows$c_hat, 4), signif(c(0.16557535, 0.29598598, 0.05106208, 0.01479457), 4))
 
    ### checking Goodman confidence intervals without continuity correction 
    expect_equal(signif(df_estimated_transmission_flows$upr_ci_goodman, 4), signif(c(0.4395994, 0.6594660, 0.2939222, 0.1604764), 4))

    expect_equal(signif(df_estimated_transmission_flows$lwr_ci_goodman, 4), signif(c(0.17032542, 0.36458588, 0.07298735, 0.01265617), 4))
   
	### checking Goodman confidence intervals with continuity correction 
    expect_equal(signif(df_estimated_transmission_flows$lwr_ci_goodman_cc, 4), signif(c(0.16554020, 0.35898561, 0.06896478, 0.01013982), 4))

    expect_equal(signif(df_estimated_transmission_flows$upr_ci_goodman_cc, 4), signif(c(0.4459219, 0.6649797, 0.3009860, 0.1689346), 4))


	### Results should also resemble S10_Table in BCPP deep-sequence phylogenetics paper
	### S10 Table. Observed and predicted proportions of HIV-1 transmissions in the trial population
	###            by randomization-intervention condition


})

