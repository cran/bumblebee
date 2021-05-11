#' Observed HIV transmission pairs
#' 
#' Counts of directed HIV transmission pairs observed within and between 
#' intervention and control communities in the 30-community BCPP/Ya Tsie 
#' HIV prevention trial in Botswana (2013-2018). The Botswana -Ya Tsie 
#' trial was a pair-matched community randomized trial that evaluated
#' the effect of a universal HIV test and treat intervention in reducing
#' population-level incidence. For further details see references and:
#' \url{https://magosil86.github.io/bumblebee/}.
#' 
#' @format A data frame:
#' \describe{
#' 	\item{H1_group}{Name of population group 1}
#' 	\item{H2_group}{Name of population group 1}
#' 	\item{num_linked_pairs_observed}{Number of observed directed transmission pairs between samples from population groups 1 and 2}}
#' 
#' @references
#' Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#' HIV transmission in the context of a universal testing and treatment
#' trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#'
#' @source \url{https://magosil86.github.io/bumblebee/}
#'
"counts_hiv_transmission_pairs"

#' Sampling fequency
#' 
#' Estimated number of individuals with HIV in intervention and control communities 
#' of the BCPP/Ya Tsie trial, and the number of individuals sampled from each for 
#' HIV viral phylogenetic analysis.
#' 
#' @format A data frame:
#' \describe{
#' \item{population_group}{Population group}
#' \item{number_sampled}{Number of individuals sampled per population group}
#' \item{number_population}{Estimated number of individuals in each population group}
#' }
#' 
#' @references
#' Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#' HIV transmission in the context of a universal testing and treatment
#' trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#'
#' @source \url{https://magosil86.github.io/bumblebee/}
#'
"sampling_frequency"

#' Estimated HIV transmission flows
#' 
#' Estimated HIV transmissions within and betweeen intervention and control communities
#' in the BCPP/Ya Tsie trial population adjusted for variability in sampling.
#' 
#' @format A data frame:
#' \describe{
#'     \item{H1_group}{Name of population group 1} 
#'     \item{H2_group}{Name of population group 2} 
#'     \item{number_hosts_sampled_group_1}{Number of individuals sampled from population group 1}
#'     \item{number_hosts_sampled_group_2}{Number of individuals sampled from population group 2}
#'     \item{number_hosts_population_group_1}{Estimated number of individuals in population group 1}
#'     \item{number_hosts_population_group_2}{Estimated number of individuals in population group 2}
#'     \item{max_possible_pairs_in_sample}{Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2}
#'     \item{max_possible_pairs_in_population}{Number of distinct possible transmission pairs between individuals in population groups 1 and 2}
#'     \item{num_linked_pairs_observed}{Number of observed directed transmission pairs between samples from population groups 1 and 2}
#'     \item{p_hat}{Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked}                                      
#'     \item{est_linkedpairs_in_population}{Estimated transmission pairs between population groups 1 and 2}
#'     \item{theta_hat}{Estimated transmission flows or relative probability of transmission within and between population groups 1 and 2 adjusted
#'           for sampling heterogeneity. More precisely, the conditional probability that a pair of pathogen sequences is from a specific population 
#'           group pairing given that the pair is linked.}
#'     \item{obs_trm_pairs_est_goodman}{Point estimate, Goodman method Confidence intervals for observed transmission pairs}
#'     \item{obs_trm_pairs_lwr_ci_goodman}{Lower bound of Goodman confidence interval} 
#'     \item{obs_trm_pairs_upr_ci_goodman}{Upper bound of Goodman confidence interval} 
#'     \item{est_goodman}{Point estimate, Goodman method Confidence intervals for estimated transmission flows}
#'     \item{lwr_ci_goodman}{Lower bound of Goodman confidence interval} 
#'     \item{upr_ci_goodman}{Upper bound of Goodman confidence interval} 
#'     \item{prob_group_pairing_and_linked}{Probability that a pair of pathogen sequences is from a specific population group pairing and is linked}
#'     \item{c_hat}{Probability that a randomly selected pathogen sequence in one population group links to at least 
#'           one pathogen sequence in another population group i.e. probability of clustering}
#'     \item{est_goodman_cc}{Point estimate, Goodman method Confidence intervals with continuity correction}
#'     \item{lwr_ci_goodman_cc}{Lower bound of Goodman confidence interval} 
#'     \item{upr_ci_goodman_cc}{Upper bound of Goodman confidence interval} 
#'     \item{est_sisonglaz}{Point estimate, Sison-Glaz method Confidence intervals}
#'     \item{lwr_ci_sisonglaz}{Lower bound of Sison-Glaz confidence interval} 
#'     \item{upr_ci_sisonglaz}{Upper bound of Sison-Glaz confidence interval} 
#'     \item{est_qhurst_acswr}{Point estimate, Queensbury-Hurst method Confidence intervals via ACSWR r package} 
#'     \item{lwr_ci_qhurst_acswr}{Lower bound of Queensbury-Hurst confidence interval} 
#'     \item{upr_ci_qhurst_acswr}{Upper bound of Queensbury-Hurst confidence interval} 
#'     \item{est_qhurst_coinmind}{Point estimate, Queensbury-Hurst method Confidence intervals via CoinMinD r package} 
#'     \item{lwr_ci_qhurst_coinmind}{Lower bound of Queensbury-Hurst confidence interval}
#'     \item{upr_ci_qhurst_coinmind}{Upper bound of Queensbury-Hurst confidence interval}
#'     \item{lwr_ci_qhurst_adj_coinmind}{Lower bound of Queensbury-Hurst confidence interval adjusted} 
#'     \item{upr_ci_qhurst_adj_coinmind}{Upper bound of Queensbury-Hurst confidence interval adjusted}
#' }
#' 
#' @references
#' Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#' HIV transmission in the context of a universal testing and treatment
#' trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#'
#' @source \url{https://magosil86.github.io/bumblebee/}
#'
"estimated_hiv_transmission_flows"
