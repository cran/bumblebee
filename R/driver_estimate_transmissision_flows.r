# Goal: Generate and export S3 class/methods to estimate transmission flows 
#       and document using roxygen tags


#' @title \code{estimate_transmission_flows_and_ci} Estimates transmission flows and corresponding confidence intervals
#' 
#' @description This function estimates transmission flows or the relative probability
#'              of transmission within and between population groups accounting for 
#'              variable sampling among population groups. 
#' 
#' Corresponding confidence intervals are provided with the following methods: Goodman,
#' Goodman with a continuity correction, Sison-Glaz and Queensbury-Hurst.
#'              
#' @aliases estimate_transmission_flows_and_ci flows_and_ci 
#'              
#' @param group_in, A character vector indicating population groups/strata (e.g. communities, age-groups, genders or trial arms) between which transmission flows will be evaluated, 
#' @param individuals_sampled_in, A numeric vector indicating the number of individuals sampled per population group, 
#' @param individuals_population_in, A numeric vector of the estimated number of individuals per population group, 
#' @param linkage_counts_in A data.frame of counts of linked pairs identified between samples of each population
#'        group pairing of interest. \cr
#'        The data.frame should contain the following three fields: 
#'  	   \itemize{
#'  	       \item H1_group (character) Name of population group 1 
#'  	       \item H2_group (character) Name of population group 2 
#'  	       \item number_linked_pairs_observed (numeric) Number of 
#'  	             observed directed transmission pairs between samples 
#'  	             from population groups 1 and 2} 
#' @param verbose_output, A boolean value to display intermediate output (Default is \code{FALSE})                       
#' @param detailed_report, A boolean value to produce detailed output of the analysis                       
#' @param \dots Further arguments.
#' 
#' @return Returns a data.frame containing:
#' \itemize{
#' 	    \item H1_group, Name of population group 1 
#' 	    \item H2_group, Name of population group 2 
#' 	    \item number_hosts_sampled_group_1, Number of individuals sampled from population group 1
#' 	    \item number_hosts_sampled_group_2, Number of individuals sampled from population group 2
#' 	    \item number_hosts_population_group_1, Estimated number of individuals in population group 1
#' 	    \item number_hosts_population_group_2, Estimated number of individuals in population group 2
#' 	    \item max_possible_pairs_in_sample, Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2
#' 	    \item max_possible_pairs_in_population, Number of distinct possible transmission pairs between individuals in population groups 1 and 2
#' 	    \item num_linked_pairs_observed, Number of observed directed transmission pairs between samples from population groups 1 and 2 
#' 	    \item p_hat, Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked                                      
#' 	    \item est_linkedpairs_in_population, Estimated transmission pairs between population groups 1 and 2
#' 	    \item theta_hat, Estimated transmission flows or relative probability of transmission within and between population groups 1 and 2 adjusted
#' 	          for sampling heterogeneity. More precisely, the conditional probability that a pair of pathogen sequences is from a specific population 
#' 	          group pairing given that the pair is linked.
#' 	    \item obs_trm_pairs_est_goodman, Point estimate, Goodman method Confidence intervals for observed transmission pairs
#' 	    \item obs_trm_pairs_lwr_ci_goodman, Lower bound of Goodman confidence interval 
#' 	    \item obs_trm_pairs_upr_ci_goodman, Upper bound of Goodman confidence interval 
#' 	    \item est_goodman, Point estimate, Goodman method Confidence intervals for estimated transmission flows
#' 	    \item lwr_ci_goodman, Lower bound of Goodman confidence interval 
#' 	    \item upr_ci_goodman, Upper bound of Goodman confidence interval 
#' }
#' 
#' The following additional fields are returned if the detailed_report flag is set
#' 
#' \itemize{
#' 	    \item prob_group_pairing_and_linked, Probability that a pair of pathogen sequences is from a specific population group pairing and is linked
#' 	    \item c_hat, Probability that a randomly selected pathogen sequence in one population group links to at least 
#' 	          one pathogen sequence in another population group i.e. probability of clustering
#' 	    \item est_goodman_cc, Point estimate, Goodman method Confidence intervals with continuity correction
#' 	    \item lwr_ci_goodman_cc, Lower bound of Goodman confidence interval 
#' 	    \item upr_ci_goodman_cc, Upper bound of Goodman confidence interval 
#' 	    \item est_sisonglaz, Point estimate, Sison-Glaz method Confidence intervals
#' 	    \item lwr_ci_sisonglaz, Lower bound of Sison-Glaz confidence interval 
#' 	    \item upr_ci_sisonglaz, Upper bound of Sison-Glaz confidence interval 
#' 	    \item est_qhurst_acswr, Point estimate, Queensbury-Hurst method Confidence intervals via ACSWR r package 
#' 	    \item lwr_ci_qhurst_acswr, Lower bound of Queensbury-Hurst confidence interval 
#' 	    \item upr_ci_qhurst_acswr, Upper bound of Queensbury-Hurst confidence interval 
#' 	    \item est_qhurst_coinmind, Point estimate, Queensbury-Hurst method Confidence intervals via CoinMinD r package 
#' 	    \item lwr_ci_qhurst_coinmind, Lower bound of Queensbury-Hurst confidence interval
#' 	    \item upr_ci_qhurst_coinmind, Upper bound of Queensbury-Hurst confidence interval
#' 	    \item lwr_ci_qhurst_adj_coinmind, Lower bound of Queensbury-Hurst confidence interval adjusted 
#' 	    \item upr_ci_qhurst_adj_coinmind, Upper bound of Queensbury-Hurst confidence interval adjusted
#' }
#' 
#' @details Counts of observed directed transmission pairs can be obtained from 
#'     deep-sequence phylogenetic data (via phyloscanner) or from known
#'     epidemiological contacts. Note: Deep-sequence data is also commonly 
#'     referred to as high-throughput or next-generation sequence data. See
#'     references to learn more about phyloscanner. 
#'
#' The \code{estimate_transmission_flows_and_ci()} function is a 
#' wrapper function that calls the following functions:
#'
#' \enumerate{
#' 
#'    \item The \code{prep_p_hat()} function to determine all possible 
#'          combinations of the population groups/strata provided by 
#'          the user. Type \code{?prep_p_hat()} at R prompt to learn 
#'          more.
#'     	  
#'    \item The \code{estimate_p_hat()} function to compute the 
#'          probability of linkage between pathogen sequences from 
#'          two individuals randomly sampled from their respective 
#'          population groups. Type \code{?estimate_p_hat()} at R 
#'          prompt to learn more.
#'    	  
#'    \item The \code{estimate_theta_hat()} function that uses 
#'          \code{p_hat} estimates to compute the conditional 
#'          probability of linkage that a pair of pathogen sequences 
#'          is from a specific population group pairing given that 
#'          the pair is linked. The conditional  probability, 
#'          \code{theta_hat} represents transmission flows or 
#'          the relative probability of transmission within and between 
#'          population groups adjusted for variable sampling among 
#'          population groups. Type \code{?estimate_theta_hat()} at R 
#'          prompt to learn more.
#'    	  
#'    \item The \code{estimate_multinom_ci()} function to estimate 
#'          corresponding confidence intervals for the computed
#'          transmission flows.
#' }
#' 
#' Further to estimating transmission flows and corresponding confidence
#' intervals the \code{estimate_transmission_flows_and_ci()} function provides 
#' estimates for:
#'          
#' \enumerate{
#'     
#'    \item \code{prob_group_pairing_and_linked}, the joint probability that a
#'          pair of pathogen sequences is from a specific population group
#'          pairing and linked. Type \code{?estimate_prob_group_pairing_and_linked()} 
#'          at R prompt to learn more.
#'               
#'    \item \code{c_hat}, the probability of clustering that a pathogen sequence 
#'          from a population group of interest is linked to one or more 
#'          pathogen sequences in another population group of interest. Type
#'          \code{?estimate_c_hat()} at R prompt to learn more.
#'      
#' }
#'         
#' @seealso \code{\link{estimate_theta_hat}} and \code{\link{estimate_multinom_ci}} to learn
#'          more about estimation of transmission flows and confidence intervals.
#' 
#' @references
#' \enumerate{
#'
#'    \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'          HIV transmission in the context of a universal testing and treatment
#'          trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#'    \item Carnegie, N.B., et al., Linkage of viral sequences among HIV-infected
#' 		    village residents in Botswana: estimation of linkage rates in the 
#' 		    presence of missing data. PLoS Computational Biology, 2014. 10(1): 
#' 		    p. e1003430.
#'
#'    \item Cherry, S., A Comparison of Confidence Interval Methods for Habitat 
#'          Use-Availability Studies. The Journal of Wildlife Management, 1996. 
#'          60(3): p. 653-658.
#' 
#'    \item Ratmann, O., et al., Inferring HIV-1 transmission networks and sources of 
#' 		    epidemic spread in Africa with deep-sequence phylogenetic analysis. 
#' 		    Nature Communications, 2019. 10(1): p. 1411.
#' 
#'    \item Wymant, C., et al., PHYLOSCANNER: Inferring Transmission from Within- and
#' 		    Between-Host Pathogen Genetic Diversity. Molecular Biology and Evolution,
#' 		    2017. 35(3): p. 719-733.
#'  
#'    \item Goodman, L. A. On Simultaneous Confidence Intervals for Multinomial Proportions 
#' 		    Technometrics, 1965. 7, 247-254.
#' 
#'    \item Sison, C.P and Glaz, J. Simultaneous confidence intervals and sample size determination
#' 		    for multinomial proportions. Journal of the American Statistical Association, 
#' 		    1995. 90:366-369.
#' 
#'    \item Glaz, J., Sison, C.P. Simultaneous confidence intervals for multinomial proportions. 
#' 		    Journal of Statistical Planning and Inference, 1999. 82:251-262.
#' 
#'    \item May, W.L., Johnson, W.D. Constructing two-sided simultaneous confidence intervals for 
#' 		    multinomial proportions for small counts in a large number of cells. 
#' 		    Journal of Statistical Software, 2000. 5(6).
#' 		    Paper and code available at https://www.jstatsoft.org/v05/i06.
#' 
#' }
#' 
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#' 
#' # Estimate transmission flows and confidence intervals
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs and ?sampling_frequency 
#' 
#' # View counts of observed directed HIV transmissions within and between intervention 
#' # and control communities
#' counts_hiv_transmission_pairs
#' 
#' # View the estimated number of individuals with HIV in intervention and control 
#' # communities and the number of individuals sampled from each
#' sampling_frequency
#' 
#' # Estimate transmission flows within and between intervention and control communities
#' # accounting for variable sampling among population groups. 
#' 
#' # Basic output
#' results_estimate_transmission_flows_and_ci <- estimate_transmission_flows_and_ci(
#'     group_in = sampling_frequency$population_group, 
#'     individuals_sampled_in = sampling_frequency$number_sampled, 
#'     individuals_population_in = sampling_frequency$number_population, 
#'     linkage_counts_in = counts_hiv_transmission_pairs)
#' 
#' # View results
#' results_estimate_transmission_flows_and_ci
#'
#' # Retrieve dataset of estimated transmission flows 
#' dframe <- results_estimate_transmission_flows_and_ci$flows_dataset
#' 
#' # Detailed output
#' results_estimate_transmission_flows_and_ci_detailed <- estimate_transmission_flows_and_ci(
#'     group_in = sampling_frequency$population_group, 
#' 	   individuals_sampled_in = sampling_frequency$number_sampled, 
#' 	   individuals_population_in = sampling_frequency$number_population, 
#' 	   linkage_counts_in = counts_hiv_transmission_pairs, 
#' 	   detailed_report = TRUE)
#' 
#' # View results
#' results_estimate_transmission_flows_and_ci_detailed
#'
#' # Retrieve dataset of estimated transmission flows 
#' dframe <- results_estimate_transmission_flows_and_ci_detailed$flows_dataset 
#' 
#' # Options:
#' # To show intermediate output set verbose_output = TRUE
#' 
#' # Basic output
#' results_estimate_transmission_flows_and_ci <- estimate_transmission_flows_and_ci(
#'     group_in = sampling_frequency$population_group, 
#' 	   individuals_sampled_in = sampling_frequency$number_sampled, 
#' 	   individuals_population_in = sampling_frequency$number_population, 
#' 	   linkage_counts_in = counts_hiv_transmission_pairs, 
#' 	   verbose_output = TRUE)
#' 
#' # View results
#' results_estimate_transmission_flows_and_ci
#' 
#' # Detailed output
#' results_estimate_transmission_flows_and_ci_detailed <- estimate_transmission_flows_and_ci(
#'     group_in = sampling_frequency$population_group, 
#' 	   individuals_sampled_in = sampling_frequency$number_sampled, 
#' 	   individuals_population_in = sampling_frequency$number_population, 
#' 	   linkage_counts_in = counts_hiv_transmission_pairs, 
#' 	   detailed_report = TRUE, 
#' 	   verbose_output = TRUE)
#' 
#' # View results
#' results_estimate_transmission_flows_and_ci_detailed
#' 
#' @export
estimate_transmission_flows_and_ci <- function(group_in, individuals_sampled_in, individuals_population_in, linkage_counts_in, ...) UseMethod("estimate_transmission_flows_and_ci")

#' @describeIn estimate_transmission_flows_and_ci Estimates transmission flows and accompanying confidence intervals

#' @export
estimate_transmission_flows_and_ci.default <- function(group_in, 
                                                       individuals_sampled_in, 
                                                       individuals_population_in, 
                                                       linkage_counts_in,
                                                       detailed_report = FALSE,
                                                       verbose_output = FALSE,
                                                       ...) {

	results_estimate_transmission_flows <- estimate_transmission_flows(group_in, individuals_sampled_in, individuals_population_in, linkage_counts_in, detailed_report, verbose_output)

    results_estimate_transmission_flows$call <- match.call()
    
    class(results_estimate_transmission_flows) <- "estimate_transmission_flows_and_ci"
    
              
    results_estimate_transmission_flows
                                               
}


#' @title \code{prep_p_hat} Prepares input data to estimate \code{p_hat}
#' 
#' @description This function generates variables required for estimating
#'              \code{p_hat}, the probability that pathogen sequences from
#'              two individuals randomly sampled from their respective
#'              population groups are linked. For a population group 
#'              pairing \eqn{(u,v)}, \code{prep_p_hat} determines all possible 
#'              group pairings i.e. \eqn{(uu, uv, vu, vv)}.              
#' 
#' @aliases pp prep_phat prep_p_hat 
#' 
#' @param group_in A character vector indicating population groups/strata 
#'        (e.g. communities, age-groups, genders or trial arms) between which 
#'        transmission flows will be evaluated, 
#' @param individuals_sampled_in A numeric vector indicating the number of 
#'        individuals sampled per population group, 
#' @param individuals_population_in A numeric vector of the estimated number 
#'        of individuals per population group, 
#' @param linkage_counts_in A data.frame of counts of linked pairs identified 
#'        between samples of each population group pairing of interest. \cr
#'        The data.frame should contain the following three fields: 
#' 	   \itemize{
#' 	       \item H1_group (character) Name of population group 1 
#' 	       \item H2_group (character) Name of population group 2 
#' 	       \item number_linked_pairs_observed (numeric) Number of 
#' 	             observed directed transmission pairs between samples 
#' 	             from population groups 1 and 2
#' 	   } 
#' @param verbose_output A boolean value to display intermediate output.
#'        (Default is \code{FALSE})                        
#' @param \dots Further arguments.
#' 
#' @return Returns a data.frame containing:
#' \itemize{
#' \item H1_group, Name of population group 1 
#' \item H2_group, Name of population group 2 
#' \item number_hosts_sampled_group_1, Number of individuals sampled 
#'       from population group 1
#' \item number_hosts_sampled_group_2, Number of individuals sampled 
#'       from population group 2
#' \item number_hosts_population_group_1, Estimated number of individuals 
#'       in population group 1
#' \item number_hosts_population_group_2, Estimated number of individuals 
#'       in population group 2
#' \item max_possible_pairs_in_sample, Number of distinct possible transmission 
#'       pairs between individuals sampled from population groups 1 and 2
#' \item max_possible_pairs_in_population, Number of distinct possible 
#'       transmission pairs between individuals in population groups 1 and 2
#' \item num_linked_pairs_observed, Number of observed directed transmission 
#'       pairs between samples from population groups 1 and 2
#' }
#' 
#' @seealso \code{\link{estimate_p_hat}}
#' 
#' @details Counts of observed directed transmission pairs can be obtained from 
#'          deep-sequence phylogenetic data (via phyloscanner) or from known
#'          epidemiological contacts. Note: Deep-sequence data is also commonly 
#'          referred to as high-throughput or next-generation sequence data. See
#'          references to learn more about phyloscanner.
#' 
#' @references
#' \enumerate{
#'
#' \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'       HIV transmission in the context of a universal testing and treatment
#'       trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#' \item Ratmann, O., et al., Inferring HIV-1 transmission networks and 
#'       sources of epidemic spread in Africa with deep-sequence phylogenetic 
#'       analysis. Nature Communications, 2019. 10(1): p. 1411.
#' 
#' \item Wymant, C., et al., PHYLOSCANNER: Inferring Transmission from Within 
#'       and Between-Host Pathogen Genetic Diversity. Molecular Biology and 
#'       Evolution, 2017. 35(3): p. 719-733.
#'   
#' }
#' 
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#'
#' # Prepare input to estimate p_hat
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs and ?sampling_frequency 
#' 
#' # View counts of observed directed HIV transmissions within and between intervention 
#' # and control communities
#' counts_hiv_transmission_pairs
#'
#' # View the estimated number of individuals with HIV in intervention and control 
#' # communities and the number of individuals sampled from each
#' sampling_frequency
#'
#' results_prep_p_hat <- prep_p_hat(group_in = sampling_frequency$population_group, 
#'                                  individuals_sampled_in = sampling_frequency$number_sampled, 
#'                                  individuals_population_in = sampling_frequency$number_population, 
#'                                  linkage_counts_in = counts_hiv_transmission_pairs, 
#'                                  verbose_output = TRUE)
#'
#' # View results
#' results_prep_p_hat
#' 
#' @export
prep_p_hat <- function(group_in, individuals_sampled_in, individuals_population_in, linkage_counts_in, ...) UseMethod("prep_p_hat")

#' @describeIn prep_p_hat Prepares input data to estimate \code{p_hat}
#' @export
prep_p_hat.default <- function(group_in, 
							   individuals_sampled_in, 
							   individuals_population_in, 
							   linkage_counts_in,
							   verbose_output = FALSE, ...) {

    # Check whether all required variables are present
    
    if (missing(group_in))
        stop("Group value missing, need to specify a character vector of population groups/strata between which to estimate transmission flows.")
 
    if (missing(individuals_sampled_in))
        stop("Individuals sampled missing, need to specify a numeric vector of the number of individuals sampled per population group.")

    if (missing(individuals_population_in))
        stop("Individuals in population missing, need to specify a numeric vector of the estimated number of individuals per population group.")

    if (missing(linkage_counts_in))
        stop("Identified transmission pairs missing, need to specify a data.frame containing the number of directed transmission pairs identified/observed in the sample per population group pairing.")


    # Verify data types of required variables
    if (!is.numeric(c(individuals_sampled_in, individuals_population_in))) {

        stop("individuals_sampled_in and individuals_population_in should be of type, numeric.")

    }


    if (!is.character(c(group_in))) {

        stop("group_in should be of type, character.")

    }


    if (!(ncol(linkage_counts_in) == 3)) {

        stop("linkage_counts_in should be a data.frame with three fields: H1_counts, H2_counts and number_linked_pairs_observed.")
    
    }
    
	results_prepare_input_for_get_p_hat <- prepare_input_for_get_p_hat(group_in, 
																	   individuals_sampled_in, 
																	   individuals_population_in, 
																	   linkage_counts_in,
																	   verbose_output)
          
    results_prepare_input_for_get_p_hat
    

}


      
#' @title \code{estimate_p_hat} Estimates probability of linkage between two individuals
#' 
#' @description This function computes the probability that pathogen sequences from 
#'              two individuals randomly sampled from their respective population groups 
#'              (e.g. communities) are linked.
#' 
#' @details For a population group pairing \eqn{(u,v)}, \code{p_hat} is computed as the
#' 		 fraction of distinct possible pairs between samples from groups \eqn{u} and 
#' 		 \eqn{v} that are linked. Note: The number of distinct possible \eqn{(u,v)} 
#' 		 pairs in the sample is the product of sampled individuals in groups \eqn{u} 
#' 		 and \eqn{u}. If \eqn{u = v}, then the distinct possible pairs is the number 
#' 		 of individuals sampled in population group \eqn{u} choose 2. See bumblebee
#' 		 website for more details \url{https://magosil86.github.io/bumblebee/}.
#'          
#' @aliases p_hat phat est_p_hat estim_p_hat estimate_p_hat
#' 
#' @param df_counts A data.frame returned by the function: [prep_p_hat()]
#' @param \dots Further arguments.
#' 
#' @return Returns a data.frame containing:
#' \itemize{
#'     \item H1_group, Name of population group 1 
#'     \item H2_group, Name of population group 2 
#' 	   \item number_hosts_sampled_group_1, Number of individuals sampled from population group 1
#' 	   \item number_hosts_sampled_group_2, Number of individuals sampled from population group 2
#' 	   \item number_hosts_population_group_1, Estimated number of individuals in population group 1
#' 	   \item number_hosts_population_group_2, Estimated number of individuals in population group 2
#' 	   \item max_possible_pairs_in_sample, Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2
#' 	   \item max_possible_pairs_in_population, Number of distinct possible transmission pairs between individuals in population groups 1 and 2
#' 	   \item num_linked_pairs_observed, Number of observed directed transmission pairs between samples from population groups 1 and 2
#' 	   \item p_hat, Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked                                      
#' }
#' 
#' @seealso See \code{\link{prep_p_hat}} to prepare input data to estimate \code{p_hat} 
#' 
#' @references
#' \enumerate{
#'
#' \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'       HIV transmission in the context of a universal testing and treatment
#'       trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#' \item Carnegie, N.B., et al., Linkage of viral sequences among HIV-infected
#' 	     village residents in Botswana: estimation of linkage rates in the 
#' 	     presence of missing data. PLoS Computational Biology, 2014. 10(1): 
#' 	     p. e1003430.
#'
#'}
#'
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#' 
#' # Estimate the probability of linkage between two individuals randomly sampled from
#' # two population groups of interest.
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs and ?sampling_frequency
#' 
#' # Prepare input to estimate p_hat
#' 
#' # View counts of observed directed HIV transmissions within and between intervention 
#' # and control communities
#' counts_hiv_transmission_pairs
#'
#' # View the estimated number of individuals with HIV in intervention and control 
#' # communities and the number of individuals sampled from each
#' sampling_frequency
#'
#' results_prep_p_hat <- prep_p_hat(group_in = sampling_frequency$population_group, 
#'                                  individuals_sampled_in = sampling_frequency$number_sampled, 
#'                                  individuals_population_in = sampling_frequency$number_population, 
#'                                  linkage_counts_in = counts_hiv_transmission_pairs, 
#'                                  verbose_output = FALSE)
#'
#' # View results
#' results_prep_p_hat
#'
#' # Estimate p_hat
#' results_estimate_p_hat <- estimate_p_hat(df_counts = results_prep_p_hat)
#'
#' # View results 
#' results_estimate_p_hat
#' 
#' 
#' @export
estimate_p_hat <- function(df_counts, ...) UseMethod("estimate_p_hat")

#' @describeIn estimate_p_hat Estimates probability of linkage between two individuals
#' @export
estimate_p_hat.default <- function(df_counts, ...) {

	# Check whether the input is a data.frame
	if (!is.data.frame(df_counts)) {
		
		stop("Check input! Expecting a data.frame produced by function: prep_p_hat(), type ?prep_p_hat() for help")

	}

	# Check whether all required variables are present
	
	provided_field_names <- names(df_counts)
	expected_field_names <- c("H1_group", "H2_group", 
	                          "number_hosts_sampled_group_1", "number_hosts_sampled_group_2", 
	                          "number_hosts_population_group_1", "number_hosts_population_group_2", 
	                          "max_possible_pairs_in_sample", "max_possible_pairs_in_population", 
	                          "num_linked_pairs_observed")   

	if (base::isFALSE(base::all(expected_field_names %in% provided_field_names))) {

	    message <- list('Check input! One or more required fields is missing. \nExpecting a data.frame produced by function: prep_p_hat(), type ?prep_p_hat() for help')
	    stop(message)

	}
	
	results_get_p_hat <- get_p_hat(df_counts)
	
	results_get_p_hat
      
}


#' @title \code{estimate_prob_group_pairing_and_linked} Estimates joint probability of linkage
#' 
#' @description This function computes the joint probability that a pair 
#'              of pathogen sequences is from a specific population group
#'              pairing and linked.
#' 
#' @details For a population group pairing \eqn{(u,v)}, the joint probability that a pair
#'          is from groups \eqn{(u,v)} and is linked is computed as
#' 
#' \deqn{(N_uv / N_choose_2) * p_hat_uv ,}
#'
#' where,
#' \itemize{
#'     \item N_uv = N_u * N_v: maximum distinct possible \eqn{(u,v)} pairs in population
#'     \item p_hat_uv: probability of linkage between two individuals randomly sampled 
#'       from groups \eqn{u} and \eqn{v}
#'     \item N choose 2 or (N * (N - 1))/2 : all distinct possible pairs in population.
#' }
#'
#' See bumblebee website for more details \url{https://magosil86.github.io/bumblebee/}.
#' 
#' @aliases prob_group_pairing_and_linked estimate_prob_group_pairing_and_linked
#' 
#' @param df_counts_and_p_hat A data.frame returned by function: [estimate_p_hat()]
#' @param individuals_population_in A numeric vector of the estimated number of individuals per population group
#' @param verbose_output A boolean value to display intermediate output.
#'        (Default is \code{FALSE})                        
#' @param \dots Further arguments.
#' 
#' @return Returns a data.frame containing:
#' \itemize{
#'     \item H1_group, Name of population group 1 
#' 	   \item H2_group, Name of population group 2 
#' 	   \item number_hosts_sampled_group_1, Number of individuals sampled from population group 1
#' 	   \item number_hosts_sampled_group_2, Number of individuals sampled from population group 2
#' 	   \item number_hosts_population_group_1, Estimated number of individuals in population group 1
#' 	   \item number_hosts_population_group_2, Estimated number of individuals in population group 2
#' 	   \item max_possible_pairs_in_sample, Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2
#' 	   \item max_possible_pairs_in_population, Number of distinct possible transmission pairs between individuals in population groups 1 and 2
#' 	   \item num_linked_pairs_observed, Number of observed directed transmission pairs between samples from population groups 1 and 2
#' 	   \item p_hat, Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked                                      
#' 	   \item prob_group_pairing_and_linked, Probability that a pair of pathogen sequences is from a specific population group pairing and is linked
#' }
#' 
#' @seealso See \code{\link{estimate_p_hat}} to prepare input data to estimate \code{prob_group_pairing_and_linked}
#' 
#' @references
#' \enumerate{
#'
#' \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'       HIV transmission in the context of a universal testing and treatment
#'       trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#' \item Carnegie, N.B., et al., Linkage of viral sequences among HIV-infected
#' 	     village residents in Botswana: estimation of linkage rates in the 
#' 	     presence of missing data. PLoS Computational Biology, 2014. 10(1): 
#' 	     p. e1003430.
#'
#'}
#'
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#' 
#' # Estimate joint probability that a pair is from a specific group pairing and linked
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs and ?sampling_frequency
#'
#' 
#' # Load and view data
#' #
#' # The input data comprises counts of observed directed HIV transmission pairs 
#' # within and between intervention and control communities in the BCPP/Ya Tsie 
#' # trial, sampling information and the probability of linkage between individuals
#' # sampled from intervention and control communities (i.e. \code{p_hat})
#' #
#' # See ?estimate_p_hat() for details on estimating p_hat
#' results_estimate_p_hat <- estimated_hiv_transmission_flows[, c(1:10)]
#' 
#' results_estimate_p_hat
#' 
#' # Estimate prob_group_pairing_and_linked
#' results_prob_group_pairing_and_linked <- estimate_prob_group_pairing_and_linked(
#'     df_counts_and_p_hat = results_estimate_p_hat, 
#'     individuals_population_in = sampling_frequency$number_population)
#' 
#' # View results
#' results_prob_group_pairing_and_linked
#'
#' @export
estimate_prob_group_pairing_and_linked <- function(df_counts_and_p_hat, individuals_population_in, ...) UseMethod("estimate_prob_group_pairing_and_linked")
 
#' @describeIn estimate_prob_group_pairing_and_linked Estimates joint probability of linkage 
#'
#' @export
estimate_prob_group_pairing_and_linked.default <- function(df_counts_and_p_hat, 
                                                           individuals_population_in,
                                                           verbose_output = FALSE, ...) { 
                                                          
	# Check whether the input is a data.frame
	if (!is.data.frame(df_counts_and_p_hat)) {
		
		stop("Check input! df_counts_and_p_hat should be a data.frame produced by function: estimate_p_hat(), type ?estimate_p_hat() for help")

	}


	# Check whether all required variables are present
	
	provided_field_names <- names(df_counts_and_p_hat)
	expected_field_names <- c("H1_group", "H2_group", "number_hosts_sampled_group_1", 
	                          "number_hosts_sampled_group_2", "number_hosts_population_group_1", 
	                          "number_hosts_population_group_2", "max_possible_pairs_in_sample", 
	                          "max_possible_pairs_in_population", "num_linked_pairs_observed",
	                          "p_hat")   

	if (base::isFALSE(base::all(expected_field_names %in% provided_field_names))) {

	    message <- list('Check input! One or more required fields is missing. \nExpecting a data.frame produced by function: estimate_p_hat(), type ?estimate_p_hat() for help')
	    stop(message)

	}

    if (missing(individuals_population_in))
        stop("Individuals in population missing, need to specify a numeric vector of the estimated number of individuals per population group.")


    results_get_prob_group_pairing_and_linked <- get_prob_group_pairing_and_linked(df_counts_and_p_hat, individuals_population_in, verbose_output)

    results_get_prob_group_pairing_and_linked
    
}


#' @title \code{estimate_theta_hat} Estimates conditional probability of linkage (transmission flows)
#'  
#' @description This function estimates \code{theta_hat}, the relative probability of 
#'              transmission within and between population groups accounting for variable
#'              sampling rates among population groups. This relative probability is also
#'              refferred to as transmission flows.
#' 
#' @details For a population group pairing \eqn{(u,v)}, the estimated transmission flows
#'          within and between population groups \eqn{u} and \eqn{v}, are represented by
#'          the vector theta_hat,
#'
#' \deqn{\hat{\theta} = ( \hat{\theta}_{uu}, \hat{\theta}_{uv}, \hat{\theta}_{vu}, \hat{\theta}_{vv} ) ,}        
#' 
#' and are computed as \deqn{\hat{\theta_{ij}} = Pr(pair from groups (i,j) | pair is linked), where i = u,v and j = u,v ,}
#'                         
#' \deqn{\hat{\theta_{ij}} = \frac{N_{ij}p_{ij}}{ \sum_m \sum_{n \ge m}N_{mn}p_{mn}}, where i = u,v and j = u,v ,}
#' 
#' See bumblebee website for more details \url{https://magosil86.github.io/bumblebee/}.
#' 
#' @aliases theta_hat est_theta_hat estim_theta_hat estimate_theta_hat
#' 
#' @param df_counts_and_p_hat  A data.frame returned by the function: [estimate_p_hat()]
#' @param \dots Further arguments.
#'
#' @return Returns a data.frame containing:
#' \itemize{
#' 	    \item H1_group, Name of population group 1 
#' 	    \item H2_group, Name of population group 2 
#' 	    \item number_hosts_sampled_group_1, Number of individuals sampled from population group 1
#' 	    \item number_hosts_sampled_group_2, Number of individuals sampled from population group 2
#' 	    \item number_hosts_population_group_1, Estimated number of individuals in population group 1
#' 	    \item number_hosts_population_group_2, Estimated number of individuals in population group 2
#' 	    \item max_possible_pairs_in_sample, Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2
#' 	    \item max_possible_pairs_in_population, Number of distinct possible transmission pairs between individuals in population groups 1 and 2
#' 	    \item num_linked_pairs_observed, Number of observed directed transmission pairs between samples from population groups 1 and 2 
#' 	    \item p_hat, Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked                                      
#' 	    \item est_linkedpairs_in_population, Estimated transmission pairs between population groups 1 and 2
#' 	    \item theta_hat, Estimated transmission flows or relative probability of transmission within and between population groups 1 and 2 adjusted
#' 		      for sampling heterogeneity. More precisely, the conditional probability that a pair of pathogen sequences is from a specific population 
#' 		      group pairing given that the pair is linked.}
#' 
#' @seealso See \code{\link{estimate_p_hat}} to prepare input data to estimate \code{theta_hat}
#' 
#' @references
#' \enumerate{
#'
#' \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'       HIV transmission in the context of a universal testing and treatment
#'       trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#' \item Carnegie, N.B., et al., Linkage of viral sequences among HIV-infected
#' 	     village residents in Botswana: estimation of linkage rates in the 
#' 	     presence of missing data. PLoS Computational Biology, 2014. 10(1): 
#' 	     p. e1003430.
#'
#'}
#'
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#' 
#' # Estimate transmission flows within and between population groups accounting for variable
#' # sampling among population groups
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs and ?sampling_frequency
#' 
#' 
#' # Load and view data
#' #
#' # The input data comprises counts of observed directed HIV transmission pairs within 
#' # and between intervention and control communities in the BCPP/Ya Tsie trial,
#' # sampling information and the probability of linkage between individuals sampled
#' # from intervention and control communities (i.e. \code{p_hat})
#' #
#' # See ?estimate_p_hat() for details on estimating p_hat
#' results_estimate_p_hat <- estimated_hiv_transmission_flows[, c(1:10)]
#' 
#' results_estimate_p_hat
#' 
#' # Estimate theta_hat
#' results_estimate_theta_hat <- estimate_theta_hat(df_counts_and_p_hat = results_estimate_p_hat)
#' 
#' # View results
#' results_estimate_theta_hat
#' 
#' @export
estimate_theta_hat <- function(df_counts_and_p_hat, ...) UseMethod("estimate_theta_hat")

#' @describeIn estimate_theta_hat Estimates conditional probability of linkage (transmission flows)
#' 
#' @export
estimate_theta_hat.default <- function(df_counts_and_p_hat, ...) {

	# Check whether the input is a data.frame
	if (!is.data.frame(df_counts_and_p_hat)) {
		
	    stop("Check input! df_counts_and_p_hat should be a data.frame produced by function: estimate_p_hat(), type ?estimate_p_hat() for help")

	}

    results_get_theta_hat <- get_theta_hat(df_counts_and_p_hat) 

    results_get_theta_hat

}


#' @title \code{estimate_c_hat} Estimates probability of clustering 
#'   
#' @description This function estimates \code{c_hat}, the probability that a randomly
#'              selected pathogen sequence in one population group links to at least
#'              one pathogen sequence in another population group.
#'                 
#' @aliases c_hat est_c_hat estim_c_hat estimate_c_hat
#' 
#' @param df_counts_and_p_hat  A data.frame returned by the function: [estimate_p_hat()]
#' @param \dots Further arguments.
#' 
#' @return Returns a data.frame containing:
#' \itemize{
#'     \item H1_group, Name of population group 1 
#'     \item H2_group, Name of population group 2 
#'     \item number_hosts_sampled_group_1, Number of individuals sampled from population group 1
#'     \item number_hosts_sampled_group_2, Number of individuals sampled from population group 2
#'     \item number_hosts_population_group_1, Estimated number of individuals in population group 1
#'     \item number_hosts_population_group_2, Estimated number of individuals in population group 2
#'     \item max_possible_pairs_in_sample, Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2
#'     \item max_possible_pairs_in_population, Number of distinct possible transmission pairs between individuals in population groups 1 and 2
#'     \item num_linked_pairs_observed, Number of observed directed transmission pairs between samples from population groups 1 and 2 
#'     \item p_hat, Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked                                      
#'     \item c_hat, Probability that a randomly selected pathogen sequence in one population group links to at least 
#'           one pathogen sequence in another population group i.e. probability of clustering
#' }
#' 
#' @seealso See \code{\link{estimate_p_hat}} to prepare input data to estimate \code{c_hat}
#'  
#' @references
#' \enumerate{
#'
#' \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'       HIV transmission in the context of a universal testing and treatment
#'       trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#' \item Carnegie, N.B., et al., Linkage of viral sequences among HIV-infected
#' 	     village residents in Botswana: estimation of linkage rates in the 
#' 	     presence of missing data. PLoS Computational Biology, 2014. 10(1): 
#' 	     p. e1003430.
#'
#'}
#' 
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#' 
#' # Estimate the probability of clustering between individuals from two population groups of interest
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs, ?sampling_frequency and ?estimated_hiv_transmission_flows
#' 
#' 
#' # Load and view data
#' #
#' # The input data comprises counts of observed directed HIV transmission pairs within and
#' # between intervention and control communities in the BCPP/Ya Tsie trial, sampling 
#' # information and the probability of linkage between individuals sampled from
#' # intervention and control communities (i.e. \code{p_hat})
#' #
#' # See ?estimate_p_hat() for details on estimating p_hat
#' results_estimate_p_hat <- estimated_hiv_transmission_flows[, c(1:10)]
#' 
#' results_estimate_p_hat
#' 
#' # Estimate c_hat
#' results_estimate_c_hat <- estimate_c_hat(df_counts_and_p_hat = results_estimate_p_hat)
#' 
#' # View results
#' results_estimate_c_hat
#' 
#' @export
estimate_c_hat <- function(df_counts_and_p_hat, ...) UseMethod("estimate_c_hat")

#' @describeIn estimate_c_hat Estimates probability of clustering
#' 
#' @export
estimate_c_hat.default <- function(df_counts_and_p_hat, ...) {

	# Check whether the input is a data.frame
	if (!is.data.frame(df_counts_and_p_hat)) {
		
	    stop("Check input! df_counts_and_p_hat should be a data.frame produced by function: estimate_p_hat(), type ?estimate_p_hat() for help")

	}

    results_get_c_hat <- get_c_hat(df_counts_and_p_hat) 

    results_get_c_hat

}


#' @title \code{estimate_multinom_ci} Estimates confidence intervals for transmission flows
#' 
#' @description This function computes simultaneous confidence intervals at the
#'              5% significance level for estimated transmission flows. Available
#'              methods for computing confidence intervals are: Goodman, Goodman
#'              with a continuity correction, Sison-Glaz and Queensbury-Hurst.
#'              
#' @aliases estimate_multinom_ci estim_multinom_ci est_multinom_ci
#' 
#' @param df_theta_hat A data.frame returned by the function: [estimate_theta_hat()]
#' @param detailed_report A boolean value to produce detailed output of the analysis. (Default is \code{FALSE})                       
#' @param \dots Further arguments.
#' 
#' @return Returns a data.frame containing:
#' \itemize{
#' 	\item H1_group, Name of population group 1 
#' 	\item H2_group, Name of population group 2 
#' 	\item number_hosts_sampled_group_1, Number of individuals sampled from population group 1
#' 	\item number_hosts_sampled_group_2, Number of individuals sampled from population group 2
#' 	\item number_hosts_population_group_1, Estimated number of individuals in population group 1
#' 	\item number_hosts_population_group_2, Estimated number of individuals in population group 2
#' 	\item max_possible_pairs_in_sample, Number of distinct possible transmission pairs between individuals sampled from population groups 1 and 2
#' 	\item max_possible_pairs_in_population, Number of distinct possible transmission pairs between individuals in population groups 1 and 2
#' 	\item num_linked_pairs_observed, Number of observed directed transmission pairs between samples from population groups 1 and 2 
#' 	\item p_hat, Probability that pathogen sequences from two individuals randomly sampled from their respective population groups are linked                                      
#' 	\item est_linkedpairs_in_population, Estimated transmission pairs between population groups 1 and 2
#' 	\item theta_hat, Estimated transmission flows or relative probability of transmission within and between population groups 1 and 2 adjusted
#' 	      for sampling heterogeneity. More precisely, the conditional probability that a pair of pathogen sequences is from a specific population 
#' 	      group pairing given that the pair is linked.
#' 	\item obs_trm_pairs_est_goodman, Point estimate, Goodman method Confidence intervals for observed transmission pairs
#' 	\item obs_trm_pairs_lwr_ci_goodman, Lower bound of Goodman confidence interval 
#' 	\item obs_trm_pairs_upr_ci_goodman, Upper bound of Goodman confidence interval 
#' 	\item est_goodman, Point estimate, Goodman method Confidence intervals for estimated transmission flows
#' 	\item lwr_ci_goodman, Lower bound of Goodman confidence interval 
#' 	\item upr_ci_goodman, Upper bound of Goodman confidence interval 
#' }
#' 
#' The following additional fields are returned if the detailed_report flag is set
#' 
#' \itemize{
#' 	\item est_goodman_cc, Point estimate, Goodman method Confidence intervals with continuity correction
#' 	\item lwr_ci_goodman_cc, Lower bound of Goodman confidence interval 
#' 	\item upr_ci_goodman_cc, Upper bound of Goodman confidence interval 
#' 	\item est_sisonglaz, Point estimate, Sison-Glaz method Confidence intervals
#' 	\item lwr_ci_sisonglaz, Lower bound of Sison-Glaz confidence interval 
#' 	\item upr_ci_sisonglaz, Upper bound of Sison-Glaz confidence interval 
#' 	\item est_qhurst_acswr, Point estimate, Queensbury-Hurst method Confidence intervals via ACSWR r package 
#' 	\item lwr_ci_qhurst_acswr, Lower bound of Queensbury-Hurst confidence interval 
#' 	\item upr_ci_qhurst_acswr, Upper bound of Queensbury-Hurst confidence interval 
#' 	\item est_qhurst_coinmind, Point estimate, Queensbury-Hurst method Confidence intervals via CoinMinD r package 
#' 	\item lwr_ci_qhurst_coinmind, Lower bound of Queensbury-Hurst confidence interval
#' 	\item upr_ci_qhurst_coinmind, Upper bound of Queensbury-Hurst confidence interval
#' 	\item lwr_ci_qhurst_adj_coinmind, Lower bound of Queensbury-Hurst confidence interval adjusted 
#' 	\item upr_ci_qhurst_adj_coinmind, Upper bound of Queensbury-Hurst confidence interval adjusted
#' }
#' 
#' @seealso See \code{\link{estimate_theta_hat}} to prepare input data to estimate confidence intervals.
#' 
#'          To learn more about the Goodman and Sison-Glaz confidence interval methods 
#'          see \code{\link[DescTools]{MultinomCI}}. For Queensbury-Hurst confidence
#'          intervals see \code{\link[ACSWR]{QH_CI}} and \code{\link[CoinMinD]{QH}}
#' 
#' @references
#' \enumerate{
#'
#'    \item Magosi LE, et al., Deep-sequence phylogenetics to quantify patterns of 
#'       HIV transmission in the context of a universal testing and treatment
#'       trial – BCPP/ Ya Tsie trial. To submit for publication, 2021.
#' 
#'    \item Goodman, L. A. On Simultaneous Confidence Intervals for Multinomial Proportions 
#' 		 Technometrics, 1965. 7, 247-254.
#' 
#'    \item Cherry, S., A Comparison of Confidence Interval Methods for Habitat 
#'       Use-Availability Studies. The Journal of Wildlife Management, 1996. 
#'       60(3): p. 653-658.
#' 
#'    \item Sison, C.P and Glaz, J. Simultaneous confidence intervals and sample size determination
#' 		 for multinomial proportions. Journal of the American Statistical Association, 
#' 		 1995. 90:366-369.
#' 
#'    \item Glaz, J., Sison, C.P. Simultaneous confidence intervals for multinomial proportions. 
#' 		 Journal of Statistical Planning and Inference, 1999. 82:251-262.
#' 
#'    \item May, W.L., Johnson, W.D. Constructing two-sided simultaneous confidence intervals for 
#' 		 multinomial proportions for small counts in a large number of cells. 
#' 		 Journal of Statistical Software, 2000. 5(6).
#' 		 Paper and code available at https://www.jstatsoft.org/v05/i06.
#' 
#'    \item Carnegie, N.B., et al., Linkage of viral sequences among HIV-infected
#' 		 village residents in Botswana: estimation of linkage rates in the 
#' 		 presence of missing data. PLoS Computational Biology, 2014. 10(1): 
#' 		 p. e1003430.
#' 
#'    \item Ratmann, O., et al., Inferring HIV-1 transmission networks and sources of 
#' 		 epidemic spread in Africa with deep-sequence phylogenetic analysis. 
#' 		 Nature Communications, 2019. 10(1): p. 1411.
#' 
#'    \item Wymant, C., et al., PHYLOSCANNER: Inferring Transmission from Within- and
#' 		 Between-Host Pathogen Genetic Diversity. Molecular Biology and Evolution,
#' 		 2017. 35(3): p. 719-733.
#' }
#' 
#' @examples 
#' library(bumblebee)
#' library(dplyr)
#' 
#' # Compute confidence intervals for estimated transmission flows
#' 
#' # We shall use the data of HIV transmissions within and between intervention and control
#' # communities in the BCPP/Ya Tsie HIV prevention trial. To learn more about the data 
#' # ?counts_hiv_transmission_pairs and ?sampling_frequency
#' 
#' 
#' # Load and view data
#' #
#' # The data comprises counts of observed directed HIV transmission pairs between individuals 
#' # sampled from intervention and control communities (i.e. num_linked_pairs_observed); 
#' # and the estimated HIV transmissions within and between intervention and control 
#' # communities in the BCPP/Ya Tsie trial population adjusted for sampling heterogneity
#' # (i.e. \code{est_linkedpairs_in_population}). See ?estimate_theta_hat() for details on 
#' # computing \code{est_linkedpairs_in_population} and \code{theta_hat}.
#' 
#' results_estimate_theta_hat <- estimated_hiv_transmission_flows[, c(1:13)]
#'  
#' results_estimate_theta_hat
#'  
#' # Compute Goodman confidence intervals (Default)
#' results_estimate_multinom_ci <- estimate_multinom_ci(
#'     df_theta_hat = results_estimate_theta_hat, 
#'     detailed_report = FALSE)
#' 
#' # View results
#' results_estimate_multinom_ci
#' 
#' # Compute Goodman, Sison-Glaz and Queensbury-Hurst confidence intervals
#' results_estimate_multinom_ci_detailed <- estimate_multinom_ci(
#'     df_theta_hat = results_estimate_theta_hat, 
#'     detailed_report = TRUE)
#' 
#' # View results
#' results_estimate_multinom_ci_detailed
#' 
#' @export
estimate_multinom_ci <- function(df_theta_hat, ...) UseMethod("estimate_multinom_ci")

#' @describeIn estimate_multinom_ci Estimates confidence intervals for transmission flows

#' @export
estimate_multinom_ci.default <- function(df_theta_hat, detailed_report = FALSE, ...) {

	# Check whether the input is a data.frame
	if (!is.data.frame(df_theta_hat)) {
		
	    stop("Check input! df_theta_hat should be a data.frame produced by function: estimate_theta_hat(), type ?estimate_theta_hat() for help")

	}

    # check input fields
	provided_field_names <- names(df_theta_hat)
	expected_field_names <- c("est_linkedpairs_in_population", "num_linked_pairs_observed")   

	if (base::isFALSE(base::all(expected_field_names %in% provided_field_names))) {

	    message <- list('Check input! One or more required fields is missing. \nExpecting a data.frame produced by function: estimate_theta_hat(), type ?estimate_theta_hat() for help. \nTip! Minimum required fields: "est_linkedpairs_in_population" and "num_linked_pairs_observed"')
	    stop(message)

	}

    results_get_multinomial_proportion_conf_ints_extended <- get_multinomial_proportion_conf_ints_extended(df_theta_hat, detailed_report) 

    results_get_multinomial_proportion_conf_ints_extended
    
}
