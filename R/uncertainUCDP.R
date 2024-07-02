

#' Parameter extraction for uncertainUCDP-functions
#'
#' @description
#' Extracting parameters for the reported-value inflated Gumbel mixture distribution for UCDP events. Primarily intended for internal use by the uncertainUCDP-functions, but can be used to extract parameters for the distribution manually.
#'
#' @param fatalities A vector of non-negative integers representing the number of fatalities of the UCDP event. Non-integer values are allowed but should be considered experimental
#' @param tov A character string or integer value representing the type of violence of the UCDP. Must be one of "sb", "ns", "os", or "any" or their numeric equivalent The options are:
#'
#' * "sb" or 1 for state-based violence
#' * "ns" or 2 for non-state violence
#' * "os" or 3 for one-sided violence
#' * "any" or 4 for parameters estimated across all type of violence. This is somewhat experimental and should be used with caution. This is possibly useful when the type of violence is unknown or when the user wants to combine all types of violence into a single category.
#'
#' @return A list with three elements: loc, scale, and w. loc and scale are the location and scale parameters of the Gumbel distribution, respectively. w is the weight parameter for the reported-value inflation
#'
#' @export
#'
#'
uncertainUCDP_parameters <- function(fatalities, tov){

	if(is.numeric(tov)){
		tov <- switch(tov,
					  '1' = 'sb',
					  '2' = 'ns',
					  '3' = 'os',
					  '4' = 'any')
	}
	tov <- match.arg(tov, c('sb','ns','os','any'))

	if(!(tov %in% c('sb','ns','os','any'))){
		stop('tov must be one of "sb", "ns", "os", or "any"')
	}

	if(tov == 'any'){
		warning("tov = 'any' is somewhat experimental as it combines all types of violence into a single category. Use with caution")
	}

	if(any(fatalities < 0)){
		warning('fatalities must be non-negative integers, replacing negative values with NA')
		fatalities[fatalities < 0] <- NA
	}

	if(any(fatalities != round(fatalities))){
		warning('fatalities should be an interger. Using non-integers may lead to unexpected results and should be considered experimental')
	}

	data <- tibble::tibble(fatlev = fatalities,
								 logfatlev = base::log1p(.data$fatlev),
								 fl_1 = dplyr::case_when(.data$fatlev == 1 ~ 1,
								 								 T~0),
								 fl_2 = dplyr::case_when(.data$fatlev == 2 ~ 1,
								 								 T~0),
								 fl_3 = dplyr::case_when(.data$fatlev == 3 ~ 1,
								 								 T~0),
								 fl_13 = dplyr::case_when(.data$fatlev == 13 ~ 1,
								 									T~0),
								 fl_20 = dplyr::case_when(.data$fatlev == 20 ~ 1,
								 									T~0),
								 fl_24 = dplyr::case_when(.data$fatlev == 24 ~ 1,
								 									T~0),
								 fl_40 = dplyr::case_when(.data$fatlev == 40 ~ 1,
								 									T~0),
								 fl_101 = dplyr::case_when(.data$fatlev == 101 ~ 1,
								 									 T~0),
								 fl_200 = dplyr::case_when(.data$fatlev == 200 ~ 1,
								 									 T~0),
								 fl_1001 = dplyr::case_when(.data$fatlev == 1001 ~ 1,
								 										T~0),
								 fl_2000 = dplyr::case_when(.data$fatlev == 2000 ~ 1,
								 										T~0))


	if(tov == 'sb'){
	loc <- expm1(stats::predict(internal_models$sb[[1]], newdata = data))
	scale <- exp(stats::predict(internal_models$sb[[2]], newdata = data))
	w <- stats::predict(internal_models$sb[[3]], newdata = data, type = 'r')
	}else if(tov == 'ns'){
	loc <- expm1(stats::predict(internal_models$ns[[1]], newdata = data))
	scale <- exp(stats::predict(internal_models$ns[[2]], newdata = data))
	w <- stats::predict(internal_models$ns[[3]], newdata = data, type = 'r')
	}else if(tov == 'os'){
	loc <- expm1(stats::predict(internal_models$os[[1]], newdata = data))
	scale <- exp(stats::predict(internal_models$os[[2]], newdata = data))
	w <- stats::predict(internal_models$os[[3]], newdata = data, type = 'r')
	}else if(tov == 'any'){
	loc <- expm1(stats::predict(internal_models$any[[1]], newdata = data))
	scale <- exp(stats::predict(internal_models$any[[2]], newdata = data))
	w <- stats::predict(internal_models$any[[3]], newdata = data, type = 'r')
	}

	return(list(loc = loc, scale = scale, w = w))


}




#' Parametric uncertainty distributions for UCDP events
#'
#' @description
#' Density, distribution, quantile and random number generation functions for the parametric reported-value inflated Gumbel mixture distribution for UCDP events. The functions estimate the parameters of the distribution based on the number of fatalities and the type of violence of the UCDP event.
#'
#' @details
#' The reported-value inflated Gumbel mixture distribution is a parametric distribution for modeling the uncertainty in the number of fatalities of UCDP events. The distribution is a mixture of a Gumbel distribution and a point mass at the reported number of fatalities. The distribution is estimated based on the number of fatalities and the type of violence of the UCDP event. The distribution is estimated using a set of regression models that estimate the location, scale, and weight parameters of the distribution based on the number of fatalities and the type of violence of the UCDP event.
#'
#' @param n Number of observations to generate random values for
#' @param x,q Vector of quantiles
#' @param p Vector of probabilities
#' @param fatalities A vector of non-negative integers representing the number of fatalities of the UCDP events. Non-integer values are allowed but should be considered experimental.
#' @param tov A character string representing the type of violence of the UCDP. Must be one of "sb", "ns", "os", or "any". The options are:
#'
#' * "sb" for state-based violence
#' * "ns" for non-state violence
#' * "os" for one-sided violence
#' * "any" for parameters estimated across all type of violence. This is somewhat experimental and should be used with caution. This is possibly useful when the type of violence is unknown or when the user wants to combine all types of violence into a single category.
#'
#' @return
#'
#' * \code{duncertainUCDP} gives the density function
#' * \code{puncertainUCDP} gives the distribution function
#' * \code{quncertainUCDP} gives the quantile function
#' * \code{runcertainUCDP} generates random values as a vector of length \code{n}
#' @export
#'
#' @examples
#'
#' data(ucdpged)
#'
#' # Generate 10 random values for an arbitrary UCDP event
#' runcertainUCDP(n = 10, fatalities = 100, tov = 'sb')
#'
#' # Generate 10 random values for the first event in the GED sample
#' runcertainUCDP(n = 10, fatalities = ucdpged$best[1], tov = ucdpged$type_of_violence[1])
#'
#' # Obtaining the probability that an arbitrary UCDP event has at least 150 fatalities
#' puncertainUCDP(q = 150, fatalities = 100, tov = 'ns')
#'
#' # Obtaining the probability that the for the first event in the GED sample has at least 5 fatalities
#' puncertainUCDP(q = 5, fatalities = ucdpged$best[1], tov = ucdpged$type_of_violence[1])
#'
#' # Obtaining the 90th percentile for an arbitrary UCDP event and one-sided violence
#' quncertainUCDP(p = 0.9, fatalities = 100, tov = 'os')
#'
#' # Obtaining the 90th percentile for the first event in the GED sample
#' quncertainUCDP(p = 0.9, fatalities = ucdpged$best[1], tov = ucdpged$type_of_violence[1])
#'
#' # Obtaining the density for an arbitrary UCDP event and state-based violence
#' duncertainUCDP(x = seq(from = 0, to = 500), fatalities = 100, tov = 'sb')
#'
#' # Obtaining the density for the first event in the GED sample
#' duncertainUCDP(x = seq(0, 50), fatalities = ucdpged$best[1], tov = ucdpged$type_of_violence[1])
#'
runcertainUCDP <- function(n, fatalities, tov = c('sb','ns','os','any')){

	if(length(n) != 1){
		stop('n must be a single integer')
	}

	params <- uncertainUCDP_parameters(fatalities, tov)

	gumbels <- rgumbel(n, params$loc, params$scale)
	rv_inflation <- stats::rbinom(n, 1, params$w)

	return(gumbels * (1-rv_inflation) + fatalities * rv_inflation)
}

#' @rdname runcertainUCDP
#' @export
puncertainUCDP <- function(q, fatalities, tov = c('sb','ns','os','any')){

	params <- uncertainUCDP_parameters(fatalities, tov)

	return(mistr::pgumbel(q, params$loc, params$scale) * (1-params$w) + q>=fatalities * params$w)

}

#' @rdname runcertainUCDP
#' @export
duncertainUCDP <- function(x, fatalities, tov = c('sb','ns','os','any')){

	params <- uncertainUCDP_parameters(fatalities, tov)

	return(mistr::dgumbel(x, params$loc, params$scale) * (1-params$w) + x==fatalities * params$w)
}

#' @rdname runcertainUCDP
#' @export
quncertainUCDP <- function(p, fatalities, tov = c('sb','ns','os','any')){

	params <- uncertainUCDP_parameters(fatalities, tov)

	infliction_point <- mistr::pgumbel(fatalities, params$loc, params$scale) * (1-params$w)

	return(mistr::qgumbel(p, params$loc, params$scale) * (1-params$w) + (p>=infliction_point) * params$w)

}


#' Mean, median, and quantiles of the parametric uncertainty distributions for UCDP events
#'
#' @description
#' Mean, median, and quantiles of the parametric uncertainty distributions for UCDP events. The parametric uncertainty distributions are based on the reported-value inflation Gumbel mixture distribution. The \code{median} and \code{quantile} functions are shortcuts for the \code{quncertainUCDP} function.
#'
#'
#' @param fatalities A vector of non-negative integers representing the number of fatalities of the UCDP events. Non-integer values are allowed but should be considered experimental.
#' @param tov A character string representing the type of violence of the UCDP. Must be one of "sb", "ns", "os", or "any". The options are:
#'
#' * "sb" for state-based violence
#' * "ns" for non-state violence
#' * "os" for one-sided violence
#' * "any" for parameters estimated across all type of violence. This is somewhat experimental and should be used with caution. This is possibly useful when the type of violence is unknown or when the user wants to combine all types of violence into a single category.
#' @param probs A numeric vector of probabilities with values in [0,1]. The quantiles to calculate.
#'
#'
#' @return A numeric vector of the same length as the input vector of fatalities representing the means, medians, and quantiles of the parametric uncertainty distribution for each UCDP event.
#' @export
#'
#' @examples
#'
#' data(ucdpged)
#'
#' # Calculate the mean for an arbitrary UCDP event
#' mean_uncertainUCDP(fatalities = 100, tov = 'sb')
#'
#' # Calculate the mean for the first event in the UCDP GED sample
#' mean_uncertainUCDP(ucdpged$best[1], tov = ucdpged$type_of_violence[1])
#'
#' # Calculate the median for an arbitrary UCDP event
#' median_uncertainUCDP(fatalities = 100, tov = 'sb')
#'
#' # Calculate the median for the first event in the UCDP GED sample
#' median_uncertainUCDP(ucdpged$best[1], tov = ucdpged$type_of_violence[1])
#'
#' # Calculate the 90th percentile for an arbitrary UCDP event
#' quantiles_unceartainUCDP(probs = 0.9, fatalities = 100, tov = 'sb')
#'
#' # Calculate the 90th percentile for the first event in the UCDP GED sample
#' quantiles_unceartainUCDP(ucdpged$best[1], 0.9, tov = ucdpged$type_of_violence[1])
#'
mean_uncertainUCDP <- function(fatalities, tov = c('sb','ns','os','any')){

	params <- uncertainUCDP_parameters(fatalities, tov)

	return((params$loc + params$scale * -digamma(1))*(1-params$w) + fatalities * params$w)
}

#' @rdname mean_uncertainUCDP
#' @export
median_uncertainUCDP <- function(fatalities, tov = c('sb','ns','os','any')){

	params <- uncertainUCDP_parameters(fatalities, tov)

	puncertainUCDP(0.5, fatalities, tov)
}

#' @rdname mean_uncertainUCDP
#' @export
quantiles_unceartainUCDP <- function(probs, fatalities, tov = c('sb','ns','os','any')){

	params <- uncertainUCDP_parameters(fatalities, tov)

	return(quncertainUCDP(probs, fatalities, tov))
}


rgumbel <- function (n, loc, scale){
	if (any(scale <= 0)) {
		warning("NaNs produced")
		return(rep.int(NaN, n))
	}
	loc - scale * log(-log(stats::runif(n)))
}









