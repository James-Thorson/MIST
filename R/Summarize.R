#' Summarize model outputs
#'
#' \code{Summarize} calculates useful outputs
#'
#' @export
Summarize = function( Report, SD, InputList, SimList=NULL, species_names=NULL ){

  # Local function to extract standard errors in useful way
  SE_hat_fn = function( SD, Report, Map=NULL, parname){
    Return = array(NA, dim=dim(Report[[parname]]))
    if( !is.null(Map) ){
      if( parname %in% names(Map) ) Return[which(!is.na(Map[[parname]]))] = summary(SD)[which(parname==rownames(summary(SD))),'Std. Error']
      if( !(parname %in% names(Map)) ) Return[] = summary(SD)[which(parname==rownames(summary(SD))),'Std. Error']
    }
    if( is.null(Map) ) Return[] = summary(SD)[which(parname==rownames(summary(SD))),'Std. Error']
    return(Return)
  }
  DerivedQuants = NULL

  # Stationary covariance
  DerivedQuants[["Vinf_pp_hat"]] = Report$Vinf_pp  # Equals matrix( solve(diag(Nspecies^2) - kronecker(Report$B_pp,Report$B_pp)) %*% as.vector(Report$Cov_pp), nrow=Nspecies, ncol=Nspecies)
  #if( !is.null("SimList")) DerivedQuants[["Vinf_pp"]] = SimList$Vinf_pp

  # Interaction matrix
  DerivedQuants[["B_pp_hat"]] = Report$B_pp
  colnames(DerivedQuants[["B_pp_hat"]]) = rownames(DerivedQuants[["B_pp_hat"]]) = levels(DF$spp)
  if( !is.null(SD) ){
    DerivedQuants[["B_pp_hat_SE"]] = SE_hat_fn( SD=SD, Report=Report, parname="B_pp")
  }
  if( !is.null("SimList")) DerivedQuants[["B_pp"]] = SimList$B_pp

  # Process error variance
  DerivedQuants[["Cov_pp_hat"]] = Report$Cov_pp
  colnames(DerivedQuants[["Cov_pp_hat"]]) = rownames(DerivedQuants[["Cov_pp_hat"]]) = levels(DF$spp)
  if( !is.null(SimList)) DerivedQuants[["Cov_pp"]] = SimList$Cov_pp

  # Resilience
  DerivedQuants[["MeanPropVar_hat"]] = Report$MeanPropVar
  if( !is.null("SimList")) DerivedQuants[["MeanPropVar"]] = SimList$MeanPropVar

  # Harmonic mean of eigenvalue (i.e., "average rho")
  DerivedQuants[["Harmonic_Mean_Rho"]] = sqrt(Report$MeanPropVar)
  if( !is.null(SimList)) DerivedQuants[["MeanPropVar"]] = sqrt(SimList$MeanPropVar)

  # Reactivity
  DerivedQuants[["Reactivity_hat"]] = Report$Reactivity
  if( !is.null(SimList)) DerivedQuants[["Reactivity"]] = SimList$Reactivity

  # Density dependence
  DerivedQuants[["Eigen_B_hat"]] = eigen( Report$B_pp )
  if( !is.null(SimList)) DerivedQuants[["Eigen_B"]] = eigen( SimList$B_pp )

  # Density dependence
  DerivedQuants[["Eigen_BminusI_hat"]] = eigen( Report$B_pp - diag(InputList$TmbData$n_p))
  if( !is.null(SimList)) DerivedQuants[["Eigen_BminusI"]] = eigen( SimList$B_pp - diag(InputList$TmbData$n_p))

  # Terms
  DerivedQuants[["Alpha_pr"]] = Report$Alpha_pr
  DerivedQuants[["Beta_pr"]] = Report$Beta_pr
  if( !is.null(SD) ){
    DerivedQuants[["Alpha_pr_SE"]] = SE_hat_fn( SD=SD, Report=Report, Map=InputList$Map, parname="Alpha_pr")
    DerivedQuants[["Beta_pr_SE"]] = SE_hat_fn( SD=SD, Report=Report, Map=InputList$Map, parname="Beta_pr")
  }

  # Add names
  if( !is.null(DerivedQuants) ){
    for(i in 1:length(DerivedQuants)){
      if( is.matrix(DerivedQuants[[i]]) && nrow(DerivedQuants[[i]])==length(species_names) && ncol(DerivedQuants[[i]])==length(species_names) ){
        rownames(DerivedQuants[[i]]) = species_names
        colnames(DerivedQuants[[i]]) = species_names
      }
    }
  }

  # Return
  return( DerivedQuants )
}
