#' @import xtable
#' @import R2WinBUGS
#' @import forestplot




#' Fixed effect NMA BUGS script for binary outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k])			# binomial likelihood
      logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]]	# model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k]			# expected value of the numerators 
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k]))     #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))   
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  d[1]<-0						# treatment effect is zero for reference treatment
  for (k in 2:nt){d[k]~dnorm(0,.001)%_%I(-4,4)}		# vague priors for treatment effects
  
  
  #Output
  
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
  
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
}











#' Fixed effect inconsistency NMA BUGS script for binary outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k])			# binomial likelihood
      logit(p[i,k]) <- mu[i] +  d[t[i,1],t[i,k]]	# model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k]			# expected value of the numerators 
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k]))     #Deviance contribution
                     +  (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))   
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for(k in 1:nt) { d[k,k] <- 0 } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)) { # priors for all mean treatment effects 
    for(k in (c+1):nt){ d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }  	
  
}






#' Fixed effect NMA BUGS script for binary outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Binary_Covariate <- function(){
  # Binomial likelihood, logit link, subgroup
  # Fixed effects model with one covariate
  for(i in 1:ns){# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1
      logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4)
  
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}






#' Fixed effect inconsistency NMA BUGS script for binary outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Binary_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k])			# binomial likelihood
      logit(p[i,k]) <- mu[i] +  d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*x[i]	# model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k]			# expected value of the numerators 
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k]))     #Deviance contribution
                     +  (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))   
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)){ # priors for all mean treatment effects 
    for (k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  } 
  B~dnorm(0,.001)%_%I(-4,4)
  
}






#' Fixed effect NMA BUGS script for binary outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Continuous_Covariate <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model with continuous covariate
  for(i in 1:ns){# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1
      logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
      
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}





#' Fixed effect inconsistency NMA BUGS script for binary outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Continuous_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k])			# binomial likelihood
      logit(p[i,k]) <- mu[i] +  d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)	# model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k]			# expected value of the numerators 
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k]))     #Deviance contribution
                     +  (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))   
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects 
    for (k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  } 
  B~dnorm(0,.001)%_%I(-4,4)
  
}






#' Fixed effect NMA BUGS script for binary outcomes with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model with continuous covariate
  for(i in 1:ns){# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1
      logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
      
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}





#' Fixed effect inconsistency NMA BUGS script for binary outcomes with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Inconsist_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-3.5,3.5)			# vague priors for all trial baselines
    for(k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k])			# binomial likelihood
      logit(p[i,k]) <- mu[i] +  d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)	# model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k]			# expected value of the numerators 
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k]))     #Deviance contribution
                     +  (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))   
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-3.5,3.5)}
  } 
  B~dnorm(0,.001)%_%I(-3.5,3.5)
  
}





#' Random effect NMA BUGS script for binary outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k]~dnorm(0,.001)%_%I(-4,4)} # vague priors for treatment effects
  #sd~dunif(0,2) # vague prior for between-trial SD. ALTERNATIVES BELOW
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
  
  #Output
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]}
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
  
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
}






#' Random effect NMA BUGS script for binary outcomes (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k]~dnorm(0,.001)%_%I(-4,4)} # vague priors for treatment effects
  sd~dunif(0,5) # vague prior for between-trial SD. ALTERNATIVES BELOW
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  #Output
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]}
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
  
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
}







#' Random effect inconsistency NMA BUGS script for binary outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1] <- 0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
}







#' Random effect inconsistency NMA BUGS script for binary outcomes (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Inconsist_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1] <- 0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial precision
}







#' Random effect NMA BUGS script for binary outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Binary_Covariate <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k], taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  #sd~dunif(0,2) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}





#' Random effect NMA BUGS script for binary outcomes with a binary covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Binary_Covariate_Noninfo <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k], taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  sd~dunif(0,2) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]}
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}






#' Random effect inconsistency NMA BUGS script for binary outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Binary_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i] # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]],tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
}





#' Random effect inconsistency NMA BUGS script for binary outcomes with a binary covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Binary_Covariate_Inconsist_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i] # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]],tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial precision
}





#' Random effect NMA BUGS script for binary outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Continuous_Covariate <- function(){
  # Binomial likelihood, logit link, continuous covariate
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1 (centring)
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  #sd~dunif(0,2) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]}
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}





#' Random effect NMA BUGS script for binary outcomes with a continuous covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Continuous_Covariate_Noninfo <- function(){
  # Binomial likelihood, logit link, continuous covariate
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1 (centring)
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  sd~dunif(0,2) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}







#' Random effect inconsistency NMA BUGS script for binary outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Continuous_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
}





#' Random effect inconsistency NMA BUGS script for binary outcomes with a continuous covariate
#' (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Continuous_Covariate_Inconsist_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial precision
}







#' Random effect NMA BUGS script for binary outcomes with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Base_Risk <- function(){
  # Binomial likelihood, logit link, continuous covariate
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1 (centring)
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k], taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  #sd~dunif(0,2) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]}
  A <- sum(mu1[])/sum(nt1[])
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}





#' Random effect NMA BUGS script for binary outcomes with adjustment to baseline risk (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Base_Risk_Noninfo <- function(){
  # Binomial likelihood, logit link, continuous covariate
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k]~dbin(p[i,k],n[i,k]) # binomial likelihood
      # model for linear predictor, covariate effect relative to treat in arm 1 (centring)
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      rhat[i,k] <- p[i,k]*n[i,k] # expected value of the numerators
      dev[i,k] <- 2*(r[i,k]*(log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                     + (n[i,k]-r[i,k])*(log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k], taud[i,k])%_%I(-3.5,3.5) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-3.5,3.5) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-3.5,3.5) # vague prior for covariate effect
  sd~dunif(0,2) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      OR[c,k] <- exp(d[k] - d[c])
      RR[c,k] <-T[k]/T[c] 
    }  
  }
  
  
  # Provide estimates of treatment effects T[k] on the natural (probability) scale 
  # Absolute log odds on placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]}
  A <- sum(mu1[])/sum(nt1[])
  
  for (k in 1:nt){logit(T[k]) <- A + d[k]}
}







#' Random effect inconsistency NMA BUGS script for binary outcomes with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Inconsist_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx) # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for(k in 1:nt){ 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  V ~ dnorm(-2.34, 1.62)
  tau <- exp(-V)
}





#' Random effect inconsistency NMA BUGS script for binary outcomes with adjustment to baseline risk
#' (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_binomial_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Binomial_Inconsist_Base_Risk_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx) # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution
                       + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for(k in 1:nt){ 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial precision
}







#' Fixed effect NMA BUGS script for continuous outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] + d[t[i,k]] - d[t[i,1]] 	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1] <- 0 						# treatment effect is zero for reference treatment
  for (k in 2:nt){d[k]~dnorm(0,.001)} 		# vague priors for treatment effects
  
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){ 
    for(k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for (i in 1:ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  #A~dnorm(meanA, precA)
  A <- meanA
  
  for (k in 1:nt){T[k] <- A + d[k]}
}






#' Fixed effect NMA BUGS script for continuous outcomes in a SMD format
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_SMD <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    for (k in 2:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k], 2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k] ~ dnorm(theta[i,k], prec[i,k]) 		# normal likelihood
      theta[i,k] <- d[t[i,k]] - d[t[i,1]] 	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i, 2:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1]<-0 						# treatment effect is zero for reference treatment
  for (k in 2:nt){d[k] ~ dnorm(0, .001)} 		# vague priors for treatment effects
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt - 1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k] - d[c])
    }  
  }
  
}





#' Fixed effect inconsistency NMA BUGS script for continuous outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  d[t[i,1],t[i,k]]	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for(k in 1:nt) {d[k,k] <- 0} 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  }  	
  
}




#' Fixed effect inconsistency NMA BUGS script for continuous outcomes in a SMD format
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Inconsist_SMD <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    for (k in 2:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k], prec[i,k]) 		# normal likelihood
      theta[i,k] <- d[t[i,1],t[i,k]]	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { d[k,k] <- 0 } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects 
    for (k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)}
  }  	
  
}





#' Fixed effect NMA BUGS script for continuous outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Binary_Covariate <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*x[i] 	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1] <- 0 						# treatment effect is zero for reference treatment
  beta[1] <- 0        # covariate effect is zero for reference treatment  
  for (k in 2:nt){ 
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B
  } 		
  B~dnorm(0,.001)
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){ 
    for(k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  #A~dnorm(meanA, precA)
  A <- meanA
  
  for (k in 1:nt){T[k] <- A + d[k]}
}





#' Fixed effect NMA BUGS script for continuous outcomes in a SMD format with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Binary_Covariate_SMD <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    for (k in 2:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*x[i] 	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1]<-0 						# treatment effect is zero for reference treatment
  beta[1] <- 0        # covariate effect is zero for reference treatment  
  for (k in 2:nt){ 
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B
  } 		
  B~dnorm(0,.001)
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){ 
    for(k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
}





#' Fixed effect inconsistency NMA BUGS script for continuous outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Binary_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  d[t[i,1],t[i,k]]	 + (beta[t[i,k]]-beta[t[i,1]])*x[i]# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  } 
  B~dnorm(0,.001)
}






#' Fixed effect inconsistency NMA BUGS script for continuous outcomes in a SMD format
#'  with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Binary_Covariate_Inconsist_SMD <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    for (k in 2:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*x[i]# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  } 
  B~dnorm(0,.001)
  
}






#' Fixed effect NMA BUGS script for continuous outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Continuous_Covariate <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k] ~ dnorm(theta[i,k], prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]] - beta[t[i,1]])*(x[i] - mx) 	# model for linear predictor
      dev[i,k] <- (y[i,k] - theta[i,k])*(y[i,k] - theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1]<-0 						# treatment effect is zero for reference treatment
  beta[1] <- 0        # covariate effect is zero for reference treatment  
  for (k in 2:nt){ 
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B
  } 		
  B~dnorm(0,.001)
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  #A~dnorm(meanA, precA)
  A <- meanA
  
  for (k in 1:nt){T[k] <- A + d[k]}
}







#' Fixed effect NMA BUGS script for continuous outcomes in a SMD format with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Continuous_Covariate_SMD <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    for (k in 2:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) 	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1]<-0 						# treatment effect is zero for reference treatment
  beta[1] <- 0        # covariate effect is zero for reference treatment  
  for (k in 2:nt){ 
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B
  } 		
  B~dnorm(0,.001)
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
}






#' Fixed effect NMA BUGS script for continuous outcomes with adjustment for study baseline value
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Base_Risk <- function(){
  # Normal likelihood, identity link
  # Fixed effects model for multi-arm trials
  for(i in 1:ns){					#   LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx) 	# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) 		 #  summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[]) 			#Total Residual Deviance
  d[1]<-0 						# treatment effect is zero for reference treatment
  beta[1] <- 0        # covariate effect is zero for reference treatment  
  for(k in 2:nt){ 
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B
  } 		
  B~dnorm(0,.001)
  
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){ 
    for(k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for(i in 1: ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  A~dnorm(meanA, precA)
  
  for (k in 1:nt){T[k] <- A + d[k]}
}







#' Fixed effect inconsistency NMA BUGS script for continuous outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Continuous_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0, 0.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  d[t[i,1],t[i,k]]	 + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  } 
  B~dnorm(0,.001)
  
}






#' Fixed effect inconsistency NMA BUGS script for continuous outcomes in a SMD format with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Continuous_Covariate_Inconsist_SMD <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    for (k in 2:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- d[t[i,1],t[i,k]]	 + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  } 
  B~dnorm(0,.001)
  
}





#' Fixed effect inconsistency NMA BUGS script for continuous outcomes with adjustment to study
#' baseline value.
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Normal_Inconsist_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  d[t[i,1],t[i,k]]	 + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)# model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]     #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt) { 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)){ # priors for all mean treatment effects 
    for (k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  } 
  B~dnorm(0,.001)
  
}






#' Random effect NMA BUGS script for continuous outcomes
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- mu[i] + delta[i,k] # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k]~dnorm(0,.001)} # vague priors for treatment effects
  sd~dunif(0,12) # vague prior for between-trial SD.
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for (i in 1: ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  A~dnorm(meanA, precA)
  
  for (k in 1:nt){T[k] <- A + d[k]}
}






#' Random effect NMA BUGS script for continuous outcomes in a SMD format
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_SMD <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- delta[i,k] # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k]~dnorm(0,.001)} # vague priors for treatment effects
  V ~ dgamma(0.94, .00005)%_%I(0,325)
  U <- 1/V
  Var <- U*step(2-U) + 2*step(U-2)
  tau <- 1/Var # between-trial precision
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
}





RE_NMA_Normal_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  delta[i,k]
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)}
  }
  sd ~ dunif(0,5) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial precision
}





RE_NMA_Normal_Inconsist_SMD <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1] <- 0 # treatment effect is zero in control arm
    for(k in 2:na[i]){ # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- delta[i,k]
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial
    for(k in 2:na[i]){ # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  }
  V ~ dgamma(0.21, .0003)%_%I(0,350)
  U <- 1/V
  Var <- U*step(2-U) + 2*step(U-2)
  tau <- 1/Var # between-trial precision
}






#' Random effect NMA BUGS script for continuous outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Binary_Covariate <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i] # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  } 
  B~dnorm(0,.001) # vague prior for covariate effect
  sd~dunif(0,5) # vague prior for between-trial SD.
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for (i in 1:ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  A~dnorm(meanA, precA)
  
  for (k in 1:nt){T[k] <- A + d[k]}
}






#' Random effect NMA BUGS script for continuous outcomes in a SMD format with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Binary_Covariate_SMD <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i] # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  } 
  B~dnorm(0,.001) # vague prior for covariate effect
  V ~ dgamma(0.94, .00005)%_%I(0,350)
  U <- 1/V
  Var <- U*step(2-U) + 2*step(U-2)
  tau <- 1/Var # between-trial precision
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
}







#' Random effect NMA BUGS script for continuous outcomes with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Binary_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  B~dnorm(0,.001)
  sd ~ dunif(0,5) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial variance
}






#' Random effect NMA BUGS script for continuous outcomes in a SMD format with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Binary_Covariate_Inconsist_SMD <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for (k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  }
  B~dnorm(0,.001)
  #sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  V ~ dgamma(0.94, .00005)%_%I(0,350)
  U <- 1/V
  Var <- U*step(4-U) + 4*step(U-4)
  tau <- 1/Var # between-trial precision
}







#' Random effect NMA BUGS script for continuous outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Continuous_Covariate <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  } 
  B~dnorm(0,.001) # vague prior for covariate effect
  sd~dunif(0,5) # vague prior for between-trial SD.
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for(i in 1: ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  A~dnorm(meanA, precA)
  
  for (k in 1:nt){T[k] <- A + d[k]}
}






#' Random effect NMA BUGS script for continuous outcomes in a SMD format with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Continuous_Covariate_SMD <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  } 
  B~dnorm(0,.001) # vague prior for covariate effect
  V ~ dgamma(0.94, .00005)%_%I(0,350)
  U <- 1/V
  Var <- U*step(2-U) + 2*step(U-2)
  tau <- 1/Var # between-trial precision
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){ 
    for (k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
}






#' Random effect NMA BUGS script for continuous outcomes with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Continuous_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1] <- 0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  B~dnorm(0,.001)
  sd ~ dunif(0,5) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial precision
}





#' Random effect NMA BUGS script for continuous outcomes in a SMD format with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_SMD}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Continuous_Covariate_Inconsist_SMD <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for(k in 1:nt){ 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)}
  }
  B~dnorm(0,.001)
  #sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  V ~ dgamma(0.94, .00005)%_%I(0,350)
  U <- 1/V
  Var <- U*step(2-U) + 2*step(U-2)
  tau <- 1/Var # between-trial precision
}





#' Random effect NMA BUGS script for continuous outcomes with adjustment to study baseline value 
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Base_Risk <- function(){
  # Normal likelihood, identity link
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001) # vague priors for all trial baselines
    for(k in 1:na[i]){# LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) # calculate variances
      prec[i,k] <- 1/var[i,k] # set precisions
      y[i,k]~dnorm(theta[i,k],prec[i,k]) # normal likelihood
      theta[i,k] <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx) # model for linear predictor
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){
    d[k]~dnorm(0,.001) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  } 
  B~dnorm(0,.001) # vague prior for covariate effect
  sd~dunif(0,5) # vague prior for between-trial SD.
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  #Output
  # pairwise treatment effect for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){ 
    for(k in (c+1):nt){
      D[c,k] <- (d[k]-d[c])
    }  
  }
  
  
  # Provide estimates of effects T[k]  
  # Absolute effects with placebo treatment based on number of placebo controlled trials
  for(i in 1: ns){
    prec1[i] <- 1/pow(se[i,1], 2)*equals(t[i,1],1)
    mu1[i] <- mu[i]*prec1[i] 
  }
  
  precA <- sum(prec1[])
  meanA <- sum(mu1[])/precA
  A~dnorm(meanA, precA)
  
  for (k in 1:nt){T[k] <- A + d[k]}
}





#' Random effect NMA BUGS script for continuous outcomes with adjustment to study baseline value
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_normal_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Normal_Inconsist_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1] <- 0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2) 			# calculate variances
      prec[i,k] <- 1/var[i,k] 			 # set precisions
      y[i,k] ~ dnorm(theta[i,k],prec[i,k]) 		# normal likelihood
      theta[i,k] <- mu[i] +  delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- (y[i,k] - theta[i,k])*(y[i,k] - theta[i,k])*prec[i,k] #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)}
  }
  B~dnorm(0,.001)
  sd ~ dunif(0,5) # vague prior for between-trial standard deviation
  tau <- pow(sd, -2) # between-trial precision
}







#' Fixed effect NMA BUGS script for count data
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#' @return None
#' @examples
#'
#' @export 
FE_NMA_Poisson <- function(){ # *** PROGRAM STARTS
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k] ~ dnorm(0,.001)%_%I(-4,4)} # vague priors for treatment effects
  
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt){ 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for (k in 1:nt){log(T[k]) <- A + d[k]}
} 






#' Fixed effect inconsistency NMA BUGS script for count data
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Poisson_Inconsist <- function(){
  # Poisson likelihood, log link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      log(lambda[i,k]) <- mu[i] + d[t[i,1],t[i,k]] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution   
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for(k in 1:nt){d[k,k] <- 0} 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)){ # priors for all mean treatment effects 
    for (k in (c+1):nt){
      d[c,k] ~ dnorm(0,0.001)%_%I(-4,4)
    }
  }  	
  
}






#' Fixed effect NMA BUGS script for count data with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Poisson_Binary_Covariate <- function(){
  # Binomial likelihood, logit link, subgroup
  # Fixed effects model with one covariate
  for(i in 1:ns){# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,0.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,0.001)%_%I(-4,4)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)) {  
    for (k in (c+1):nt)  { 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for (k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Fixed effect inconsistency NMA BUGS script for count data with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Binomial_Binary_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      log(lambda[i,k]) <- mu[i] + d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for (k in 1:nt){ 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects 
    for (k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  } 
  B~dnorm(0,.001)%_%I(-4,4)
  
}







#' Fixed effect NMA BUGS script for count data with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Poisson_Continuous_Covariate <- function(){
  # Binomial likelihood, logit link, subgroup
  # Fixed effects model with one covariate
  for(i in 1:ns){# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)) {  
    for (k in (c+1):nt)  { 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for (k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Fixed effect NMA BUGS script for count data with adjustment for baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Poisson_Base_Risk <- function(){
  # Binomial likelihood, logit link, subgroup
  # Fixed effects model with one covariate
  for(i in 1:ns){# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  
  d[1] <- 0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4)
  
  
  #Output
  # pairwise OR, LORs, RR, RD for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)){  
    for (k in (c+1):nt)  { 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}







#' Fixed effect inconsistency NMA BUGS script for count data with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Poisson_Continuous_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      log(lambda[i,k]) <- mu[i] + d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for(k in 1:nt){ 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)) { # priors for all mean treatment effects 
    for(k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  } 
  B~dnorm(0,.001)%_%I(-4,4)
  
}





#' Fixed effect inconsistency NMA BUGS script for count data with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
FE_NMA_Poisson_Inconsist_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Fixed effects model
  for(i in 1:ns){					# LOOP THROUGH STUDIES
    mu[i]~dnorm(0,0.001)%_%I(-4,4)			# vague priors for all trial baselines
    for (k in 1:na[i]){				#  LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # event rate * exposure
      log(lambda[i,k]) <- mu[i] + d[t[i,1],t[i,k]] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]])		# summed residual deviance contribution for this trial
  }   
  totresdev <- sum(resdev[])			#Total Residual Deviance
  for(k in 1:nt){ 
    d[k,k] <- 0 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)) { # priors for all mean treatment effects 
    for(k in (c+1):nt) {d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  } 
  B~dnorm(0,.001)%_%I(-4,4)
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Random effect NMA BUGS script for count data
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson <- function(){ # *** PROGRAM STARTS
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      log(lambda[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k] ~ dnorm(0,.001)%_%I(-4,4)} # vague priors for treatment effects
  #sd ~ dunif(0,5) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 0.62)
  tau <- 1/exp(V)
  
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for (c in 1:(nt-1)) {  
    for (k in (c+1):nt)  { 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}





#' Random effect NMA BUGS script for count data (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Noninfo <- function(){ # *** PROGRAM STARTS
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      log(lambda[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) #Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  for (k in 2:nt){d[k] ~ dnorm(0,.001)%_%I(-4,4)} # vague priors for treatment effects
  sd ~ dunif(0,5) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt) { 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}





#' Random effect inconsistency NMA BUGS script for count data
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      log(lambda[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
}





#' Random effect inconsistency NMA BUGS script for count data (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Inconsist_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      log(lambda[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial variance
}






#' Random effect NMA BUGS script for count data with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Binary_Covariate <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  #sd~dunif(0,2) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
  
  #Output
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for (k in (c+1):nt){ 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}





#' Random effect NMA BUGS script for count data with a binary covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Binary_Covariate_Noninfo <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  sd~dunif(0,2) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  #Output
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt){ 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Random effect inconsistency NMA BUGS script for count data with a binary covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Binary_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]],tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for(k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
}






#' Random effect inconsistency NMA BUGS script for count data with a binary covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Binary_Covariate_Inconsist_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*x[i]
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]],tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  B~dnorm(0,.001)%_%I(-4,4)
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial variance
}







#' Random effect NMA BUGS script for count data with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Continuous_Covariate <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  #sd~dunif(0,2) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
  
  #Output
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt) { 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}





#' Random effect NMA BUGS script for count data with a continuous covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Continuous_Covariate_Noninfo <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  sd~dunif(0,2) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  #Output
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt){ 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Random effect NMA BUGS script for count data with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Base_Risk <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  #sd~dunif(0,2) # vague prior for between-trial SD
  #tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
  
  #Output
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt){ 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Random effect NMA BUGS script for count data with adjustment to baseline risk (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Base_Risk_Noninfo <- function(){
  # Binomial likelihood, logit link, subgroup
  # Random effects model for multi-arm trials
  for(i in 1:ns){# LOOP THROUGH STUDIES
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0 # treatment effect is zero for control arm
    mu[i]~dnorm(0,.001)%_%I(-4,4) # vague priors for all trial baselines
    for (k in 1:na[i]){# LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]){# LOOP THROUGH ARMS
      delta[i,k]~dnorm(md[i,k],taud[i,k])%_%I(-4,4) # trial-specific LOR distributions
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction)
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs
      sw[i,k] <- sum(w[i,1:(k-1)])/(k-1) # cumulative adjustment for multi-arm trials
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0 # covariate effect is zero for reference treatment
  for (k in 2:nt){# LOOP THROUGH TREATMENTS
    d[k]~dnorm(0,.001)%_%I(-4,4) # vague priors for treatment effects
    beta[k] <- B # common covariate effect
  }
  B~dnorm(0,.001)%_%I(-4,4) # vague prior for covariate effect
  sd~dunif(0,2) # vague prior for between-trial SD
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)
  
  #Output
  # pairwise RRs for all possible pair-wise comparisons, if nt>2
  for(c in 1:(nt-1)){  
    for(k in (c+1):nt){ 
      RR[c,k] <- exp(d[k]-d[c])
    }  
  }
  
  
  for(i in 1:ns){
    nt1[i] <- equals(t[i,1],1)
    mu1[i] <- mu[i]*nt1[i]
  }
  A <- sum(mu1[])/sum(nt1[])
  
  # Alternative: Given a Mean Effect, mean A, for 'standard' treatment 1, with precision (1/variance) precA
  #A~dnorm(meanA,precA)
  for(k in 1:nt){log(T[k]) <- A + d[k]}
}






#' Random effect inconsistency NMA BUGS script for count data with a continuous covariate
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Continuous_Covariate_Inconsist <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  B~dnorm(0,.001)%_%I(-4,4)
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
}





#' Random effect inconsistency NMA BUGS script for count data with a continuous covariate (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Continuous_Covariate_Inconsist_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for (c in 1:(nt-1)) { # priors for all mean treatment effects
    for (k in (c+1):nt) { d[c,k] ~ dnorm(0,.001)%_%I(-4,4) }
  }
  B~dnorm(0,.001)%_%I(-4,4)
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial variance
}






#' Random effect inconsistency NMA BUGS script for count data with adjustment to baseline risk
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Inconsist_Base_Risk <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for (k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  V ~ dnorm(-2.34, 0.62)
  tau <- exp(-V)
}






#' Random effect inconsistency NMA BUGS script for count data with adjustment to baseline risk (flat prior)
#'
#' Written to a temporary location as a '.bug' file and run from \code{\link{NMA_run_poisson_base_risk}}.
#'
#' @return None
#' @examples
#'
#' @export
RE_NMA_Poisson_Inconsist_Base_Risk_Noninfo <- function(){
  # Binomial likelihood, logit link
  # Random effects model for multi-arm trials
  for(i in 1:ns){ # LOOP THROUGH STUDIES
    delta[i,1]<-0 # treatment effect is zero in control arm
    mu[i] ~ dnorm(0,.001)%_%I(-4,4) # vague priors for trial baselines
    for (k in 1:na[i]) { # LOOP THROUGH ARMS
      r[i,k] ~ dpois(theta[i,k]) # Poisson likelihood
      theta[i,k] <- lambda[i,k]*E[i,k] # failure rate * exposure
      # model for linear predictor, covariate effect relative to treat in arm 1
      log(lambda[i,k]) <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)
      dev[i,k] <- 2*((theta[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/theta[i,k])) #Deviance contribution
    }
    resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial
    for (k in 2:na[i]) { # LOOP THROUGH ARMS
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau)%_%I(-4,4) # trial-specific LOR distributions
    }
  }
  totresdev <- sum(resdev[]) # Total Residual Deviance
  for (k in 1:nt) { 
    beta[k] <- B
  } 					# treatment effect is zero for reference treatment
  for(c in 1:(nt-1)){ # priors for all mean treatment effects
    for(k in (c+1):nt){d[c,k] ~ dnorm(0,.001)%_%I(-4,4)}
  }
  B~dnorm(0,.001)%_%I(-4,4)
  sd ~ dunif(0,2) # vague prior for between-trial standard deviation
  tau <- pow(sd,-2) # between-trial variance
}







#' Renumbering of treatments
#'
#' Renumbering of treatments/interventions in the evidence matrix successively.
#'
#' @param t The evidence matrix; each row stands for a single study.
#' @param trt_names Treatment names and numbers; a matrix consisting of 2 columns.
#'
#' @return A list with the following elements:
#' \item{t}{The renumbered evidence matrix.}
#' \item{trt_names}{The renumbered treatment names.}
#' @examples
#'
#' @export
renumbering = function(t, trt_names){
  exist = sort(unique(c(t[!is.na(t)])))
  for(i in 1:length(exist)){
    for(j in 1:nrow(t)){
      for(k in 1:ncol(t)){
        if(!is.na(t[j,k]) & t[j,k]==exist[i]) t[j,k] = i 
      }
    }
  }
  
  inds = which(trt_names[,2] %in% exist)
  trt_names = trt_names[inds,]
  row.names(trt_names) = c()
  return(list(t=t, trt_names=trt_names))
}



#' Handling zero cases
#'
#' Handling zero cases in studies with binary outcomes, by adding half a case to the \code{r} column and 1 case to the \code{n} column.
#'
#' @param r A matrix consisting of the number of cases found in the different treatment groups.
#' @param n A matrix consisting of the number of participants in the different treatment groups.
#'
#' @return A list with the following elements:
#' \item{r}{The adjusted cases matrix.}
#' \item{n}{The adjusted participants matrix.}
#' @examples
#'
#' @export
r_n_correct = function(r, n){
  inds = which(r==0)
  r[inds] = .5
  n[inds] = n[inds] + 1
  
  return(list(r=r, n=n))
}






#' Identifying covariate type
#'
#' Determining whether a covariate selected for meta-regression is binary or continuous.
#'
#' @param BUGS_data The original data set with headers in WinBUGS format.
#' @param Covar_Name The covariate name.
#'
#' @return A list with the following elements:
#' \item{BUGS_data}{The original data set with headers in WinBUGS format.}
#' \item{Type}{The covariate type: "Binary" or "Continuous".}
#' @examples
#'
#' @export
Covar_Type_Identify = function(BUGS_data, Covar_Name){
  ind = which(colnames(BUGS_data) == Covar_Name)
  v = as.numeric(BUGS_data[, ind])
  a = unique(v)
  
  if(length(a)==2){
    u = v
    u[v==a[1]] = 0
    u[v==a[2]] = 1
    
    BUGS_data[, ind] = u
    names(BUGS_data)[ind] = Covar_Name
    Type = "Binary"
  } else if(class(v) != "numeric"){
    Type = "Error"
  } else{
    Type = "Continuous"
  }
  
  return(list(Type=Type, data=BUGS_data))
}







#' Reading NMA data
#'
#' Reading the NMA evidence data, in order to pass it to WinBUGS.
#'
#' @param BUGS_data A text file consisting of columns entitled \code{t[,1]}, \code{t[,2]} etc. for directly compared treatments, as well as columns headed \code{r}'s and \code{n}'s (for binary outcomes) or \code{y}'s ans \code{se}'s (for continuous outcomes).
#' @param trt_names A matrix consisting of the different treatments/interventions and their numbers.
#' @param Correct Logical: do we or do we not wish to correct zero cases.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param ontinuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param response_type One of "binary", "binary_base_risk", "normal", "normal_base_risk", "nomal_SMD", "poisson" and "poisson_base_risk".
#' 
#' @return A list with the following elements:
#' \item{t}{The evidence matrix; each row stands for a single study.}
#' \item{r}{The (adjusted) number of events matrix (for binary or count outcomes).}
#' \item{E}{The exposure time matrix (for count outcomes).}
#' \item{n}{The adjusted number of patients matrix (for binary outcomes).}
#' \item{y}{The mean difference from baseline (for continuous outcomes).}
#' \item{se}{The standard error matrix (for continuous outcomes).}
#' \item{na}{A vector consisting of the numbers of treatments directly compared per study.}
#' \item{ns}{The overall number of studies.}
#' \item{nt}{The overall number of different treatments.}
#' \item{trt_names}{A matrix consisting of the different treatments/interventions and their numbers.}
#' \item{x}{Covariate values: \code{0/1} for binary or numeric for continuous.}
#' \item{mx}{Mean of covariate values (for continuous covariate values) or study baseline values (for baseline adjusted models).}
#' @examples
#'
#' @export
read_data = function(BUGS_data, trt_names, Correct, Binary_Covariate=NA, 
                     Continuous_Covariate=NA, response_type){
  if(!is.na(Continuous_Covariate) & !(Continuous_Covariate %in% names(BUGS_data))){
    cat(paste("\nError: covariate '", Continuous_Covariate, "' not found!\n\n", sep=''))
  } else if(!is.na(Binary_Covariate) & !(Binary_Covariate %in% names(BUGS_data))){
    cat(paste("\nError: covariate '", Binary_Covariate, "' not found!\n\n", sep=''))
  } else if(!is.na(Binary_Covariate) & !is.na(Continuous_Covariate)){
    cat("\nError: the software allows for one covariate at most!\n\n")
  } else{
    t_ind = grep("t[[]", colnames(BUGS_data))
    t = as.matrix(BUGS_data[,t_ind])
    temp = renumbering(t, trt_names)
    t = temp$t
    trt_names = temp$trt_names
    
    na = apply(t, 1, function(v) sum(!is.na(v)))
    ns = nrow(t)
    nt = length(unique(t[!is.na(t)]))
    
    if(response_type == "binomial"){
      r_ind = grep("r[[]", colnames(BUGS_data))
      r = as.matrix(BUGS_data[,r_ind])
      
      n_ind = grep("n[[]", colnames(BUGS_data))
      n = as.matrix(BUGS_data[,n_ind])
      
      if(Correct){
        temp = r_n_correct(r,n)
        r = temp$r
        n = temp$n
      }
      
      if(!is.na(Binary_Covariate)){
        x = as.vector(BUGS_data[,Binary_Covariate])
        return(list(t=t, r=r, n=n, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    x=x, Binary_Covariate=Binary_Covariate, Continuous_Covariate=NA))
      } else if(!is.na(Continuous_Covariate)){
        x = as.vector(BUGS_data[,Continuous_Covariate])
        x = (x-mean(x))/sd(x)
        mx = 0
        return(list(t=t, r=r, n=n, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    x=x, mx=mx, Binary_Covariate=NA, Continuous_Covariate=Continuous_Covariate))
      } else{
        return(list(t=t, r=r, n=n, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    Binary_Covariate=NA, Continuous_Covariate=NA))
      }
    } else if(response_type %in% c("normal", "normal_SMD")){
      y_ind = grep("y[[]", colnames(BUGS_data))
      se_ind = grep("se[[]", colnames(BUGS_data))
      y = as.data.frame(BUGS_data[,y_ind])
      se = as.data.frame(BUGS_data[,se_ind])
      
      if(response_type == "normal_SMD"){
        y = cbind(NA, y)
        names(y)[1] = "y[,1]"
        se = cbind(NA, se)
        names(se)[1] = "se[,1]"
      }
      
      y = as.matrix(y)
      se = as.matrix(se)
      
      if(!is.na(Binary_Covariate)){
        x = as.vector(BUGS_data[,Binary_Covariate])
        return(list(t=t, y=y, se=se, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    x=x, Binary_Covariate=Binary_Covariate, Continuous_Covariate=NA))
      } else if(!is.na(Continuous_Covariate)){
        x = as.vector(BUGS_data[,Continuous_Covariate])
        x = (x-mean(x))/sd(x)
        mx = 0
        return(list(t=t, y=y, se=se, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    x=x, mx=mx, Binary_Covariate=NA, Continuous_Covariate=Continuous_Covariate))
      } else{
        return(list(t=t, y=y, se=se, na=na, ns=ns, nt=nt, trt_names=trt_names,
                    Binary_Covariate=NA, Continuous_Covariate=NA))
      }
    } else if(response_type == "normal_base_risk"){
      y_ind = grep("y[[]", colnames(BUGS_data))
      se_ind = grep("se[[]", colnames(BUGS_data))
      y = as.matrix(BUGS_data[,y_ind])
      se = as.matrix(BUGS_data[,se_ind])
      one_inds = y[which(t==1)]
      mx = mean(one_inds)
      
      return(list(t=t, y=y, se=se, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                  mx=mx, Binary_Covariate=NA, Continuous_Covariate=NA))
    } else if(response_type == "binomial_base_risk"){
      r_ind = grep("r[[]", colnames(BUGS_data))
      r = as.matrix(BUGS_data[,r_ind])
      
      n_ind = grep("n[[]", colnames(BUGS_data))
      n = as.matrix(BUGS_data[,n_ind])
      
      if(Correct){
        temp = r_n_correct(r,n)
        r = temp$r
        n = temp$n
      }
      
      one_inds = which(t==1)
      mx = mean(log(r[one_inds]/(n[one_inds]-r[one_inds])))
      
      return(list(t=t, r=r, n=n, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                  mx=mx, Binary_Covariate=NA, Continuous_Covariate=NA))
    } else if(response_type == "poisson"){
      E_ind = grep("E[[]", colnames(BUGS_data))
      E = as.matrix(BUGS_data[,E_ind])
      
      r_ind = grep("r[[]", colnames(BUGS_data))
      r = as.matrix(BUGS_data[,r_ind])
      
      if(!is.na(Binary_Covariate)){
        x = as.vector(BUGS_data[,Binary_Covariate])
        return(list(t=t, E=E, r=r, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    x=x, Binary_Covariate=Binary_Covariate, Continuous_Covariate=NA))
      } else if(!is.na(Continuous_Covariate)){
        x = as.vector(BUGS_data[,Continuous_Covariate])
        x = (x-mean(x))/sd(x)
        mx = 0
        return(list(t=t, E=E, r=r, na=na, ns=ns, nt=nt, trt_names=trt_names, 
                    x=x, mx=mx, Binary_Covariate=NA, Continuous_Covariate=Continuous_Covariate))
      } else{
        return(list(t=t, E=E, r=r, na=na, ns=ns, nt=nt, trt_names=trt_names,
                    Binary_Covariate=NA, Continuous_Covariate=NA))
      }
    } else if(response_type == "poisson_base_risk"){
      E_ind = grep("E[[]", colnames(BUGS_data))
      E = as.matrix(BUGS_data[,E_ind])
      
      r_ind = grep("r[[]", colnames(BUGS_data))
      r = as.matrix(BUGS_data[,r_ind])
      
      one_inds = which(t==1)
      mx = mean(log(r[one_inds]/E[one_inds]))
      
      return(list(t=t, E=E, r=r, na=na, ns=ns, nt=nt, trt_names=trt_names,
                  mx=mx, Binary_Covariate=NA, Continuous_Covariate=NA))
    }
    
  }
}






#' Reading NMA data and identifying covariate type
#'
#' A wrapper function for \code{read_data} that identifies the meta-regressor type automatically.
#'
#' @param BUGS_data A text file consisting of columns entitled \code{t[,1]}, \code{t[,2]} etc. for directly compared treatments, as well as columns headed \code{r}'s and \code{n}'s (for binary outcomes) or \code{y}'s ans \code{se}'s (for continuous outcomes).
#' @param trt_names A matrix consisting of the different treatments/interventions and their numbers.
#' @param Correct Logical: do we or do we not wish to correct zero cases.
#' @param Covar_Name The covariate name for meta-regression (default is \code{NA}).
#' @param response_type One of "binary", "binary_base_risk", "normal",  "normal_base_risk", "nomal_SMD", "poisson" and "poisson_base_risk"..
#' 
#' @return The output of \code{read_data}.
#' 
#' @export
read_data_wrapper = function(BUGS_data, trt_names, Correct,
                             Covar_Name=NA, response_type){
  if(is.na(Covar_Name)){
    read_data(BUGS_data=BUGS_data, trt_names=trt_names, 
              Correct=Correct, response_type=response_type)
  } else{
    Binary_Covariate = Continuous_Covariate = NA
    data_and_covtype = Covar_Type_Identify(BUGS_data, Covar_Name)
    BUGS_data = data_and_covtype$data
    type = data_and_covtype$Type
    if(type == "Binary") Binary_Covariate = Covar_Name
    else if(type == "Continuous") Continuous_Covariate = Covar_Name
    
    read_data(BUGS_data=BUGS_data, trt_names=trt_names, Correct=Correct, 
              Binary_Covariate=Binary_Covariate, 
              Continuous_Covariate=Continuous_Covariate, 
              response_type=response_type)
  }
}







#' Running NMA for binary outcomes
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' @param Noninfo Logical: is a noninformative heterogeneity prior distribution used? Default is 'FALSE'.
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#' 
#' @export
NMA_run_binomial = function(RE, Binary_Covariate=NA, 
                            Continuous_Covariate=NA, 
                            data, Inconsist,
                            N_chains=3, N_iter=100000, 
                            burnin=50000, Sim_Size=50000, 
                            Noninfo=FALSE){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "rhat")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "rhat", "OR", "RR", "T")
  }
  
  if(RE){
    params = c(params, "tau")
    if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
      if(!Noninfo){
        inits = function(){
          list(d=d, mu=rep(0, ns), V=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Binomial_Inconsist.bug"),
                          file.path(tempdir(), "RE_NMA_Binomial.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Binomial_Inconsist,
                     RE_NMA_Binomial
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits,
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      } else{
        inits = function(){
          list(d=d, mu=rep(0, ns), sd=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Binomial_Inconsist_Noninfo.bug"),
                          file.path(tempdir(), "RE_NMA_Binomial_Noninfo.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Binomial_Inconsist_Noninfo,
                     RE_NMA_Binomial_Noninfo
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits,
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      }
    } else if(!is.na(Binary_Covariate)){
      params = c(params, "B")
      if(!Noninfo){
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, V=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Binomial_Binary_Covariate_Inconsist.bug"),
                          file.path(tempdir(), "RE_NMA_Binomial_Binary_Covariate.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Binomial_Binary_Covariate_Inconsist,
                     RE_NMA_Binomial_Binary_Covariate
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      } else{
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, sd=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Binomial_Binary_Covariate_Inconsist_Noninfo.bug"),
                          file.path(tempdir(), "RE_NMA_Binomial_Binary_Covariate_Noninfo.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Binomial_Binary_Covariate_Inconsist_Noninfo,
                     RE_NMA_Binomial_Binary_Covariate_Noninfo
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      }
    } else if(!is.na(Continuous_Covariate)){
      params = c(params, "B")
      if(!Noninfo){
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, V=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Binomial_Continuous_Covariate_Inconsist.bug"),
                          file.path(tempdir(), "RE_NMA_Binomial_Continuous_Covariate.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Binomial_Continuous_Covariate_Inconsist,
                     RE_NMA_Binomial_Continuous_Covariate
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      } else{
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, sd=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Binomial_Continuous_Covariate_Inconsist_Noninfo.bug"),
                          file.path(tempdir(), "RE_NMA_Binomial_Continuous_Covariate_Noninfo.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Binomial_Continuous_Covariate_Inconsist_Noninfo,
                     RE_NMA_Binomial_Continuous_Covariate_Noninfo
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      }
    } else if(!is.na(Binary_Covariate) & !is.na(continuous_Covariate)){
      cat("\nThe software allows for one covariate at most!\n\n")
    }
  } else if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
    inits = function(){
      list(d=d, mu=rep(0, ns))
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Binomial_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Binomial.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Binomial_Inconsist,
                 FE_NMA_Binomial
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Binary_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Binomial_Binary_Covariate_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Binomial_Binary_Covariate.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Binomial_Binary_Covariate_Inconsist,
                 FE_NMA_Binomial_Binary_Covariate
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Continuous_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Binomial_Continuous_Covariate_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Binomial_Continuous_Covariate.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Binomial_Continuous_Covariate_Inconsist,
                 FE_NMA_Binomial_Continuous_Covariate
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Binary_Covariate) & !is.na(continuous_Covariate)){
    cat("\nThe software allows for one covariate at most!\n\n")
  }
  
}






#' Running NMA for binary outcomes with adjustment for study baseline risk 
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' @param Noninfo Logical: is a noninformative heterogeneity prior distribution used? Default is 'FALSE'.
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#' 
#' @export
NMA_run_binomial_base_risk = function(RE, Binary_Covariate=NA, 
                                      Continuous_Covariate=NA, data, 
                                      Inconsist, N_chains=3, N_iter=100000, 
                                      burnin=50000, Sim_Size=50000, 
                                      Noninfo=FALSE){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "rhat")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "rhat", "OR", "RR", "T", "B")
  }
  
  if(RE){
    params = c(params, "tau")
    if(!Noninfo){
      inits = function(){
        list(d=d, mu=rep(0, ns), V=1, B=0)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Binomial_Inconsist_Base_Risk.bug"),
                        file.path(tempdir(), "RE_NMA_Binomial_Base_Risk.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Binomial_Inconsist_Base_Risk,
                   RE_NMA_Binomial_Base_Risk
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits, 
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    } else{
      inits = function(){
        list(d=d, mu=rep(0, ns), sd=1, B=0)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Binomial_Inconsist_Base_Risk_Noninfo.bug"),
                        file.path(tempdir(), "RE_NMA_Binomial_Base_Risk_Noninfo.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Binomial_Inconsist_Base_Risk_Noninfo,
                   RE_NMA_Binomial_Base_Risk_Noninfo
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits,
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    }
  } else{
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Binomial_Inconsist_Base_Risk.bug"),
                      file.path(tempdir(), "FE_NMA_Binomial_Base_Risk.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Binomial_Inconsist_Base_Risk,
                 FE_NMA_Binomial_Base_Risk
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  }
}







#' Running NMA for count data
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' @param Noninfo Logical: is a noninformative heterogeneity prior distribution used? Default is 'FALSE'.
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#' 
#' @export
NMA_run_poisson = function(RE, Binary_Covariate=NA, 
                           Continuous_Covariate=NA, 
                           data, Inconsist,
                           N_chains=3, N_iter=100000, 
                           burnin=50000, Sim_Size=50000, 
                           Noninfo=FALSE){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "theta")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "theta", "RR", "T")
  }
  
  if(RE){
    params = c(params, "tau")
    if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
      if(!Noninfo){
        inits = function(){
          list(d=d, mu=rep(0, ns), V=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Poisson_Inconsist.bug"),
                          file.path(tempdir(), "RE_NMA_Poisson.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Poisson_Inconsist,
                     RE_NMA_Poisson
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits,
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      } else{
        inits = function(){
          list(d=d, mu=rep(0, ns), sd=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Poisson_Inconsist_Noninfo.bug"),
                          file.path(tempdir(), "RE_NMA_Poisson_Noninfo.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Poisson_Inconsist_Noninfo,
                     RE_NMA_Poisson_Noninfo
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits,
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      }
    } else if(!is.na(Binary_Covariate)){
      params = c(params, "B")
      if(Noninfo){
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, V=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Poisson_Binary_Covariate_Inconsist.bug"),
                          file.path(tempdir(), "RE_NMA_Poisson_Binary_Covariate.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Poisson_Binary_Covariate_Inconsist,
                     RE_NMA_Poisson_Binary_Covariate
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      } else{
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, sd=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Poisson_Binary_Covariate_Inconsist_Noninfo.bug"),
                          file.path(tempdir(), "RE_NMA_Poisson_Binary_Covariate_Noninfo.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Poisson_Binary_Covariate_Inconsist_Noninfo,
                     RE_NMA_Poisson_Binary_Covariate_Noninfo
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      }
    } else if(!is.na(Continuous_Covariate)){
      params = c(params, "B")
      if(!Noninfo){
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, V=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Poisson_Continuous_Covariate_Inconsist.bug"),
                          file.path(tempdir(), "RE_NMA_Poisson_Continuous_Covariate.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Poisson_Continuous_Covariate_Inconsist,
                     RE_NMA_Poisson_Continuous_Covariate
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      } else{
        inits = function(){
          list(d=d, mu=rep(0, ns), B=0, sd=1)
        }
        
        filename = ifelse(Inconsist,
                          file.path(tempdir(), "RE_NMA_Poisson_Continuous_Covariate_Inconsist_Noninfo.bug"),
                          file.path(tempdir(), "RE_NMA_Poisson_Continuous_Covariate_Noninfo.bug")
        )
        
        mod = ifelse(Inconsist,
                     RE_NMA_Poisson_Continuous_Covariate_Inconsist_Noninfo,
                     RE_NMA_Poisson_Continuous_Covariate_Noninfo
        )
        write.model(mod, filename)
        
        meta.sim = bugs(data, inits, 
                        model.file = filename,
                        parameters = params,
                        n.chains=N_chains, n.iter=N_iter, 
                        n.burnin = burnin, n.sims=Sim_Size)
      }
    } else if(!is.na(Binary_Covariate) & !is.na(continuous_Covariate)){
      cat("\nThe software allows for one covariate at most!\n\n")
    }
  } else if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
    inits = function(){
      list(d=d, mu=rep(0, ns))
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Poisson_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Poisson.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Poisson_Inconsist,
                 FE_NMA_Poisson
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Binary_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Poisson_Binary_Covariate_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Poisson_Binary_Covariate.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Poisson_Binary_Covariate_Inconsist,
                 FE_NMA_Poisson_Binary_Covariate
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Continuous_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Poisson_Continuous_Covariate_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Poisson_Continuous_Covariate.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Poisson_Continuous_Covariate_Inconsist,
                 FE_NMA_Poisson_Continuous_Covariate
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Binary_Covariate) & !is.na(continuous_Covariate)){
    cat("\nThe software allows for one covariate at most!\n\n")
  }
  
}






#' Running NMA for count data
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' @param Noninfo Logical: is a noninformative heterogeneity prior distribution used? Default is 'FALSE'.
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#' 
#' @export
NMA_run_poisson_base_risk = function(RE, Binary_Covariate=NA, 
                                     Continuous_Covariate=NA, data, 
                                     Inconsist, N_chains=3, 
                                     N_iter=100000, burnin=50000, 
                                     Sim_Size=50000, Noninfo=FALSE){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "theta")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "theta", "RR", "T", "B")
  }
  
  if(RE){
    params = c(params, "tau")
    if(!Noninfo){
      inits = function(){
        list(d=d, mu=rep(0, ns), V=1, B=0)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Poisson_Inconsist_Base_Risk.bug"),
                        file.path(tempdir(), "RE_NMA_Poisson_Base_Risk.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Poisson_Inconsist_Base_Risk,
                   RE_NMA_Poisson_Base_Risk
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits,
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    } else{
      inits = function(){
        list(d=d, mu=rep(0, ns), sd=1, B=0)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Poisson_Inconsist_Base_Risk_Noninfo.bug"),
                        file.path(tempdir(), "RE_NMA_Poisson_Base_Risk_Noninfo.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Poisson_Inconsist_Base_Risk_Noninfo,
                   RE_NMA_Poisson_Base_Risk_Noninfo
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits,
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    }
  } else{
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Poisson_Inconsist_Base_Risk.bug"),
                      file.path(tempdir(), "FE_NMA_Poisson_Base_Risk.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Poisson_Inconsist_Base_Risk,
                 FE_NMA_Poisson_Base_Risk
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  }
}






#' Running NMA for continuous outcomes
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#'
#' @export
NMA_run_normal = function(RE, Binary_Covariate, 
                          Continuous_Covariate, 
                          data, Inconsist,
                          N_chains=3, N_iter=100000, 
                          burnin=50000, Sim_Size=50000){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "theta")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "D", "theta", "T")
  }
  
  if(RE){
    params = c(params, "tau")
    if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
      inits = function(){
        list(d=d, mu=rep(0, ns), sd=1)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Normal_Inconsist.bug"),
                        file.path(tempdir(), "RE_NMA_Normal.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Normal_Inconsist,
                   RE_NMA_Normal
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits,
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    } else if(!is.na(Binary_Covariate)){
      params = c(params, "B")
      inits = function(){
        list(d=d, mu=rep(0, ns), B=0, sd=1)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Normal_Binary_Covariate_Inconsist.bug"),
                        file.path(tempdir(), "RE_NMA_Normal_Binary_Covariate.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Normal_Binary_Covariate_Inconsist,
                   RE_NMA_Normal_Binary_Covariate
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits, 
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    } else if(!is.na(Continuous_Covariate)){
      params = c(params, "B")
      inits = function(){
        list(d=d, mu=rep(0, ns), B=0, sd=1)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Normal_Continuous_Covariate_Inconsist.bug"),
                        file.path(tempdir(), "RE_NMA_Normal_Continuous_Covariate.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Normal_Continuous_Covariate_Inconsist,
                   RE_NMA_Normal_Continuous_Covariate
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits, 
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims=Sim_Size)
    }
  } else if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
    inits = function(){
      list(d=d, mu=rep(0, ns))
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Normal.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Inconsist,
                 FE_NMA_Normal
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Binary_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Binary_Covariate_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Normal_Binary_Covariate.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Binary_Covariate_Inconsist,
                 FE_NMA_Normal_Binary_Covariate
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  } else if(!is.na(Continuous_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Continuous_Covariate_Inconsist.bug"),
                      file.path(tempdir(), "FE_NMA_Normal_Continuous_Covariate.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Continuous_Covariate_Inconsist,
                 FE_NMA_Normal_Continuous_Covariate
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims=Sim_Size)
  }
}







#' Running NMA for continuous outcomes in a SMD format
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#'
#' @export
NMA_run_normal_SMD = function(RE, Binary_Covariate, 
                              Continuous_Covariate, 
                              data, Inconsist,
                              N_chains=3, N_iter=100000, 
                              burnin=50000, Sim_Size=50000){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "theta")
    if(RE) params = c(params, "Var")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "D", "theta")
    if(RE) params = c(params, "Var")
  }
  
  if(RE){
    if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
      inits = function(){
        list(d=d, V=1)
      }
      
      filename = ifelse(Inconsist, 
                        file.path(tempdir(), "RE_NMA_Normal_Inconsist_SMD.bug"),
                        file.path(tempdir(), "RE_NMA_Normal_SMD.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Normal_Inconsist_SMD,
                   RE_NMA_Normal_SMD
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits,
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims = Sim_Size)
    } else if(!is.na(Binary_Covariate)){
      params = c(params, "B")
      inits = function(){
        list(d=d, B=0, V=1)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Normal_Binary_Covariate_Inconsist_SMD.bug"),
                        file.path(tempdir(), "RE_NMA_Normal_Binary_Covariate_SMD.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Normal_Binary_Covariate_Inconsist_SMD,
                   RE_NMA_Normal_Binary_Covariate_SMD
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits, 
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims = Sim_Size)
    } else if(!is.na(Continuous_Covariate)){
      params = c(params, "B")
      inits = function(){
        list(d=d, B=0, V=1)
      }
      
      filename = ifelse(Inconsist,
                        file.path(tempdir(), "RE_NMA_Normal_Continuous_Covariate_Inconsist_SMD.bug"),
                        file.path(tempdir(), "RE_NMA_Normal_Continuous_Covariate_SMD.bug")
      )
      
      mod = ifelse(Inconsist,
                   RE_NMA_Normal_Continuous_Covariate_Inconsist_SMD,
                   RE_NMA_Normal_Continuous_Covariate_SMD
      )
      write.model(mod, filename)
      
      meta.sim = bugs(data, inits, 
                      model.file = filename,
                      parameters = params,
                      n.chains=N_chains, n.iter=N_iter, 
                      n.burnin = burnin, n.sims = Sim_Size)
    }
  } else if(is.na(Binary_Covariate) & is.na(Continuous_Covariate)){
    inits = function(){
      list(d=as.vector(d))
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Inconsist_SMD.bug"),
                      file.path(tempdir(), "FE_NMA_Normal_SMD.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Inconsist_SMD,
                 FE_NMA_Normal_SMD
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains = N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims = Sim_Size)
  } else if(!is.na(Binary_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Binary_Covariate_Inconsist_SMD.bug"),
                      file.path(tempdir(), "FE_NMA_Normal_Binary_Covariate_SMD.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Binary_Covariate_Inconsist_SMD,
                 FE_NMA_Normal_Binary_Covariate_SMD
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains = N_chains, n.iter = N_iter, 
                    n.burnin = burnin, n.sims = Sim_Size)
  } else if(!is.na(Continuous_Covariate)){
    params = c(params, "B")
    inits = function(){
      list(d=d, B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Continuous_Covariate_Inconsist_SMD.bug"),
                      file.path(tempdir(), "FE_NMA_Normal_Continuous_Covariate_SMD.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Continuous_Covariate_Inconsist_SMD,
                 FE_NMA_Normal_Continuous_Covariate_SMD
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains = N_chains, n.iter = N_iter, 
                    n.burnin = burnin, n.sims = Sim_Size)
  }
}







#' Running NMA for continuous outcomes with adjustment for study baseline value
#'
#' Generating an MCMC sample, using WinBUGS.
#'
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate.
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate.
#' @param data a list of inputs delivered to WinBUGS.
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' 
#' @return A \code{\link[R2WinBUGS]{bugs}} object.
#'
#' @export
NMA_run_normal_base_risk = function(RE, Binary_Covariate, 
                                    Continuous_Covariate, 
                                    data, Inconsist,
                                    N_chains=3, N_iter=100000, 
                                    burnin=50000, Sim_Size=50000){
  nt = data$nt
  ns = data$ns
  
  if(Inconsist){
    d_mat = matrix(0, ncol=nt, nrow=nt)
    for(i in 1:nt){
      for(j in 1:i){
        d_mat[i,j] = NA
      }
    }
    d = structure(.Data=d_mat)
    params = c("d", "dev", "theta")
  } else{
    d=c(NA, rep(0, nt-1))
    params = c("d", "dev", "D", "theta", "T", "B")
  }
  
  if(RE){
    params = c(params, "tau")
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0, sd=1)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "RE_NMA_Normal_Inconsist_Base_Risk.bug"),
                      file.path(tempdir(), "RE_NMA_Normal_Base_Risk.bug")
    )
    
    mod = ifelse(Inconsist,
                 RE_NMA_Normal_Inconsist_Base_Risk,
                 RE_NMA_Normal_Base_Risk
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits, 
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, n.burnin = burnin, n.sims = Sim_Size)
  } else{
    inits = function(){
      list(d=d, mu=rep(0, ns), B=0)
    }
    
    filename = ifelse(Inconsist,
                      file.path(tempdir(), "FE_NMA_Normal_Inconsist_Base_Risk.bug"),
                      file.path(tempdir(), "FE_NMA_Normal_Base_Risk.bug")
    )
    
    mod = ifelse(Inconsist,
                 FE_NMA_Normal_Inconsist_Base_Risk,
                 FE_NMA_Normal_Base_Risk
    )
    write.model(mod, filename)
    
    meta.sim = bugs(data, inits,
                    model.file = filename,
                    parameters = params,
                    n.chains=N_chains, n.iter=N_iter, 
                    n.burnin = burnin, n.sims = Sim_Size)
  }
}








#' Running NMA
#'
#' Creating an NMA object
#'
#' @param RE Logical: use \code{TRUE} for a random effect model (default is \code{FALSE})
#' @param BUGS_data A text file consisting of columns entitled \code{t1}, \code{t2} etc. for directly compared treatments, as well as columns headed \code{r}'s and \code{n}'s (for binary outcomes) or \code{y}'s ans \code{se}'s (for continuous outcomes).
#' @param trt_names A matrix consisting of the different treatments/interventions and their numbers. 
#' @param Correct Logical: do we or do we not wish to correct zero cases (default is \code{TRUE}).
#' @param Binary_Covariate Logical: does the model contain meta-regression on a binary covariate (default is \code{FALSE}).
#' @param Continuous_Covariate Logical: does the model contain meta-regression on a continuous covariate (default is \code{FALSE}).
#' @param response_type One of "binary", "binary_base_risk", "normal",  "normal_base_risk", "nomal_SMD", "poisson" and "poisson_base_risk".
#' @param Inconsist Logical: 'TRUE' if an inconsistency model is assumed. Otherwise 'FALSE'.
#' @param N_chains number of Markov chains (default: 3).
#' @param N_iter number of total iterations per chain (including burn in; default: 100000).
#' @param burnin length of burn in, i.e. number of iterations to discard at the beginning (default: 50000).
#' @param Sim_Size The approximate number of simulations to keep after thinning (default: 50000).
#' @param Noninfo Logical: is a noninformative heterogeneity prior distribution used? Default is 'FALSE'.
#' 
#' @return A list with the following elements:
#' \item{meta.sim}{A \code{\link[R2WinBUGS]{bugs}} object.}
#' \item{r}{The (adjusted) number of events matrix (for binary or count outcomes).}
#' \item{n}{The adjusted participants matrix (for binary outcomes).}
#' \item{y}{The mean differences from baseline (for continuous outcomes).}
#' \item{se}{The standard errors matrix (for continuous outcomes).}
#' \item{E}{The exposure time matrix (for count outcomes).}
#' \item{nt}{The overall number of different treatments.}
#' \item{trt_names}{A vector consisting of the different treatments/interventions.}
#' \item{leverages}{See \code{\link{http://scharr.dept.shef.ac.uk/nicedsu/wp-content/uploads/sites/7/2017/05/TSD2-General-meta-analysis-corrected-2Sep2016v2.pdf}} for details.}
#' \item{deviance_residuals}{See \code{\link{http://scharr.dept.shef.ac.uk/nicedsu/wp-content/uploads/sites/7/2017/05/TSD2-General-meta-analysis-corrected-2Sep2016v2.pdf}} for details.}
#' \item{DIC}{The Deviance Information Criterion value for the fitted model.}
#' \item{response_type}{One of "binary", "binary_base_risk", "normal",  "normal_base_risk", "nomal_SMD", "poisson" and "poisson_base_risk".}
#' \item{RE}{Logical: use "TRUE" for a random effect model.}
#' \item{Thin}{The thinning factor for the MCMC sampling.}
#' \item{N_chains}{Number of Markov chains.}
#' \item{N_iter}{Number of total iterations per chain.}
#' \item{burnin}{Length of burn in, i.e. number of iterations to discard at the beginning.}
#' \item{Hetero_Prior}{'Informative' if an informative prior is employed for the random effect variance, otherwise 'Noninformative'. If a fixed effect model is used - ''.}
#' \item{Runtime}{The duration of the MCMC sampling in a hh:mm:ss format (returned as a string).}
#' \item{base_risk_adjust}{Logical: are treatment effects adjusted for study baseline values?}
#' \item{Noninfo}{Logical: is a noninformative heterogeneity prior distribution used?}
#' \item{BUGS_data}{A text file consisting of columns entitled \code{t1}, \code{t2} etc. for directly compared treatments, as well as columns headed \code{r}'s and \code{n}'s (for binary outcomes) or \code{y}'s ans \code{se}'s (for continuous outcomes).}
#' 
#' 
#' @examples
#' 
#' @export
NMA = function(RE=FALSE, BUGS_data, trt_names, 
               Correct=TRUE, Covar_Name=NA, response_type, 
               Inconsist=FALSE, N_chains=3, N_iter=100000, 
               burnin=50000, Sim_Size=50000,
               Noninfo=FALSE){
  Thin = round((N_iter - burnin)/Sim_Size)
  str = strsplit(response_type, "_")[[1]]
  base_risk_adjust = ifelse(length(str) >= 3, TRUE, FALSE)
  Hetero_Prior = ifelse(RE, "Noninformative", "")
  if((response_type %in% c("normal", "normal_base_risk")) & RE) Noninfo = TRUE
  if(RE & ((!Noninfo)|(response_type == "normal_SMD"))) Hetero_Prior = "Informative"
  
  
  data.read = read_data_wrapper(BUGS_data, trt_names, Correct, Covar_Name, response_type)
  t = as.matrix(data.read$t)
  na = as.vector(data.read$na)
  ns = as.numeric(data.read$ns)
  nt = as.numeric(data.read$nt)
  Binary_Covariate = data.read$Binary_Covariate
  Continuous_Covariate = data.read$Continuous_Covariate
  
  BUGS_data[, grep("t[[]", colnames(BUGS_data))] = t
  
  trt_names = as.vector(data.read$trt_names[,1])
  if(!is.na(Binary_Covariate) | !is.na(Continuous_Covariate)) x = data.read$x
  if(!is.na(Continuous_Covariate)) mx = data.read$mx
  
  if(response_type == "binomial"){
    r = data.read$r
    n = data.read$n
    
    data = list(ns=ns, nt=nt, na=na, t=t, r=r, n=n)
    if(!is.na(Binary_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, r=r, n=n, x=x)
    if(!is.na(Continuous_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, r=r, n=n, x=x, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_binomial(RE, Binary_Covariate, Continuous_Covariate, data,
                                Inconsist, N_chains, N_iter, burnin, Sim_Size, Noninfo)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    rt = rtime - 3600*hrs
    mnts = rt%/%60
    scds = round(rt - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    rhat_inds = grep("rhat", rownames(meta.sim$summary))
    rhat = sims[,rhat_inds]
    rhat.tilde = apply(rhat, 2, mean)
    
    r.flat = c(t(r))
    r.flat = r.flat[!is.na(r.flat)]
    
    n.flat = c(t(n))
    n.flat = n.flat[!is.na(n.flat)]
    
    devs.tilde = 2*(r.flat*log(r.flat/rhat.tilde) + (n.flat-r.flat)*log((n.flat-r.flat)/(n.flat-rhat.tilde)))
    devs.tilde[is.na(devs.tilde)] = 0
    leverages = devs.bar - devs.tilde
    signs = sign(r.flat - rhat.tilde)
    w = sqrt(devs.bar)*signs
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                r=r, n=n, nt=nt, Runtime=str,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin, 
                Hetero_Prior=Hetero_Prior,
                base_risk_adjust=base_risk_adjust,
                Noninfo=Noninfo))
    
  } else if(response_type == "binomial_base_risk"){
    r = data.read$r
    n = data.read$n
    mx = data.read$mx
    
    data = list(ns=ns, nt=nt, na=na, t=t, r=r, n=n, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_binomial_base_risk(RE, Binary_Covariate, Continuous_Covariate, data,
                                          Inconsist, N_chains, N_iter, burnin, Sim_Size, Noninfo)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    rt = rtime - 3600*hrs
    mnts = rt%/%60
    scds = round(rt - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    rhat_inds = grep("rhat", rownames(meta.sim$summary))
    rhat = sims[,rhat_inds]
    rhat.tilde = apply(rhat, 2, mean)
    
    r.flat = c(t(r))
    r.flat = r.flat[!is.na(r.flat)]
    
    n.flat = c(t(n))
    n.flat = n.flat[!is.na(n.flat)]
    
    devs.tilde = 2*(r.flat*log(r.flat/rhat.tilde) + (n.flat-r.flat)*log((n.flat-r.flat)/(n.flat-rhat.tilde)))
    devs.tilde[is.na(devs.tilde)] = 0
    leverages = devs.bar - devs.tilde
    signs = sign(r.flat - rhat.tilde)
    w = sqrt(devs.bar)*signs
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                r=r, n=n, nt=nt, Runtime=str,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin,
                base_risk_adjust=base_risk_adjust,
                Hetero_Prior=Hetero_Prior,
                Noninfo=Noninfo))
    
  } else if(response_type == "normal"){
    y = data.read$y
    se = data.read$se
    
    data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se)
    if(!is.na(Binary_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se, x=x)
    if(!is.na(Continuous_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se, x=x, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_normal(RE, Binary_Covariate, Continuous_Covariate, data,
                              Inconsist, N_chains, N_iter, burnin, Sim_Size)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    r = rtime - 3600*hrs
    mnts = r%/%60
    scds = round(r - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    theta_inds = grep("theta", rownames(meta.sim$summary))
    thetas = sims[,theta_inds]
    theta.tilde = apply(thetas, 2, mean)
    
    
    y.flat = c(t(y))
    y.flat = y.flat[!is.na(y.flat)]
    se.flat = c(t(se))
    se.flat = se.flat[!is.na(se.flat)]
    devs.tilde = (y.flat-theta.tilde)^2/se.flat^2
    leverages = devs.bar - devs.tilde
    w = (y.flat-theta.tilde)/se.flat
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                y=y, se=se, nt=nt,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC, 
                Runtime=str,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin,
                Hetero_Prior=Hetero_Prior,
                base_risk_adjust=base_risk_adjust,
                Noninfo=Noninfo))
  } else if(response_type == "normal_SMD"){
    y = data.read$y
    se = data.read$se
    
    data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se)
    if(!is.na(Binary_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se, x=x)
    if(!is.na(Continuous_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se, x=x, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_normal_SMD(RE, Binary_Covariate, Continuous_Covariate, data,
                                  Inconsist, N_chains, N_iter, burnin, Sim_Size)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    r = rtime - 3600*hrs
    mnts = r%/%60
    scds = round(r - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    theta_inds = grep("theta", rownames(meta.sim$summary))
    thetas = sims[,theta_inds]
    theta.tilde = apply(thetas, 2, mean)
    
    
    y.flat = c(t(y))
    y.flat = y.flat[!is.na(y.flat)]
    se.flat = c(t(se))
    se.flat = se.flat[!is.na(se.flat)]
    devs.tilde = (y.flat-theta.tilde)^2/se.flat^2
    leverages = devs.bar - devs.tilde
    w = (y.flat-theta.tilde)/se.flat
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                y=y, se=se, nt=nt,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC, 
                Runtime=str,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin,
                Hetero_Prior=Hetero_Prior,
                base_risk_adjust=base_risk_adjust,
                Noninfo=Noninfo))
  } else if(response_type == "normal_base_risk"){
    y = as.matrix(data.read$y)
    se = as.matrix(data.read$se)
    mx = as.numeric(data.read$mx)
    
    data = list(ns=ns, nt=nt, na=na, t=t, y=y, se=se, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_normal_base_risk(RE, Binary_Covariate, Continuous_Covariate, data,
                                        Inconsist, N_chains, N_iter, burnin, Sim_Size)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    r = rtime - 3600*hrs
    mnts = r%/%60
    scds = round(r - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    theta_inds = grep("theta", rownames(meta.sim$summary))
    thetas = sims[,theta_inds]
    theta.tilde = apply(thetas, 2, mean)
    
    
    y.flat = c(t(y))
    y.flat = y.flat[!is.na(y.flat)]
    se.flat = c(t(se))
    se.flat = se.flat[!is.na(se.flat)]
    devs.tilde = (y.flat-theta.tilde)^2/se.flat^2
    leverages = devs.bar - devs.tilde
    w = (y.flat-theta.tilde)/se.flat
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                y=y, se=se, nt=nt,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC, 
                Runtime=str,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin,
                Hetero_Prior=Hetero_Prior,
                base_risk_adjust=base_risk_adjust,
                Noninfo=Noninfo))
  } else if(response_type == "poisson"){
    E = data.read$E
    r = data.read$r
    
    data = list(ns=ns, nt=nt, na=na, t=t, E=E, r=r)
    if(!is.na(Binary_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, E=E, r=r, x=x)
    if(!is.na(Continuous_Covariate)) data = list(ns=ns, nt=nt, na=na, t=t, E=E, r=r, x=x, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_poisson(RE, Binary_Covariate, Continuous_Covariate, data,
                               Inconsist, N_chains, N_iter, burnin, Sim_Size, Noninfo)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    rt = rtime - 3600*hrs
    mnts = rt%/%60
    scds = round(rt - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    theta_inds = grep("theta", rownames(meta.sim$summary))
    thetas = sims[,theta_inds]
    theta.tilde = apply(thetas, 2, mean)
    
    
    r.flat = c(t(r))
    r.flat = r.flat[!is.na(r.flat)]
    
    devs.tilde = 2*(r.flat*log(r.flat/theta.tilde) + (theta.tilde-r.flat))
    devs.tilde[is.na(devs.tilde)] = 0
    leverages = devs.bar - devs.tilde
    signs = sign(r.flat - theta.tilde)
    w = sqrt(devs.bar)*signs
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                E=E, r=r, nt=nt,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC, 
                Runtime=str,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin,
                Hetero_Prior=Hetero_Prior,
                base_risk_adjust=base_risk_adjust,
                Noninfo=Noninfo))
  } else if(response_type == "poisson_base_risk"){
    E = data.read$E
    r = data.read$r
    mx = data.read$mx
    
    data = list(ns=ns, nt=nt, na=na, t=t, E=E, r=r, mx=mx)
    
    ptm = proc.time()
    meta.sim = NMA_run_poisson_base_risk(RE, Binary_Covariate, Continuous_Covariate, data,
                                         Inconsist, N_chains, N_iter, burnin, Sim_Size, Noninfo)
    runtime = proc.time() - ptm
    rtime = runtime[3]
    hrs = rtime%/%3600
    rt = rtime - 3600*hrs
    mnts = rt%/%60
    scds = round(rt - 60*mnts)
    
    hrs_str = ifelse(hrs == 0, "", paste(hrs, "hours,"))
    mnts_str = ifelse(mnts == 0, "", paste(mnts, "minutes"))
    scds_str = paste(scds, "seconds")
    str = paste(hrs_str, mnts_str, scds_str)
    
    sims = meta.sim$sims.matrix
    sims = sims[,-ncol(sims)]
    
    dev_inds = grep("dev", rownames(meta.sim$summary))
    dev_inds = dev_inds[-length(dev_inds)]
    devs = sims[,dev_inds]
    devs.bar = apply(devs, 2, mean)
    
    theta_inds = grep("theta", rownames(meta.sim$summary))
    thetas = sims[,theta_inds]
    theta.tilde = apply(thetas, 2, mean)
    
    
    r.flat = c(t(r))
    r.flat = r.flat[!is.na(r.flat)]
    
    devs.tilde = 2*(r.flat*log(r.flat/theta.tilde) + (theta.tilde-r.flat))
    devs.tilde[is.na(devs.tilde)] = 0
    leverages = devs.bar - devs.tilde
    signs = sign(r.flat - theta.tilde)
    w = sqrt(devs.bar)*signs
    DIC = sum(leverages) + sum(devs.bar)
    
    return(list(meta.sim=meta.sim, trt_names=trt_names, 
                E=E, r=r, nt=nt,
                leverages=leverages, 
                deviance_residuals=w, DIC=DIC, 
                Runtime=str,
                response_type=response_type, 
                RE=RE, BUGS_data=BUGS_data,
                Covar_Name=Covar_Name,
                N_chains=N_chains,
                N_iter=N_iter,
                burnin=burnin,
                Thin=Thin,
                Hetero_Prior=Hetero_Prior,
                base_risk_adjust=base_risk_adjust,
                Noninfo=Noninfo))
  }
}






#' Goodness of fit diagnostics
#'
#' Plots the leverages vs. the deviance residuals, to help detect outlier studies. Also prints suspected outliers to the screen.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param RE Logical: use "TRUE" for a random effect model.
#' @param Plot Logical: should a levergae plot be generated? 
#'
#' 
#' @return A list with the following elements:
#' \item{outliers}{A matrix whose rows contain the study IDs and treatment arm numbers of treatment arms whose contribution to the DIC exceeds 3.}
#' \item{extreme_outliers}{A matrix whose rows contain the study IDs and treatment arm numbers of treatment arms whose contribution to the DIC exceeds 4.}
#'
#' @export
residual_plot_new = function(sim, RE, Plot=FALSE){
  meta.sim = sim$meta.sim
  DIC = sim$DIC
  w = sim$deviance_residuals
  leverages = sim$leverages
  
  if(Plot){
    x = seq(-3, 3, length=100)
    y1 = 1 - x^2
    y2 = 2 - x^2
    y3 = 3 - x^2
    y4 = 4 - x^2
    ylims = c(0, max(max(y4), max(leverages))+.5)
    xlims = c(min(min(x), min(w)), max(max(x), max(w)))
    
    df = data.frame(unname(cbind(w, leverages)))
    names(df) = c("w", "leverages")
    
    df2 = data.frame(cbind(x, y1, y2, y3, y4))
    names(df2) = c("x", "y1", "y2", "y3", "y4")
    
    p = ggplot() +
      ggtitle(paste0("DIC = ", round(DIC,1), 
                     ",  Mean Res. Dev. = ", round(sum(w^2), 2))) +
      coord_cartesian(xlim = xlims, ylim=ylims) +
      geom_line(aes(df2$x[df2$y1>0], df2$y1[df2$y1>0]), 
                lwd=2, lty=1, col=2) +
      geom_line(aes(df2$x[df2$y2>0], df2$y2[df2$y2>0]), 
                lwd=2, lty=2, col='palegreen3') +
      geom_line(aes(df2$x[df2$y3>0], df2$y3[df2$y3>0]), 
                lwd=2, lty=3, col='slateblue') +
      geom_line(aes(df2$x[df2$y4>0], df2$y4[df2$y4>0]), 
                lwd=2, lty=4, col='deepskyblue1') +
      geom_point(aes(df$w, df$leverages), shape=21, stroke = 1.6,
                 col='dodgerblue1', bg="black", size=5) +
      theme(axis.text=element_text(size=22),
            axis.title=element_text(size=22, face="bold")) +
      xlab(expression(bold(w[ik]))) + ylab(expression(bold(leverage[ik]))) + 
      theme(axis.line = element_line(colour = "black"), 
            plot.title = element_text(size = 16, face = "bold")) +
      scale_x_continuous(expand=c(.0,0)) + 
      scale_y_continuous(expand=c(.02,0)) 
  }
  
  outliers = which((leverages+w^2)>3 & (leverages+w^2) <= 4)
  extreme_outliers = which((leverages+w^2) > 4)
  
  outliers_tabs = function(set){
    if(length(set)>0){
      parse = function(i){
        row_num = as.integer(strsplit(strsplit(unlist(names(set)), ",")[i][[1]][1], "[[]")[[1]][2])
        col_num = as.integer(strsplit(strsplit(unlist(names(set)), ",")[i][[1]][2], "[]]")[[1]][1])
        string = c(row_num, col_num)
        
        string
      }
      t(sapply(1:length(set), parse))
    }
  }
  
  outliers = outliers_tabs(outliers)
  extreme_outliers = outliers_tabs(extreme_outliers)
  
  if(Plot){
    return(p)
  } else{
    return(list(outliers=outliers, extreme_outliers=extreme_outliers))
  }
}






residual_plotly = function(sim, RE, Plot=FALSE, resid_tab){
  meta.sim = sim$meta.sim
  DIC = sim$DIC
  w = resid_tab$`Residual Deviance` = as.numeric(as.character(resid_tab$`Residual Deviance`))
  leverages = resid_tab$`Leverage` = as.numeric(as.character(resid_tab$`Leverage`))
  
  if(Plot){
    x = seq(-3, 3, length=100)
    y1 = 1 - x^2
    y2 = 2 - x^2
    y3 = 3 - x^2
    y4 = 4 - x^2
    ylims = c(0, max(max(y4), max(leverages))+.5)
    xlims = c(min(min(x), min(w)), max(max(x), max(w)))
    
    df = resid_tab
    df$strings = paste0("Study: ", df[,1], ", Arm: ", df$Arm)
    
    df2 = data.frame(cbind(x, y1, y2, y3, y4))
    names(df2) = c("x", "y1", "y2", "y3", "y4")
    
    p = ggplot() +
      ggtitle(paste0("DIC = ", round(DIC,1), 
                     ",  Mean Res. Dev. = ", round(sum(w^2), 2))) +
      coord_cartesian(xlim = xlims, ylim=ylims) +
      geom_line(aes(df2$x[df2$y1 > 0], df2$y1[df2$y1 > 0]), 
                lwd=1, lty=1, col=2) +
      geom_line(aes(df2$x[df2$y2 > 0], df2$y2[df2$y2 > 0]), 
                lwd=1, lty=2, col='palegreen3') +
      geom_line(aes(df2$x[df2$y3 > 0], df2$y3[df2$y3 > 0]), 
                lwd=1, lty=3, col='slateblue') +
      geom_line(aes(df2$x[df2$y4 > 0], df2$y4[df2$y4 > 0]), 
                lwd=1, lty=4, col='deepskyblue1') +
      geom_point(df, 
                 mapping = aes(`Residual Deviance`, 
                               Leverage, text = strings), 
                 shape=21, stroke = 1,
                 col='dodgerblue1', bg="black", size=3) +
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=18, face="bold")) +
      xlab("Residual deviance") + 
      ylab('Leverage') + 
      theme(axis.line = element_line(colour = "black"), 
            plot.title = element_text(size = 14, face = "bold")) +
      scale_x_continuous(expand=c(.0,0)) + 
      scale_y_continuous(expand=c(.02,0)) 
    
    p =  ggplotly(p, tooltip = "strings")
  }
  
  outliers = df$strings[which((leverages+w^2)>3 & (leverages+w^2) <= 4)]
  extreme_outliers = df$strings[which((leverages+w^2) > 4)]
  
  if(Plot){
    return(p)
  } else{
    return(list(outliers=outliers, extreme_outliers=extreme_outliers))
  }
}









#' Odds-Ratio cross table
#'
#' Creates the OR cross-table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A matrix containing posterior medians and credible intervals. "Significant" cells are marked by an asterisk.
#'
#' @export
OR_CrossTab = function(sim, coverage=.95, bottom){
  trt_names = sim$trt_names
  trt_names = gsub("_",  " ", trt_names)
  
  meta.sim = sim$meta.sim
  nt = sim$nt
  
  alpha = 1 - coverage 
  analysis = function(c){
    quantile(c, c(.5, alpha/2, 1-alpha/2))
  }
  
  sims = meta.sim$sims.matrix
  OR_inds = grep("OR", colnames(sims))
  OR = sims[, OR_inds]
  OR_results = unname(round(as.data.frame(t(apply(OR, 2, analysis))), 2))
  r_nm = row.names(OR_results)
  
  Cross_Table = matrix('', ncol=nt, nrow=nt)
  bottom_ind = which(trt_names == bottom)
  
  for(i in 1:nrow(OR_results)){
    all = unlist(OR_results[i,])
    nm = r_nm[i]
    str = strsplit(nm, ",")[[1]]
    
    left = as.numeric(strsplit(str[1], "\\[")[[1]][2])
    if(left > bottom_ind){
      left = left - 1
    } else if(left == bottom_ind){
      left = length(trt_names)
    }
    
    right = as.numeric(strsplit(str[2], "\\]")[[1]][1])
    if(right > bottom_ind){
      right = right - 1
    } else if(right == bottom_ind){
      right = length(trt_names)
    }
    
    if(left > right){
      temp = left
      left = right
      right = temp
      all = 1/all
      temp = all[2]
      all[2] = all[3]
      all[3] = temp
    }
    
    
    md = sprintf("%.2f", all[1])
    LL = sprintf("%.2f", all[2])
    UL = sprintf("%.2f", all[3])
    Int = paste("[", LL, ",", UL, "]", sep='')
    Cross_Table[right,left] = paste(md, Int)
    if((log(as.numeric(LL))*log(as.numeric(UL))) >= 0) 
      Cross_Table[right,left] = paste("* ", Cross_Table[right,left])
  }
  
  diag(Cross_Table) = c(trt_names[-bottom_ind], trt_names[bottom_ind])
  return(Cross_Table)
}



#' Relative Risk (or Rate-Ratio) cross table
#'
#' Creates the OR cross-table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A matrix containing posterior medians and credible intervals. "Significant" cells are marked by an asterisk.
#'
#' @export
RR_CrossTab = function(sim, coverage=.95, bottom){
  trt_names = sim$trt_names
  trt_names = gsub("_",  " ", trt_names)
  
  meta.sim = sim$meta.sim
  nt = sim$nt
  
  alpha = 1 - coverage 
  analysis = function(c){
    quantile(c, c(.5, alpha/2, 1-alpha/2))
  }
  
  sims = meta.sim$sims.matrix
  RR_inds = grep("RR", colnames(sims))
  RR = sims[, RR_inds]
  RR_results = unname(round(as.data.frame(t(apply(RR, 2, analysis))), 2))
  r_nm = row.names(RR_results)
  
  Cross_Table = matrix('', ncol=nt, nrow=nt)
  bottom_ind = which(trt_names == bottom)
  
  for(i in 1:nrow(RR_results)){
    all = unlist(RR_results[i,])
    nm = r_nm[i]
    str = strsplit(nm, ",")[[1]]
    
    left = as.numeric(strsplit(str[1], "\\[")[[1]][2])
    if(left > bottom_ind){
      left = left - 1
    } else if(left == bottom_ind){
      left = length(trt_names)
    }
    
    right = as.numeric(strsplit(str[2], "\\]")[[1]][1])
    if(right > bottom_ind){
      right = right - 1
    } else if(right == bottom_ind){
      right = length(trt_names)
    }
    
    if(left > right){
      temp = left
      left = right
      right = temp
      all = 1/all
      temp = all[2]
      all[2] = all[3]
      all[3] = temp
    }
    
    
    md = sprintf("%.2f", all[1])
    LL = sprintf("%.2f", all[2])
    UL = sprintf("%.2f", all[3])
    Int = paste("[", LL, ", ", UL, "]", sep='')
    Cross_Table[right,left] = paste(md, Int)
    if((log(as.numeric(LL))*log(as.numeric(UL))) >= 0) 
      Cross_Table[right,left] = paste("*", Cross_Table[right,left])
  }
  
  diag(Cross_Table) = c(trt_names[-bottom_ind], trt_names[bottom_ind])
  return(Cross_Table)
}





#' Mean Difference cross table
#'
#' Creates the OR cross-table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A matrix containing posterior medians and credible intervals. "Significant" cells are marked by an asterisk.
#'
#' @export
Diff_CrossTab = function(sim, coverage=.95, bottom){
  trt_names = sim$trt_names
  trt_names = gsub("_",  " ", trt_names)
  
  meta.sim = sim$meta.sim
  nt = sim$nt
  Cross_Table = matrix('', ncol=nt, nrow=nt)
  
  alpha = 1 - coverage
  analysis = function(c){
    quantile(c, c(.5, alpha/2, 1-alpha/2))
  }
  
  sims = meta.sim$sims.matrix
  D_inds = grep("D", colnames(sims))
  D = sims[, D_inds]
  D_results = unname(round(as.data.frame(t(apply(D, 2, analysis))), 2))
  r_nm = row.names(D_results)
  
  bottom_ind = which(trt_names == bottom)
  for(i in 1:nrow(D_results)){
    all = unlist(D_results[i,])
    nm = r_nm[i]
    str = strsplit(nm, ",")[[1]]
    
    left = as.numeric(strsplit(str[1], "\\[")[[1]][2])
    if(left > bottom_ind){
      left = left - 1
    } else if(left == bottom_ind){
      left = length(trt_names)
    }
    
    right = as.numeric(strsplit(str[2], "\\]")[[1]][1])
    if(right > bottom_ind){
      right = right - 1
    } else if(right == bottom_ind){
      right = length(trt_names)
    }
    
    if(left > right){
      temp = left
      left = right
      right = temp
      all = -all
      temp = all[2]
      all[2] = all[3]
      all[3] = temp
    }
    
    
    md = sprintf("%.2f", all[1])
    md = gsub("-0\\.00", "0\\.00", md)
    LL = sprintf("%.2f", all[2])
    UL = sprintf("%.2f", all[3])
    Int = paste("[", LL, ", ", UL, "]", sep='')
    Int = gsub("-0\\.00", "0\\.00", Int)
    Cross_Table[right,left] = paste(md, Int)
    if((as.numeric(LL)*as.numeric(UL)) >= 0) 
      Cross_Table[right,left] = paste("*", Cross_Table[right,left])
  }
  
  diag(Cross_Table) = c(trt_names[-bottom_ind], trt_names[bottom_ind])
  
  return(Cross_Table)
}




#' Odds-Ratio forest plot
#'
#' Creates the OR cross-table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param baseline reference treatment number.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A forest plot of all interventions vs. the selected reference intervention.
#'
#' @export
OR_forest = function(sim, baseline, coverage=.95, all_vs_one=T){
  trt_names = sim$trt_names
  trt_names = gsub("_",  " ", trt_names)
  
  meta.sim = sim$meta.sim
  nt = sim$nt
  RE = sim$RE
  
  alpha = 1 - coverage
  analysis = function(c){
    c(median(c), sd(c), quantile(c, c(.5, alpha/2, 1-alpha/2)))
  }
  
  sims = meta.sim$sims.matrix
  OR_inds = grep("OR", colnames(sims))
  OR = sims[, OR_inds]
  OR_results = unname(as.data.frame(t(apply(OR, 2, analysis))))
  names(OR_results) = c("Mean", "SD", "Median", "LL", "UL")
  inds_greater = grep(paste("OR[[]", baseline,",", sep=''), rownames(OR_results))
  inds_smaller = grep(paste(",", baseline, "[]]", sep=''), rownames(OR_results))
  temp = 1/OR_results[inds_smaller,]
  
  OR_results[inds_smaller,1] = temp[,1]
  OR_results[inds_smaller,4] = temp[,5]
  OR_results[inds_smaller,5] = temp[,4]
  inds = sort(c(inds_greater, inds_smaller))
  relevant = OR_results[inds,]
  
  if(!all_vs_one){
    D_temp = relevant
    D_temp[,1] = 1/relevant[,1]
    D_temp[,4] = 1/relevant[,5]
    D_temp[,5] = 1/relevant[,4]
    relevant = D_temp
  }
  
  signif = (as.numeric(relevant[,4])>=1 | as.numeric(relevant[,5])<=1)
  
  
  forest_data = data.frame(
    cbind(c(NA, relevant[ ,1]), 
          c(NA, relevant[, 4]), 
          c(NA, relevant[, 5])
    ))
  names(forest_data) = c("mean", "lower", "upper")
  
  R = range(forest_data[-1,2:3])
  len = ifelse(diff(R) > .5, 5, 3)
  u = c(.025, .05, seq(.1, 2, by=.1), 5, 10, 20)
  Llim = ifelse(length(which(R[1]-u >= 0))>1, max(which(R[1]-u >= 0))+1, 1)
  Rlim = ifelse(length(which(R[2]-u <= 0))>1, min(which(R[2]-u <= 0))-1, length(u))
  u_new = u[Llim:Rlim]
  my_vals = seq(range(u_new)[1], range(u_new)[2], length=len)
  attr(my_vals, "labels") = round(my_vals,1)
  
  str = paste("vs.", trt_names[baseline])
  if(!all_vs_one) str = paste(trt_names[baseline], "vs.")
  inds2 = c(1:nt)[-baseline]
  
  options(digits=2)
  tabletext = cbind(
    c(str, trt_names[inds2]),
    c(paste0("OR (", 100*coverage, "% interval)"), 
      sapply(1:nrow(relevant), 
             function(i){
               x = paste(sprintf("%.2f", relevant[i,1]), " (", 
                         sprintf("%.2f", relevant[i,4]), ", ", 
                         sprintf("%.2f", relevant[i,5]), ")", sep='')
               if(signif[i]) x = paste("*", x)
               x
             })))
  
  forestplot::forestplot(tabletext, boxsize = .1, graphwidth=unit(.2, "npc"),
                         forest_data, new_page = TRUE, cex=2,
                         is.summary = c(TRUE, rep(FALSE, length(inds2))),
                         xlog=TRUE, 
                         txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                                       gpar(fontfamily = "")),
                                          ticks = gpar(fontfamily = "", cex=.75),
                                          xlab  = gpar(fontfamily = "HersheySerif", cex = 1.25),
                                          cex = 1.25),
                         clip =c(u[Llim]*.9, u[Rlim]*1.1),
                         xticks=my_vals,
                         col=fpColors(box="royalblue",line="darkblue", summary="royalblue", 
                                      hrz_lines = "#444444"),
                         vertices = TRUE)
}




#' Relative Risk (or Rate Ratio) forest plot
#'
#' Creates the OR cross-table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param baseline reference treatment number.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A forest plot of all interventions vs. the selected reference intervention.
#'
#' @export
RR_forest = function(sim, baseline, coverage=.95, all_vs_one=T){
  trt_names = sim$trt_names
  trt_names = gsub("_",  " ", trt_names)
  
  meta.sim = sim$meta.sim
  nt = sim$nt
  RE = sim$RE
  
  alpha = 1 - coverage
  analysis = function(c){
    c(median(c), sd(c), quantile(c, c(.5, alpha/2, 1-alpha/2)))
  }
  
  sims = meta.sim$sims.matrix
  RR_inds = grep("RR", colnames(sims))
  RR = sims[, RR_inds]
  RR_results = unname(as.data.frame(t(apply(RR, 2, analysis))))
  names(RR_results) = c("Mean", "SD", "Median", "LL", "UL")
  inds_greater = grep(paste("RR[[]", baseline,",", sep=''), rownames(RR_results))
  inds_smaller = grep(paste(",", baseline, "[]]", sep=''), rownames(RR_results))
  temp = 1/RR_results[inds_smaller,]
  
  RR_results[inds_smaller,1] = temp[,1]
  RR_results[inds_smaller,4] = temp[,5]
  RR_results[inds_smaller,5] = temp[,4]
  inds = sort(c(inds_greater, inds_smaller))
  relevant = RR_results[inds,]
  signif = (as.numeric(relevant[,4])>=1 | as.numeric(relevant[,5])<=1)
  
  if(!all_vs_one){
    D_temp = relevant
    D_temp[,1] = 1/relevant[,1]
    D_temp[,4] = 1/relevant[,5]
    D_temp[,5] = 1/relevant[,4]
    relevant = D_temp
  }
  
  forest_data = data.frame(
    cbind(c(NA, relevant[ ,1]), 
          c(NA, relevant[, 4]), 
          c(NA, relevant[, 5])
    ))
  names(forest_data) = c("mean", "lower", "upper")
  
  R = range(forest_data[-1,2:3])
  len = ifelse(diff(R) > .5, 5, 3)
  u = c(.025, .05, seq(.1, 2, by=.1), 5, 10, 20)
  Llim = ifelse(length(which(R[1]-u >= 0))>1, max(which(R[1]-u >= 0))+1, 1)
  Rlim = ifelse(length(which(R[2]-u <= 0))>1, min(which(R[2]-u <= 0))-1, length(u))
  u_new = u[Llim:Rlim]
  my_vals = c(1, seq(range(u_new)[1], range(u_new)[2], length=len))
  attr(my_vals, "labels") = round(my_vals,1)
  
  str = paste("vs.", trt_names[baseline])
  if(!all_vs_one) str = paste(trt_names[baseline], "vs.")
  inds2 = c(1:nt)[-baseline]
  
  options(digits=2)
  tabletext = cbind(
    c(str, trt_names[inds2]),
    c(paste0("RR (", 100*coverage, "% interval)"), 
      sapply(1:nrow(relevant), 
             function(i){
               x = paste(sprintf("%.2f", relevant[i,1]), " (", 
                         sprintf("%.2f", relevant[i,4]), ", ", 
                         sprintf("%.2f", relevant[i,5]), ")", sep='')
               if(signif[i]) x = paste("*", x)
               x
             })))
  
  forestplot::forestplot(tabletext, boxsize = .1, graphwidth=unit(.2, "npc"),
                         forest_data, new_page = TRUE, cex=2,
                         is.summary = c(TRUE, rep(FALSE, length(inds2))),
                         xlog=TRUE, 
                         txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                                       gpar(fontfamily = "")),
                                          ticks = gpar(fontfamily = "", cex=1),
                                          xlab  = gpar(fontfamily = "HersheySerif", cex = 1.5),
                                          cex = 1.25),
                         clip = c(u[Llim]*.9, u[Rlim]*1.1),
                         xticks = my_vals,
                         col=fpColors(box="royalblue",line="darkblue", summary="royalblue", 
                                      hrz_lines = "#444444"),
                         vertices = TRUE)
}






#' Mean Difference forest plot
#'
#' Creates the OR cross-table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param baseline reference treatment number.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A forest plot of all interventions vs. the selected reference intervention.
#'
#' @export
Diff_forest = function(sim, baseline, coverage=.95, all_vs_one=T){
  trt_names = sim$trt_names
  #trt_names = gsub("-",  " ", trt_names)
  trt_names = gsub("_",  " ", trt_names)
  
  meta.sim = sim$meta.sim
  nt = sim$nt
  RE = sim$RE
  
  alpha = 1 - coverage
  analysis = function(c){
    c(median(c), sd(c), quantile(c, c(.5, alpha/2, 1-alpha/2)))
  }
  
  sims = meta.sim$sims.matrix
  D_inds = grep("D", colnames(sims))
  D = sims[, D_inds]
  D_results = unname(as.data.frame(t(apply(D, 2, analysis))))
  names(D_results) = c("Mean", "SD", "Median", "LL", "UL")
  inds_greater = grep(paste("D[[]", baseline,",", sep=''), rownames(D_results))
  inds_smaller = grep(paste(",", baseline, "[]]", sep=''), rownames(D_results))
  temp = -D_results[inds_smaller,]
  
  D_results[inds_smaller,1] = temp[,1]
  D_results[inds_smaller,4] = temp[,5]
  D_results[inds_smaller,5] = temp[,4]
  inds = sort(c(inds_greater, inds_smaller))
  relevant = D_results[inds,]
  signif = ((as.numeric(relevant[,4])*as.numeric(relevant[,5])) >= 0)
  
  if(!all_vs_one){
    D_temp = relevant
    D_temp[,1] = -relevant[,1]
    D_temp[,4] = -relevant[,5]
    D_temp[,5] = -relevant[,4]
    relevant = D_temp
  }
  
  forest_data = data.frame(
    cbind(c(NA, relevant[ ,1]), 
          c(NA, relevant[, 4]), 
          c(NA, relevant[, 5])
    ))
  names(forest_data) = c("mean", "lower", "upper")
  
  U_sort = sort(forest_data$upper[-1], 
                decreasing=TRUE)
  L_sort = sort(forest_data$lower[-1])
  L_lim = L_sort[2]
  U_lim = U_sort[2]
  
  str = paste("vs.", trt_names[baseline])
  if(!all_vs_one) str = paste(trt_names[baseline], "vs.")
  inds2 = c(1:nt)[-baseline]
  delta = L_lim - U_lim
  
  options(digits=2)
  tabletext = cbind(
    c(str, trt_names[inds2]),
    c(paste0("Mean Diff. (", 100*coverage, "% interval)"), 
      sapply(1:nrow(relevant), 
             function(i){
               x = paste(sprintf("%.2f", relevant[i,1]), " (", 
                         sprintf("%.2f", relevant[i,4]), ", ", 
                         sprintf("%.2f", relevant[i,5]), ")", sep='')
               if(signif[i]) x = paste("*", x)
               x
             })))
  tabletext = gsub("-0.00", "0.00", tabletext)
  
  forestplot::forestplot(tabletext, boxsize = .1, graphwidth=unit(.2, "npc"),
                         forest_data, new_page = TRUE, cex=2.5,
                         is.summary = c(TRUE, rep(FALSE, length(inds2))),
                         xlog=FALSE, 
                         clip =c(L_lim - .1*delta, U_lim + .1*delta),
                         txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),
                                                       gpar(fontfamily = "")),
                                          ticks = gpar(fontfamily = "", cex=1),
                                          xlab  = gpar(fontfamily = "HersheySerif", cex = 1.25),
                                          cex = 1.25),
                         xticks=round(c(0, seq(L_lim, U_lim, length=5)),1),
                         col=fpColors(box="royalblue",line="darkblue", summary="royalblue", 
                                      hrz_lines = "#444444"),
                         vertices = TRUE)
}






#' Plotting the network
#'
#' Plots the intervention network and return the edges matrix.
#'
#' @param BUGS_data A text file consisting of columns entitled \code{t1}, \code{t2} etc. for directly compared treatments, as well as columns headed \code{r}'s and \code{n}'s (for binary outcomes) or \code{y}'s ans \code{se}'s (for continuous outcomes).
#' @param Plot Logical: should a levergae plot be generated? Default is \code{TRUE}.
#' @param lay The layout used for the plot. Options include "circ" (circle), "layered" (tree), "kk" (Kamada-Kawai), "star" and "LGL" (large graph layout). Default is "circ".
#' 
#' @return A data frame whose rows consist of the different treatments and the number of studies connecting them.
#'
#' @export
NetworkPlot = function(BUGS_data, Plot=TRUE, lay="circ"){
  t_ind = grep("t[[]", colnames(BUGS_data))
  T_table = as.matrix(BUGS_data[,t_ind])
  
  aux_func = function(i){
    a = T_table[i,]
    combn(a[!is.na(a)], 2, simplify=FALSE)
  }
  M = as.data.frame(matrix(unlist(sapply(1:nrow(T_table), aux_func)), ncol=2, byrow=TRUE))
  x = rep(1, nrow(M))
  edges = data.frame(aggregate(x  ~ M[,1]+M[,2], sum, data=M))
  edges[,3] = as.integer(edges[,3])
  edges = edges[with(edges, order(edges[,1], edges[,2])),]
  names(edges) = c("Tx 1", "Tx 2", "No. of Studies")
  row.names(edges) = c()
  vertices = sort(unique(unlist(M)))
  
  net = igraph::graph_from_data_frame(d=edges, vertices=vertices, directed=F)
  E(net)$width = edges[,3]
  
  
  if(Plot & (lay == "circ")){
    layout = layout.circle(net)
  } else if(Plot & (lay == "star")){
    layout = layout_as_star(net)
  }else if(Plot & (lay == "layered")){
    layout = layout_with_sugiyama(net)
  }else if(Plot & (lay == "grid")){
    layout = layout_on_grid(net)
  } else if(Plot & (lay == "kk")){
    layout = layout_with_kk(net)
  }
  
  if(Plot){
    par(mar=c(1,1,1,1), mfrow=c(1,1))
    if(lay == "layered"){
      plot(layout$extd_graph, edge.width=E(net)$width)
    } else{
      plot(net, layout=layout)
    }
  }
  
  edges
}





#' Plotting the network
#'
#' Plots the intervention network and return the edges matrix.
#'
#' @param BUGS_data A text file consisting of columns entitled \code{t1}, \code{t2} etc. for directly compared treatments, as well as columns headed \code{r}'s and \code{n}'s (for binary outcomes) or \code{y}'s ans \code{se}'s (for continuous outcomes).
#' @param Plot Logical: should a levergae plot be generated? Default is \code{TRUE}.
#' @param lay The layout used for the plot. Options include "circ" (circle), "layered" (tree), "kk" (Kamada-Kawai), "star" and "LGL" (large graph layout). Default is "circ".
#' @param vertex.color A vector of colours for the different nodes. Entry options include "Orange", "Light Pink", "Deep Pink", "Brown", "Light Yellow", "Black" and "Blue". 
#' 
#' @return A data frame whose rows consist of the different treatments and the number of studies connecting them.
#'
#' @export
NetworkPlot_new = function(BUGS_data, Plot=TRUE, lay="circ", vertex.color){
  mycircle = function(coords, v=NULL, params) {
    vertex.color = params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color = vertex.color[v]
    }
    vertex.size  = 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size = vertex.size[v]
    }
    vertex.frame.color = params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color = vertex.frame.color[v]
    }
    vertex.frame.width = params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width = vertex.frame.width[v]
    }
    
    mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
           vertex.size, vertex.frame.width,
           FUN=function(x, y, bg, fg, size, lwd) {
             symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                     circles=size, add=TRUE, inches=FALSE)
           })
  }
  
  
  add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                   plot=mycircle, 
                   parameters=list(vertex.frame.color='black',
                                   vertex.frame.width=1))
  
  
  t_ind = grep("t[[]", colnames(BUGS_data))
  T_table = as.matrix(BUGS_data[,t_ind])
  
  n_trts = length(unique(na.omit(c(T_table))))
  vertex.label.color= rep("black", n_trts)
  vertex.label.color[vertex.color %in% c("saddlebrown", "black", "deeppink4")] = "white"
  
  aux_func = function(i){
    a = T_table[i,]
    combn(a[!is.na(a)], 2, simplify=FALSE)
  }
  M = as.data.frame(matrix(unlist(sapply(1:nrow(T_table), aux_func)), ncol=2, byrow=TRUE))
  x = rep(1, nrow(M))
  edges = data.frame(aggregate(x  ~ M[,1]+M[,2], sum, data=M))
  edges[,3] = as.integer(edges[,3])
  edges = edges[with(edges, order(edges[,1], edges[,2])),]
  names(edges) = c("Tx 1", "Tx 2", "No. of Studies")
  row.names(edges) = c()
  vertices = sort(unique(unlist(M)))
  
  net = igraph::graph_from_data_frame(d=edges, vertices=vertices, directed=F)
  E(net)$width = edges[,3]
  
  if(Plot & (lay == "circ")){
    layout = layout.circle(net)
  } else if(Plot & (lay == "star")){
    layout = layout_as_star(net)
  }else if(Plot & (lay == "layered")){
    layout = layout_with_sugiyama(net)
  }else if(Plot & (lay == "grid")){
    layout = layout_on_grid(net)
  } else if(Plot & (lay == "kk")){
    layout = layout_with_kk(net)
  } 
  
  if(Plot){
    par(mar=c(1,1,1,1), mfrow=c(1,1))
    if(lay == "layered"){
      plot(layout$extd_graph, 
           edge.width = E(net)$width, 
           vertex.label.color = vertex.label.color,
           vertex.color = vertex.color,
           vertex.label.font = 2,
           vertex.shape = "fcircle",
           vertex.frame.width = 2)
    } else if(lay == "lgl"){
      plot(net,
           edge.width = E(net)$width, 
           vertex.label.color = vertex.label.color,
           vertex.color = vertex.color,
           vertex.label.font = 2, 
           layout=layout_with_lgl,
           vertex.shape = "fcircle",
           vertex.frame.width = 2)
    } else{
      plot(net, layout = layout, 
           vertex.label.color = vertex.label.color,
           vertex.color = vertex.color,
           vertex.label.font = 2,
           vertex.shape = "fcircle",
           vertex.frame.width = 2)
    }
  }
  
  edges
}








#' MCMC diagnostic plots
#'
#' MCMC convergence diagnostics and statistics for a single parameter.
#'
#' @param MCMC_obj A A \code{\link[R2WinBUGS]{bugs}} object converted to a Coda format (by as.mcmc).
#' @param param One of the model basic parameters.
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' 
#' @return 
#'
#' @export
diagnosis_new = function(MCMC_obj, param, sim){
  if(grepl("trt. effect", param)){
    trt = strsplit(param, " trt. effect")[[1]][1]
    ind = which(sim$trt_names == trt)
    param = paste0("d[", ind, "]")
  } else if(param == "Precision"){
    param = "tau"
  } else if(param == "Regression coefficient"){
    param = "B"
  } else if(param == "Variance"){
    param = "Var"
  }
  
  mcmcplot1(MCMC_obj[, param, drop=FALSE], greek=FALSE)
}





Inconsistency_ggplot = function(consist_sim, inconsist_sim){
  inconsist_inds = grep("dev[[]", colnames(inconsist_sim$meta.sim$sims.matrix))
  inconsist_devs = inconsist_sim$meta.sim$sims.matrix[, inconsist_inds]
  post_mean_devs_inconsist = unname(apply(inconsist_devs, 2, mean))
  
  consist_inds = grep("dev[[]", colnames(consist_sim$meta.sim$sims.matrix))
  consist_devs = consist_sim$meta.sim$sims.matrix[, consist_inds]
  post_mean_devs_consist = unname(apply(consist_devs, 2, mean))
  
  col_names = colnames(inconsist_devs)
  aux_func = function(x){
    study = strsplit(x, split=",")[[1]][1]
    study = strsplit(study, split="[[]")[[1]][2]
    arm = strsplit(x, split=",")[[1]][2]
    arm = strsplit(arm, split="[]]")[[1]][1]
    
    paste("study", study, "arm", arm)
  }
  strings = unname(sapply(col_names, aux_func))
  
  df = as.data.frame(cbind(post_mean_devs_consist, 
                           post_mean_devs_inconsist,
                           strings)
  )
  df$post_mean_devs_consist = as.numeric(as.character(df$post_mean_devs_consist))
  df$post_mean_devs_inconsist = as.numeric(as.character(df$post_mean_devs_inconsist))
  xlim_low = min(min(df$post_mean_devs_consist), 0)
  ylim_low = min(min(df$post_mean_devs_inconsist), 0)
  xlim_up = max(max(df$post_mean_devs_consist), 2)
  ylim_up = max(max(df$post_mean_devs_inconsist), 2)
  
  p = ggplot(data = df, aes(x = post_mean_devs_consist, 
                            y = post_mean_devs_inconsist,
                            text=strings)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", size=1.5) + 
    geom_point(col="blue", size=3) + 
    scale_x_continuous(name="Consistency Model Deviances",
                       limits=c(xlim_low, xlim_up)) +
    scale_y_continuous(name="Inconsistency Model Deviances",
                       limits=c(ylim_low, ylim_up)) +
    theme(axis.title=element_text(size=16,face="bold"))
  
  p
}







#' NMA consistency diagnostic plot
#'
#' Consistency vs. inconsistency residual deviance interactive plot.
#'
#' @param consist_sim An NMA object: the output of the \code{\link{NMA}} function for the consistency model.
#' @param inconsist_sim An NMA object: the output of the \code{\link{NMA}} function for the inconsistency model.
#' 
#' @return 
#'
#' @export
Inconsistency_plot = function(p){
  return(ggplotly(p, tooltip = "strings"))
}





#' NMA consistency posterior comparison table
#'
#' Consistency vs. inconsistency parameter posterior comparison.
#'
#' @param consist_sim An NMA object: the output of the \code{\link{NMA}} function for the consistency model.
#' @param inconsist_sim An NMA object: the output of the \code{\link{NMA}} function for the inconsistency model.
#' 
#' @return A table comparing the consistency and the inconsistency models by the posterior medians and 95% CrI of each parameter.
#'
#' @export
Inconsist_table = function(consist_sim, inconsist_sim){
  inconsist_sum = inconsist_sim$meta.sim$summary
  inds = grep("d[[]", row.names(inconsist_sum))
  tab = round(inconsist_sum[inds, c(5,3,7)],2)
  
  aux_func = function(x){
    paste0(sprintf("%.2f", x[1]),
           " [", sprintf("%.2f", x[2]),
           ", ", sprintf("%.2f", x[3]), "]")
  }
  
  inconsist_vals = apply(tab, 1, aux_func)
  
  analysis = function(c){
    quantile(c, c(.5, .025, .975))
  }
  
  sims = consist_sim$meta.sim$sims.matrix
  if(consist_sim$response_type %in% c("binomial", "binomial_base_risk")){
    OR_inds = grep("OR", colnames(sims))
    OR = sims[, OR_inds]
    D = round(log(unname(round(as.data.frame(t(apply(OR, 2, analysis))), 2))),2)
  } else if(consist_sim$response_type %in% c("poisson", "poisson_base_risk")){
    RR_inds = grep("RR", colnames(sims))
    RR = sims[, RR_inds]
    D = round(log(unname(round(as.data.frame(t(apply(RR, 2, analysis))), 2))),2)
  } else{
    D_inds = grep("D[[]", colnames(sims))
    D_mat = sims[, D_inds]
    D = round(unname(round(as.data.frame(t(apply(D_mat, 2, analysis))), 2)),2)
  }
  consist_vals = apply(D, 1, aux_func)
  
  table = as.data.frame(cbind(consist_vals, inconsist_vals))
  names(table) = c("Consistency", "Inconsistency")
  Parameter = names(inconsist_vals)
  table = cbind(Parameter, table)
  row.names(table) = c()
  
  signif = function(u){
    x = u[2:3]
    y = u[5:6]
    min(x) > max(y) | min(y) > max(x)
  }
  
  big_diffs = apply(cbind(tab, D), 1, signif)
  for(i in 1:nrow(table)){
    if(big_diffs[i])
      table$Parameter[i] = paste0(table$Parameter[i], "*")
  }
  
  table
}





#' NMA consistency diagnostic table
#'
#' Consistency vs. inconsistency goodness of fit comparison.
#'
#' @param consist_sim An NMA object: the output of the \code{\link{NMA}} function for the consistency model.
#' @param inconsist_sim An NMA object: the output of the \code{\link{NMA}} function for the inconsistency model.
#' 
#' @return A table comparing the consistency and the inconsistency models by residual deviance, effective number of degrees of freedom and DIC.
#'
Inconsist_DIC_Summary = function(consist_sim, inconsist_sim){
  sims = consist_sim$meta.sim$sims.matrix
  dev_inds_consist = grep("dev[[]", rownames(consist_sim$meta.sim$summary))
  devs_consist = sims[,dev_inds_consist]
  devs.bar_consist = apply(devs_consist, 2, mean)
  dev_consist = sum(devs.bar_consist)
  
  sims = inconsist_sim$meta.sim$sims.matrix
  dev_inds_inconsist = grep("dev[[]", rownames(inconsist_sim$meta.sim$summary))
  devs_inconsist = sims[,dev_inds_inconsist]
  devs.bar_inconsist = apply(devs_inconsist, 2, mean)
  dev_inconsist = sum(devs.bar_inconsist)
  
  Deviance = c(dev_consist, dev_inconsist)
  
  pD_consist = consist_sim$DIC - dev_consist
  pD_inconsist = inconsist_sim$DIC - dev_inconsist
  pD = c(pD_consist, pD_inconsist)
  
  DIC = c(consist_sim$DIC, inconsist_sim$DIC)
  Quantity = c("Deviance", "pD", "DIC")
  sum_tab = as.data.frame(cbind(Quantity, t(round(cbind(Deviance, pD, DIC),1))))
  row.names(sum_tab) = c()
  names(sum_tab) = c("", "Consistency", "Inconsistency")
  
  sum_tab
}






#' Absolute effects table
#'
#' Absolute effects plot and summary table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param coverage Coverage of Bayesian credible interval (default is \code{0.95}).
#' 
#' @return A list with the following elements:
#' \item{p}{A ggplot2 object to be printed. A  bar plot with error bars denoting credible limits.}
#' \item{Table}{A summary table containing the different absolute effects and their credible intervals.}
#'
#' @export
Abs_Effects = function(sim, coverage){
  trt_names = gsub("_", " ", sim$trt_names)
  M = sim$meta.sim$sims.matrix
  trt_effs = M[, grep("T[[]", colnames(M))]
  alpha = 1 - coverage
  analysis = function(x){
    quantile(x, c(.5, alpha/2, 1-alpha/2))
  }
  stats = as.data.frame(t(apply(trt_effs, 2, analysis)))
  row.names(stats) = c()
  names(stats) = c("Est", "LL", "UL")
  nt = length(trt_names)
  stats$Tx = 1:nt
  
  
  p = ggplot(stats, aes(x=Tx, y=Est)) + 
    geom_bar(width = 0.5, 
             stat="identity",
             fill="steelblue") +
    scale_x_continuous(breaks = 1:nt, labels=c(trt_names)) +
    labs(title = paste("Absolute Effects"),
         x = "Treatment Arm",
         y = "Effect") + 
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18, face="bold"),
          plot.title = element_text(size=22)) +
    geom_errorbar(aes(ymax=UL, 
                      ymin=LL), size=1, width=.25) + 
    theme(axis.text.x = element_text(angle=45, hjust=1, size=14)) + 
    expand_limits(x = 1)
  
  concat = function(x){
    paste0(sprintf("%.2f", as.numeric(x[1])), 
           " [", 
           sprintf("%.2f", as.numeric(x[2])),
           ", ",
           sprintf("%.2f", as.numeric(x[3])),
           "]")
  }
  
  CrIs = apply(stats, 1, concat)
  tab = as.data.frame(cbind(trt_names, CrIs))
  names(tab) = c("Treatment Arm",
                 paste0("Effect (", coverage*100, "% CrI)"))
  
  return(list(p=p, Table=tab))
}






#' Residual deviance table
#'
#' Table of residual deviances and leverages.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#'
#' 
#' @return A table of the different treatment arms, along with their residual deviances and leverages.
#'
#' @export
RD_table = function(sim){
  BUGS_data = sim$BUGS_data
  
  t = BUGS_data[,grep("t[[]", colnames(BUGS_data))]
  na = apply(t, 1, function(v) sum(!is.na(v)))
  w = sim$deviance_residuals
  l = sim$leverages
  nms = names(w)
  inds = grep("ID", colnames(BUGS_data))
  IDs = as.matrix(BUGS_data[, inds])
  colnames(IDs) = colnames(BUGS_data)[inds]
  if(ncol(IDs) > 0){
    ID_nms = colnames(BUGS_data)[inds]
    IDs = apply(as.matrix(IDs), 2, function(x){gsub("_", " ", x)})
  } 
  
  rd_line = function(i){
    rd = w[i]
    lev = l[i]
    nm = nms[i]
    temp = strsplit(nm, ",")[[1]]
    study = as.numeric(strsplit(temp[1], "[[]")[[1]][2])
    arm = as.numeric(strsplit(temp[2], "[]]")[[1]][1])
    Tx = as.character(sim$trt_names[t[study,arm]])
    
    if(ncol(IDs) > 0){
      ID = IDs[study,]
    } else{
      ID = study
    }
    return(c(ID, Tx, sprintf("%.2f", rd), sprintf("%.2f", lev)))
  }
  
  df = as.data.frame(t(sapply(1:length(w), rd_line)))
  names(df)[(ncol(df)-2):ncol(df)] = c("Arm", "Residual Deviance", "Leverage")
  if(ncol(IDs) == 0){
    names(df)[1] = "Study"
  } else{
    names(df)[1:ncol(IDs)] = ID_nms
  }
  df
}





#' Posterior treatment ranking sample
#'
#' Turning the MCMC sample into ranking.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param events_are_bad Logical: are events bad or good?
#'
#' 
#' @return A matrix of posterior rankings. Each column represents a treatment. 
#'
#' @export
ranks = function(sim, events_are_bad){
  M = sim$meta.sim$sims.matrix
  D = cbind(0, M[, grep("d[[]", colnames(M))])
  if(!events_are_bad) D = -D
  
  apply(t(apply(D, 1, rank)), 2, as.integer)
}





#' SUCRA plots
#'
#' Plotting the Surface Under the Cumulative RAnking for each treatment.
#'
#' @param rank_tab A matrix of posterior rankings: the output of the \code{\link{ranks}} function.
#'
#' 
#' @return
#'
#' @export
SUCRA_plots = function(rank_tab){
  mat = apply(rank_tab[-nrow(rank_tab),-1], 2, function(x){gsub("%", "", x)})
  mat = apply(mat, 2, as.numeric)
  mat = apply(mat, 2, cumsum)
  mat[nrow(mat),] = 100
  mat = c(mat)
  df = as.data.frame(cbind(rep(1:(nrow(rank_tab)-1), nrow(rank_tab)-1), 
                           rep(colnames(rank_tab)[-1], each=nrow(rank_tab)-1), 
                           mat))
  names(df) = c("Rank", "Treatment", "Cummuative_Probability")
  df$Treatment = factor(df$Treatment, levels=colnames(rank_tab[,-1]), ordered=TRUE)
  df$Rank = as.integer(as.character(df$Rank))
  df$Cummuative_Probability = as.numeric(as.character(df$Cummuative_Probability))
  
  ggplot(df, aes(Rank, Cummuative_Probability)) +
    geom_line(color="red", size=2) +
    facet_wrap(~ Treatment, nrow=round(sqrt(nrow(rank_tab)))) + 
    scale_x_continuous(breaks= pretty_breaks()) +
    theme(strip.text = element_text(size=10),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold")) + 
    ylab("Cumulative Probability")
}





#' SUCRA table
#'
#' A table presenting the SUCRA values for all treatments in descending order.
#'
#' @param ranks_mat A matrix of posterior rankings: the output of the \code{\link{ranks}} function.
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#'
#' 
#' @return A table containing the different treatments and their respective SUCRA values, in a descending order.
#'
#' @export
SUCRA_table = function(ranks_mat, sim){
  single_sucra = function(i){
    t = table(ranks_mat[,i])
    miss = c(1:ncol(ranks_mat))[!(1:ncol(ranks_mat) %in% as.numeric(names(t)))]
    if(length(miss) > 0){
      t = c(rep(0, length(miss)), t)
      names(t)[1:length(miss)] = miss
      t = t[order(as.numeric(names(t)))]
    }
    sum(head(cumsum(t/sum(t)), -1))/(length(t)-1)
  }
  
  sucra = as.data.frame(cbind(1:length(sim$trt_names), sim$trt_names, 
                              paste0(sprintf("%.2f", sapply(1:ncol(ranks_mat), 
                                                            single_sucra)*100), "%"))
  )
  names(sucra) = c("Treatment Number", "Treatment Name", "SUCRA")
  sucra = sucra[order(as.numeric(gsub("%", "", sucra$SUCRA))/100, 
                      decreasing = TRUE),]
  
  sucra
}






#' Ranking posterior probabilities
#'
#' A table displaying the posterior ranking probabilities.
#'
#' @param ranks_mat A matrix of posterior rankings: the output of the \code{\link{ranks}} function.
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#'
#' 
#' @return A table containing the different treatments and their respective posterior ranking probabilities, as well as poterior median ranking.
#'
#' @export
rankogram = function(ranks, sim){
  single_rankogram = function(i){
    t = table(ranks[,i])
    miss = c(1:ncol(ranks))[!(1:ncol(ranks) %in% as.numeric(names(t)))]
    if(length(miss) > 0){
      t = c(rep(0, length(miss)), t)
      names(t)[1:length(miss)] = miss
      t = t[order(as.numeric(names(t)))]
    }
    
    sapply(t/sum(t), function(x){paste0(sprintf("%.1f", 100*x), "%")})
  }
  
  tab =  as.data.frame(rbind(cbind(1:ncol(ranks), 
                                   sapply(1:ncol(ranks), 
                                          single_rankogram)
  ),
  c("Median", apply(ranks, 2, median))
  )
  )
  names(tab) = c("Rank", sim$trt_names)
  tab
}





#' Plotting the rankogram
#'
#' A stacked bar plot in terrain colors illustrating the posterior ranking probabilities of each treatment.
#'
#' @param rank_tab a rankogram table, the output of the \code{\link{rankogram}} function.
#' 
#' 
#' @return 
#' 
#' @export
rankogram_plot = function(rank_tab){
  df = rank_tab[-nrow(rank_tab),-1]
  df = apply(df, 2, function(x){gsub("%", "", x)})
  df = cbind(rep(colnames(df), each=nrow(df)), 
             rep(1:nrow(df), nrow(df)),
             c(df))
  df = as.data.frame(df)
  names(df) = c("Treatment", "Rank", "Probability")
  df$Probability = as.numeric(as.character(df$Probability))
  df$Rank = factor(df$Rank, levels=1:max(unique(as.numeric(as.character(df$Rank)))), 
                   ordered=TRUE)
  df$Treatment = factor(df$Treatment, levels=colnames(rank_tab[,-1]), ordered=TRUE)
  m = nrow(rank_tab)
  pal = terrain.colors(m)
  
  ggplot() + geom_bar(aes(y = Probability, x = Treatment, fill = Rank), 
                      data = df,
                      stat="identity", 
                      position = position_fill(reverse = TRUE)) +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.text=element_text(size=ifelse(m <= 30, 12, 12*26/m)),
          axis.title=element_text(size=18, face="bold")) +
    scale_fill_manual(values=pal)
}






#' SUCRA bar plot
#'
#' A bar plot visualization of the SUCRA table.
#'
#' @param rank_tab a rankogram table, the output of the \code{\link{rankogram}} function.
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#'
#' 
#' @return 
#' 
#' @export
SUCRA_bars = function(rank_tab, sim){
  sucra = SUCRA_table(rank_tab, sim)
  sucra$SUCRA = as.numeric(gsub("%", "", sucra$SUCRA))/100
  sucra$`Treatment Name` = factor(sucra$`Treatment Name`, 
                                  levels=sucra$`Treatment Name`[nrow(sucra):1])
  
  ggplot(data=sucra, aes(x=sucra$`Treatment Name`, y=sucra$SUCRA)) +
    geom_bar(stat="identity", fill="blue", width=.25) +
    scale_y_continuous(limits = c(0, 1)) + 
    coord_flip() +
    xlab("Treatment") + ylab("SUCRA") +
    theme(axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=16, face="bold"),
          axis.text.x = element_text(hjust=1, size=12),
          axis.text.y = element_text(hjust=1, size=12),
          plot.title = element_text(size = 18, hjust=.5,
                                    face = "bold", 
                                    colour = "black")
    ) 
}







#' Model properties table
#'
#' A table detailing the different properties of the model.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#'
#' 
#' @return A data frame containing properties and their values.
#' 
#' @export
properties_table = function(sim){
  Property = c("Response", "Fixed/Random effect", "Heterogeneity Prior", "Meta-regressor",
               "Baseline adjusted", "DIC", "No. of MCMC iterations",
               "Burn in period", "Thinning factor", 
               "Number of chains", "Duration")
  str = strsplit(sim$response_type, "_")[[1]]
  sub_str = str[1]
  Value = paste0(toupper(substring(sub_str, 1, 1)),
                 substring(sub_str, 2), collapse='')
  if(sim$response_type == "normal_SMD") Value = "Normal (SMD)"
  
  Value = c(Value, ifelse(sim$RE, "Random", "Fixed"))
  Value = c(Value, sim$Hetero_Prior)
  Value = c(Value, ifelse(is.na(sim$Covar_Name), " ", sim$Covar_Name))
  Value = c(Value, ifelse(length(str) == 3, "Yes", "No"))
  
  Value = c(Value, round(sim$DIC, 1), sim$N_iter, 
            sim$burnin, sim$Thin, 
            sim$N_chains, trimws(sim$Runtime))
  
  Value = unname(Value)
  as.data.frame(cbind(Property, Value))
}






#' BUGS MCMC summary table
#'
#' Summary and basic diagnostics of the posterior sample
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#'
#' 
#' @return A data frame containing the model's basic parameters and their summary statistics, including the Gelman-Rubin diagnostic and effective sample size.
#' 
#' @export
Params_table = function(sim){
  nms = head(colnames(sim$meta.sim$sims.matrix),-1)
  nms = nms[-grep("dev[[]|D[[]|OR[[]|RR[[]|rhat[[]|theta[[]|T[[]", nms)]
  
  M = sim$meta.sim$summary
  M = M[row.names(M)%in%nms, ]
  M[,1:8] = apply(M[,1:8], 2, function(x){sprintf("%.2f", x)})
  row.names(M) = c()
  
  nms[grep("d[[]", nms)] = paste(sim$trt_names[-1], "trt. effect")
  nms[grep("tau", nms)] = "Precision"
  nms[nms == "B"] = "Regression coefficient"
  nms[nms == "Var"] = "Variance"
  
  M = as.data.frame(cbind(nms, M))
  names(M)[1] = "Parameter"
  M
}






#' MCID column 
#'
#' Calculating a single column in the MCID table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param baseline reference treatment number.
#' @param thresh The clinical threshold (minimal clinically important difference).
#' @param response_type One of "binary", "binary_base_risk", "normal", "normal_base_risk", "nomal_SMD", "poisson" and "poisson_base_risk".
#' @param events_are_bad Logical: are events bad or good?
#' 
#' @return A data frame containing the treatment names and their posterior probability of superiority with respect to the baseline treatment.
#'
#' @export
pairwise_superiority_column = function(sim, baseline, thresh, 
                                       response_type, events_are_bad){
  sims = sim$meta.sim$sims.matrix
  nms = sim$trt_names
  baseline_ind = which(nms == baseline)
  
  
  if(response_type %in% c("binomial", "poisson")){
    RR_inds = grep("RR[[]", colnames(sims))
    RR = sims[, RR_inds]
    inds_smaller = grep(paste(",", baseline_ind, "[]]", sep=''), colnames(RR))
    inds_greater = grep(paste("RR[[]", baseline_ind,",", sep=''), colnames(RR))
    
    RR[,inds_smaller] = 1/RR[,inds_smaller]
    RR = cbind(RR[,inds_smaller], 1, RR[,inds_greater])
    if(!events_are_bad) RR = 1/RR
    RRR = 1 - RR
    
    p_superior = unname(apply(as.matrix(RRR), 2, function(x){mean(x > thresh)}))
  } else{
    D_inds = grep("D[[]", colnames(sims))
    D = sims[, D_inds]
    inds_smaller = grep(paste(",", baseline_ind, "[]]", sep=''), colnames(D))
    inds_greater = grep(paste("D[[]", baseline_ind,",", sep=''), colnames(D))
    
    D[,inds_smaller] = -D[,inds_smaller]
    D = cbind(D[,inds_smaller], 0, D[,inds_greater])
    if(events_are_bad) D = -D
    
    p_superior = unname(apply(as.matrix(D), 2, function(x){mean(x > thresh)}))
  }
  
  df = data.frame(trt = nms, p_superior = p_superior)
}






#' MCID table 
#'
#' Calculating a single column in the MCID table.
#'
#' @param sim An NMA object: the output of the \code{\link{NMA}} function.
#' @param thresh The clinical threshold (minimal clinically important difference).
#' @param response_type One of "binary", "binary_base_risk", "normal", "normal_base_risk", "nomal_SMD", "poisson" and "poisson_base_risk".
#' @param events_are_bad Logical: are events bad or good?
#' 
#' @return An MCID matrix, where each entry is the posterior probability of superiority of the row treatment vs. the column teatment.
#'
#' @export
pairwise_superiority_table = function(sim, thresh, response_type, 
                                      events_are_bad){
  
  sapply(sim$trt_names, 
         function(x){pairwise_superiority_column(baseline=x,
                                                 sim=sim, thresh=thresh, 
                                                 response_type=response_type, 
                                                 events_are_bad=events_are_bad)$p_superior}
  )
}






#' MCID plot 
#'
#' A bar plot of MCID poetrior probabilities vs. a set baseline treatment.
#'
#' @param MCID_tab An MCID posterior probability matrix. The output of \code{\link{pairwise_superiority_table}}.
#' @param baseline reference treatment number.
#' @param thresh The clinical threshold (minimal clinically important difference).
#' 
#' @return 
#'
#' @export
MCID_Probs_Plot = function(MCID_tab, baseline, thresh, base_vs_all){
  if(base_vs_all){
    p_col = MCID_tab[colnames(MCID_tab) == baseline, ]
    plot_title = paste("Probability of Superiority,", baseline, "vs.")
  } else{
    p_col = MCID_tab[, colnames(MCID_tab) == baseline]
    plot_title = paste("Probability of Superiority vs.", baseline)
  }
  
  df = data.frame(trt = colnames(MCID_tab), p_superior = p_col) %>%
    filter(trt != baseline)
  
  df = df[order(df$p_superior, decreasing = TRUE),]
  df$trt = gsub("_", " ", df$trt)
  baseline = gsub("_", " ", baseline)
  df$trt = factor(df$trt, levels=df$trt[nrow(df):1])
  
  ggplot(data=df, aes(x=df$trt, y=df$p_superior)) +
    geom_bar(stat="identity", fill="blue", width=.25) +
    scale_y_continuous(limits = c(0, 1)) + 
    coord_flip() +
    ggtitle(
      plot_title,
      subtitle = paste("MCID =", thresh)
    ) +
    xlab("Treatment") + ylab("Posterior Probability") +
    theme(axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=16, face="bold"),
          axis.text.x = element_text(hjust=1, size=12),
          axis.text.y = element_text(hjust=1, size=12),
          plot.title = element_text(size = 18, hjust=.5,
                                    face = "bold", 
                                    colour = "black"),
          plot.subtitle = element_text(size = 18, hjust=.5)
    ) 
}





post_superiority_samples = function(my_NMA, baseline, events_are_bad, 
                                    base_vs_all){
  sims = my_NMA$meta.sim$sims.matrix
  nms = my_NMA$trt_names
  baseline_ind = which(nms == baseline)
  
  if(my_NMA$response_type %in% c("binomial", "poisson",
                                 "binomial_base_risk",
                                 "poisson_base_risk")){
    RR_inds = grep("RR[[]", colnames(sims))
    RR = sims[, RR_inds]
    inds_smaller = grep(paste(",", baseline_ind, "[]]", sep=''), colnames(RR))
    inds_greater = grep(paste("RR[[]", baseline_ind,",", sep=''), colnames(RR))
    
    RR[,inds_smaller] = 1/RR[,inds_smaller]
    RR = cbind(RR[,inds_smaller], 1, RR[,inds_greater])
    if(!events_are_bad) RR = 1/RR
    if(base_vs_all) RR = 1/RR
    RRR = 1 - RR
    
    return(RRR)
  } else{
    D_inds = grep("D[[]", colnames(sims))
    D = sims[, D_inds]
    inds_smaller = grep(paste(",", baseline_ind, "[]]", sep=''), colnames(D))
    inds_greater = grep(paste("D[[]", baseline_ind,",", sep=''), colnames(D))
    
    D[,inds_smaller] = -D[,inds_smaller]
    D = cbind(D[,inds_smaller], 0, D[,inds_greater])
    if(events_are_bad) D = -D
    if(base_vs_all) D = -D
    
    return(D)
  }
}




MCID_Density_Plots = function(my_NMA, M, thresh, baseline, base_vs_all){
  trt_nms = my_NMA$trt_names
  baseline_ind = which(trt_nms == baseline)
  M = M[,-baseline_ind]
  trt_nms_new = trt_nms[-baseline_ind]
  
  Better = (M > thresh)
  
  post_sup = paste0(round(apply(Better, 2, mean)*100, 1), '%')
  
  if(base_vs_all){
    text1 = paste0("vs. ", trt_nms_new)
  } else{
    text1 = paste0(trt_nms_new, " vs.")
  }
  text2 = paste0("Posterior Superiority: ", post_sup)
  xtitle = ifelse(my_NMA$response_type %in% c("binomial", "poisson"),
                  "RRR", "Mean Difference")
  
  plist = list()
  plist = vector("list", length(trt_nms_new))
  for(i in 1:length(trt_nms_new)){
    M_temp = M[,i]
    dat = with(density(M_temp), data.frame(x, y)) 
    
    p_temp = ggplot(data = dat, mapping = aes(x = x, y = y)) + 
      geom_area(mapping = aes(x = ifelse(x > thresh, x, 0)), 
                fill = "light blue") + 
      geom_area(mapping = aes(x = ifelse(x <= thresh, x, 0)), 
                fill = "red") + 
      ylim(c(0, max(dat$y))) +
      xlim(quantile(dat$x, c(.005, .995))) + 
      xlab(xtitle) + 
      ylab("Density") + 
      geom_line() + 
      labs(title = text1[i], subtitle = text2[i])
    
    plist[[i]] = p_temp
  }
  
  
  n = length(plist)
  nCol = floor(sqrt(n))
  suppressWarnings(do.call("grid.arrange", 
                           list(grobs = plist, ncol = nCol,
                                top = textGrob(trt_nms[baseline_ind],
                                               gp = gpar(fontsize = 20)))))
}
