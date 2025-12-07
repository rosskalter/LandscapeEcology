# ESPM 137 Basic Population Genetics Functions

drift <- function(pops, N, gens, plot=TRUE){
  # Initialize populations
  pop <- list()
  for(i in 1:pops){
    pop[[i]] <- rbinom(2*N, 1, 0.5)
  }
  
  # Run through generations
  a <- matrix(nrow=gens, ncol=pops)
  for(j in 1:gens){
    for(k in 1:pops){
      pop[[k]] <- sample(pop[[k]], 2*N, replace=TRUE)
      a[j, k] <- mean(pop[[k]])
    }
  }
  if(plot==TRUE){
    colors <- brewer.pal(n=pops+1, name = "Set1")
    plot(1:gens, a[,1], type="l", col=colors[1], ylim=c(0,1), xlab="Time (Generations)", ylab="Allele Frequency")
    for(z in 2:pops){
      lines(1:gens, a[,z], type="l", col=colors[z])
    }
  }
  vars <- lapply(pop, var)  
  return(vars)
}


dispersal_drift <- function(mig, N, gens, plot=TRUE){
  # Initialize populations
  pop <- list()
  newpop <- list()
  pop[[1]] <- rbinom(2*N, 1, 0.5)
  pop[[2]] <- rbinom(2*N, 1, 0.5)
  
  # Run through generations
  a <- matrix(nrow=gens, ncol=2)
  for(i in 1:gens){
    newpop[[1]] <- c(sample(pop[[1]], 2*N, replace=TRUE), sample(pop[[2]], 2*N*mig, replace=TRUE))
    newpop[[2]] <- c(sample(pop[[2]], 2*N, replace=TRUE), sample(pop[[1]], 2*N*mig, replace=TRUE))   
    pop[[1]] <- sample(newpop[[1]], 2*N, replace=FALSE)
    pop[[2]] <- sample(newpop[[2]], 2*N, replace=FALSE)
    a[i,] <- c(mean(pop[[1]]), mean(pop[[2]]))
  }
  if(plot==TRUE){
    plot(1:gens, a[,1], type="l", col="darkred", ylim=c(0,1), xlab="Time (Generations)", ylab="Allele Frequency")
    lines(1:gens, a[,2], type="l", col="darkblue")
  }
  vars <- lapply(pop, var)
  return(vars)
}


# Sum coefficients for each predictor (each has 3 splines)
coeffs <- function(gdm.model){
  coefSums <- c()
  for (i in 1:length(gdm.model$predictors)){
    j <- (i * 3) - 2
    coefSums[i] <- sum(gdm.model$coefficients[j:(j+2)])
  }
  
  # Add those values to a simple data frame
  coeffs <- data.frame(predictor = gdm.model$predictors, coefficient = coefSums)
  return(coeffs)
}
