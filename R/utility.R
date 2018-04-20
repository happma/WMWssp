####################################################################################################################################
### Filename:    utility.R
### Description: Function 'sim_power' for performing a power simulation and 'insert_row' to add a row in a data.frame
###
###
###
###
####################################################################################################################################


#' @keywords internal
sim_power <- function(x1,x2,nsim,n1,n2) {

  simpower <- 0
  cat("Simulation:\n")

  for(i in 1:nsim){
    if(i%%1000==0){
      cat(paste(i, "/", nsim, "\n"))
    }
    z1 <- sample(x1, size = ceiling(n1), prob = NULL, replace = TRUE)
    z2 <- sample(x2, size = ceiling(n2), prob = NULL, replace = TRUE)
    df = data.frame(grp = as.factor(c(rep(1,ceiling(n1)), rep(2,ceiling(n2)))), z = c(z1,z2) )

    if(rank.two.samples(z~grp, data = df, wilcoxon = "asymptotic", info = FALSE, shift.int=FALSE, alternative = "two.sided")$Wilcoxon$p.Value <= 0.05){
      simpower <- simpower + 1
    }
  }
  cat("\n")
  return (simpower)
}


#' @keywords internal
insert_row <- function(df, newrow, r) {
  df[seq(r+1,nrow(df)+1),] <- df[seq(r,nrow(df)),]
  df[r,] <- newrow
  return(df)
}
