library(foreach)
library(parallel)
library(doParallel)

#Create two-facet crossed simulated data set called sim_data
#persons x items x raters (in Brennan and Tong 2004, r = h)
#set sample sizes
np <- 100
ni <- 20
nr <- 2
#set true parameter values
sig_p <- 4
sig_i <- 2
sig_r <- 1
sig_pi <- 8
sig_pr <- sqrt(2)
sig_ir <- sqrt(3)
sig_pir <- 12
#############
#############

#Construct design matrix (only needs to be done once)
sim_data <- data.frame(matrix(0,np*ni*nr,ncol = 4))
colnames(sim_data) <- c("Person","Item","Rater","Score")
sim_data$Person <- rep(1:np,each = ni*nr)
sim_data$Item <- rep(1:ni,each = nr,times = np)
sim_data$Rater <- rep(1:nr,times = np*ni)

results_all <- NULL
registerDoParallel(cores = 4)

ptm <- proc.time()
r <- foreach(icount(50), .combine = rbind) %dopar% {
  #Create random effects (generate new for each replication)
  Zp <- rnorm(np,0,1)
  Zi <- rnorm(ni,0,1)
  Zr <- rnorm(nr,0,1)
  Zpi <- matrix(rnorm(np*ni,0,1),nrow = np,ncol = ni)
  Zpr <- matrix(rnorm(np*nr,0,1),nrow = np,ncol = nr)
  Zir <- matrix(rnorm(ni*nr,0,1),nrow = ni,ncol = nr)
  #Simulate Scores
  score_index <- 1
  for (p in 1:np) {
    for (i in 1:ni) {
      for (r in 1:nr) {
        Score_pir <- sig_p*Zp[p] + sig_i*Zi[i] + sig_r*Zr[r] +
          sig_pi*Zpi[p,i] + sig_pr*Zpr[p,r] + sig_ir*Zir[i,r] + sig_pir*rnorm(1)
        sim_data$Score[score_index] <- Score_pir
        score_index <- score_index + 1
      }
    }
  }
  result <- c(summaryCI(CalcGTheoryCI(sim_data,100,type = "o"),0.8,4)$Gstudy_Estimates)
  result
}
results_all <- rbind(results_all,r)

colMeans(results_all)
proc.time() - ptm
