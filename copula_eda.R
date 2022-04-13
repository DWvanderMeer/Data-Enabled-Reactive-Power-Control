# Author: Dennis van der Meer
# E-mail: dennis.van_der_meer[at]minesparis.psl.eu

# This script runs the copula EDA and saves the results in the working directory

library(copulaedas)
library(R.matlab)

WORKING_DIR <- "~/Google Drive/My Drive/PhD-Thesis/My-papers/CopulaEDAs/Data-Enabled-Reactive-Power-Control/"
source(file.path(WORKING_DIR,"functions.R"))
####################################### Construct a copula EDA ####################################### 
source(file.path(WORKING_DIR,"functions.R"))
setMethod("edaSeed", "EDA", edaSeedConstant)
setMethod("edaTerminate", "EDA", edaTerminateCombined(edaTerminateEval, edaTerminateMaxGen))
setMethod("edaReport", "EDA", edaReportDumpPop)
setMethod("edaOptimize","EDA",localRepair)
gceda <- CEDA(copula = "gaussian", margin = "kernel", popSize = 75,
              fEval = 0, fEvalTol = 1e-01, maxGen = 15)
popSize <- gceda@parameters$popSize
gceda@name <- "Gaussian Copula Estimation of Distribution Algorithm"
##################################### Load data and start Matlab server ##################################### 
source(file.path(WORKING_DIR,"StartMatlabAndLoadData.R"))
idx <- getVariable(matlab, "idx"); idx <- idx$idx
################################################ Data allocation ############################################
N = 288 # Number of time steps to evaluate
Voltage <- matrix(data = NA, nrow = N, ncol = length(idx)) # length(idx) is the number of nodes of the LV grid
Q_opt <- matrix(data = NA, nrow = N, ncol = sum(idx)) # sum(idx) is the number of PV nodes of the LV grid
lower_bounds <-  matrix(data = NA, nrow = N, ncol = sum(idx))
upper_bounds <-  matrix(data = NA, nrow = N, ncol = sum(idx))
timed <- matrix(NA, nrow = N, ncol = 5) # To store the run times in
results <- list()
###################################### Run the optimization over time #######################################
ptm <- proc.time()
for(n in 1:N){ # Loop over the time steps
  tic <- proc.time() # Start timer
  newdir <- paste0("Populations/N",n)
  dir.create(newdir)      # should test for error
  cwd <- getwd()          # CURRENT dir
  setwd(newdir)
  
  setVariable(matlab, n=n)
  ############# Set lower and upper bounds #############
  up_bound <- getVariable(matlab,'maxReactiveGeneration')
  up_bound <- up_bound$maxReactiveGeneration[n,as.logical(idx)]
  lo_bound <- -up_bound
  PV_power <- getVariable(matlab,'activeGeneration')
  PV_power <- PV_power$activeGeneration[n,as.logical(idx)]
  if(sum(PV_power)==0){
    lo_bound <- rep(0,sum(idx))
    setMethod("edaSeed", "EDA", edaSeedConstant)
  } else if(sum(PV_power)!=0){
    setMethod("edaSeed", "EDA", edaSeedUniform)
  }
  ####### Solve the optimization ####### 
  if(n==1){
    result <- initialEdaRun(gceda, obj, lo_bound, up_bound)
    Voltage[n,] <- PF_ret(Q = result$bestSol, n = n) # Get the optimal voltages back
    Q_opt[n,] <- result$bestSol
  } else if(n>1){
    # The following bit is to select the part of the population that falls within bounds
    prevPop <- result$pop
    prevPopEval <- result$popEval
    tmp <- data.frame(which(prevPop < up_bound & prevPop > lo_bound, arr.ind = T))
    x <- rep(1:nrow(prevPop),each=sum(idx))
    tmp <- tmp[order(match(tmp$row, x)),]
    mat <- matrix(F, nrow = nrow(prevPop), ncol = 1)
    for(i in 1:nrow(prevPop)){mat[i,1] <- sum(tmp$row %in% i)==sum(idx)} # Check if all DGs are represented
    if(sum(mat==TRUE)>5){ # Check survival rate of the population with these conditions
      prevPop <- prevPop[mat,] # Disregard populations that do not represent all DGs
      prevPopEval <- prevPopEval[mat]
    } else if(sum(mat==TRUE)<=5){ # In case not enough of the population survives, sample new
      prevPop <- t(matrix(runif(popSize*sum(idx),lo_bound,up_bound),ncol = popSize))
      prevPopEval <- runif(popSize, 0.025, 0.075)
    }
    result <- myEdaRun(gceda, obj, lo_bound, up_bound, n, prevPop, prevPopEval)
    Voltage[n,] <- PF_ret(Q = result$bestSol, n = n) # Get the optimal voltages back
    Q_opt[n,] <- result$bestSol
  }
  # # The standard (without collecting 50% of the samples):
  # result <- edaRun(gceda, obj, lo_bound, up_bound) 
  # Voltage[n,] <- PF_ret(Q = result@bestSol, n = n) # Get the optimal voltages back
  # Q_opt[n,] <- result@bestSol
  lower_bounds[n,] <- lo_bound
  upper_bounds[n,] <- up_bound
  results[[n]] <- result
  setwd(cwd)
  timed[n,] <- proc.time() - tic # toc
}
proc.time() - ptm
close(matlab)

# Save the results
write.table(Voltage, file = "busVoltage.txt", row.names = F, col.names = F)
write.table(Q_opt, file = "optimalQ.txt", row.names = F, col.names = F)
save(results, file = "Results.RData")
write.table(timed, file = "runTime.txt",col.names = F, row.names = F)
# # Plot the voltages
# matplot(Voltage,type = "l",ylim = c(0.98,1.02))
# # Check to see whether there are no more constraint violations (TRUE if this is so)
# sum(Q_opt>lower_bounds)+sum(Q_opt==lower_bounds)==nrow(Q_opt)*ncol(Q_opt)

