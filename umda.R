# Author: Dennis van der Meer
# E-mail: dennis.van_der_meer[at]minesparis.psl.eu

# This script runs the UMDA and saves the results in the working directory

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
umda <- CEDA(copula = "indep", margin = "kernel", popSize = 75,
              fEval = 0, fEvalTol = 1e-01, maxGen = 15)
popSize <- umda@parameters$popSize
umda@name <- "Univariate Marginal Distribution Algorithm"
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
  result <- umdaRun(umda, obj, lo_bound, up_bound)
  Voltage[n,] <- PF_ret(Q = result$bestSol, n = n) # Get the optimal voltages back
  Q_opt[n,] <- result$bestSol
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



