# Author: Dennis van der Meer and Joakim Widén
# E-mail: dennis.van_der_meer[at]minesparis.psl.eu

# This script contains the objective function and other helper functions
# necessary to complete the run the optimization scripts. Importantly, it
# contains also the power flow scripts (in Matlab language) with which the
# candidate solutions are evaluated.

obj <- function(Q){
  V <- PF_opt(Q,n)                    # Where PF_opt is a function that calls Matlab to solve
  if (sum(V > 1.05 | V < 0.95) > 0) V_deviation <- Inf  
  else V_deviation <- sum(abs(1-V))                     
  V_deviation <- sum(abs(1-V)) # the power flow equations
  return(V_deviation)
}
localRepair=function(eda, gen, pop, popEval, f, lower, upper) {
  for(i in 1:nrow(pop)){ 
    x=pop[i,]
    x=ifelse(x<lower,lower,x) # assure x within
    x=ifelse(x>upper,upper,x) # bounds
    pop[i,]=x
    popEval[i]=f(x)
    if(sum(x==pop[i,])==length(lower)){
      popEval[i] <- popEval[i]
    } else {
      popEval[i]=f(x) # If a population member is changed, re-evaluate the function
    }
  }
  return(list(pop=pop,popEval=popEval))
}
edaSeedConstant <- function (eda, lower, upper) {
  popSize <- eda@parameters$popSize
  sapply(seq(along = lower),
         function (i) runif(popSize, lower[i], 0.5*upper[i]))
}

edaSelectWithinBounds <- function (eda, gen, pop, popEval, lower, upper) {
  truncFactor <- eda@parameters$truncFactor
  popSize <- eda@parameters$popSize

  if (is.null(truncFactor)) truncFactor <- 0.3

  tmp <- data.frame(which(pop < upper & pop > lower, arr.ind = T))
  x <- rep(1:nrow(pop),each=length(lower))
  tmp <- tmp[order(match(tmp$row, x)),]
  mat <- matrix(F, nrow = nrow(pop), ncol = 1) # Default assumption is FALSE
  # write.table(tmp, file = "tmp.txt", row.names = F, col.names = F)
  for(i in 1:nrow(pop)){mat[i,1] <- sum(tmp$row %in% i)==length(lower)} # Check if all DGs are represented
  # write.table(mat, file = "mat.txt", row.names = F, col.names = F)
  if(sum(mat==TRUE)>5){ # Check survival rate of the population with these conditions
    pop <- pop[mat,] # Disregard populations that do not represent all DGs
    popEval <- popEval[mat]
  } else if(sum(mat==TRUE)<=5){ # In case not enough of the population survives, sample new
    pop <- t(matrix(runif(popSize*sum(idx),lower,upper),ncol = popSize))
    popEval <- runif(popSize, 0, 0.05)
  }
  # And proceed as normal
  popOrder <- order(popEval)
  popOrder[seq(ceiling(truncFactor * length(popOrder)))]
}


PF_ret <- function(Q,n) {
  setVariable(matlab, iterations=0)
  setVariable(matlab, n=n)
  setVariable(matlab, Qopt=Q)
  evaluate(matlab, # I need to insert the optimized Qopt in Qgen for DG nodes only.
           "
           Qgen = zeros(size(idx'));
           Qgen(idx) = Qopt;
           ")
  evaluate(matlab,
           "
           % Slack bus voltage magnitude and angle
           Vmag1  = mySettings.slackNodeVoltage(n);
           delta1 = mySettings.slackNodePhaseAngle(n);
           
           % Initial values for voltage magnitude and angle
           Vmag  = ones(NUM_OF_BUSES, 1) * Vmag1;
           delta = ones(NUM_OF_BUSES, 1) * delta1;
           
           % Initialize number of iterations
           iterations = 0;
           while true
           
           % Update number of iterations
           iterations = iterations + 1;
           
           % Calculated real and reactive power at each node
           P_calc = zeros(NUM_OF_BUSES,1);
           Q_calc = zeros(NUM_OF_BUSES,1);
           for i = 2:NUM_OF_BUSES
           P_calc(i) = sum( Vmag(i) * Vmag .* Ymag(:,i) .* cos( theta(:,i) + delta - delta(i) ) );
           Q_calc(i) = -sum( Vmag(i) * Vmag .* Ymag(:,i) .* sin( theta(:,i) + delta - delta(i) ) );
           end
           
           % Scheduled real and reactive power at each node
           P_sched = Pgen(:,n) - Pload(:,n);
           % Q_sched = Qgen(:,n) - Qload(:,n);
           Q_sched = Qgen - Qload(:,n);        
           
           % Mismatch active and reactive power
           deltaP = P_sched(2:NUM_OF_BUSES) - P_calc(2:NUM_OF_BUSES);
           deltaQ = Q_sched(2:NUM_OF_BUSES) - Q_calc(2:NUM_OF_BUSES); 
           
           % Stop iteration if tolerance is reached or if maximum number of
           % iterations is reached
           if abs(sum(deltaP)) < mySettings.tolerance
           %exitFlags(n) = 1;
           break;
           end
           if iterations > mySettings.maxIterations
           %exitFlags(n) = -1;
           break;        
           end
           % Basic submatrix elements
           M = - (Vmag*ones(1,NUM_OF_BUSES)) .* (ones(NUM_OF_BUSES,1)*Vmag') .* Ymag .* sin(theta - delta*ones(1,NUM_OF_BUSES) + ones(NUM_OF_BUSES,1)*delta');
           for i = 2:NUM_OF_BUSES
           M(i,i) = -Q_calc(i) - Vmag(i)^2 * B(i,i);
           end
           N = - (Vmag*ones(1,NUM_OF_BUSES)) .* (ones(NUM_OF_BUSES,1)*Vmag') .* Ymag .* cos(theta - delta*ones(1,NUM_OF_BUSES) + ones(NUM_OF_BUSES,1)*delta');
           for i = 2:NUM_OF_BUSES
           N(i,i) = P_calc(i) - Vmag(i)^2 * G(i,i);
           end
           N2 = -N;
           for i = 2:NUM_OF_BUSES
           N2(i,i) = N(i,i) + 2*Vmag(i)^2 * G(i,i);
           end
           M2 = M;
           for i = 2:NUM_OF_BUSES
           M2(i,i) = -M(i,i) - 2*Vmag(i)^2 * B(i,i);
           end
           
           % Jacobian submatrices
           J11 = M(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J21 = N(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J12 = N2(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J22 = M2(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           
           % Jacobian
           J = [J11 J12; J21 J22];
           
           % Right-hand side
           b = [deltaP; deltaQ];
           
           % Turn ill-conditioned matrix warning temporarily off
           warning('off','MATLAB:illConditionedMatrix')
           warning('off','MATLAB:singularMatrix')
           
           % Solve equation system
           x = mldivide(J,b);
           
           % Turn warning on again
           warning('on','MATLAB:illConditionedMatrix')
           warning('on','MATLAB:singularMatrix')
           
           % Improved estimate changes
           d_delta = x(1:NUM_OF_BUSES-1);
           d_Vmag  = x(NUM_OF_BUSES:length(x));
           
           % Adjust estimates
           delta(2:NUM_OF_BUSES) = delta(2:NUM_OF_BUSES) + d_delta;
           Vmag(2:NUM_OF_BUSES)  = Vmag(2:NUM_OF_BUSES) .* (1 + d_Vmag);
           
           % Slack bus voltage
           delta(1) = delta1;
           Vmag(1)  = Vmag1;   
           
           end % This concludes the power flow analysis.
           
           % Store voltage values
           busVoltage(:,n) = complex(Vmag.*cos(delta), Vmag.*sin(delta) );
           currentV = complex(Vmag.*cos(delta), Vmag.*sin(delta) );
           ") # The last line means: only take the voltage at DG nodes
  V <- getVariable(matlab, "currentV"); V <- Re(V$currentV) # This is now a vector of real numbers
  V <- V/4160 
  return(V)
}


PF_opt <- function(Q,n) {
  setVariable(matlab, n=n)
  setVariable(matlab, Qopt=Q)
  evaluate(matlab, # Insert the optimized Qopt in Qgen for DG nodes only.
           "
           Qgen = zeros(size(idx'));
           Qgen(idx) = Qopt;
           ")
  evaluate(matlab,
           "
           % Slack bus voltage magnitude and angle
           Vmag1  = mySettings.slackNodeVoltage(n);
           delta1 = mySettings.slackNodePhaseAngle(n);
           
           % Initial values for voltage magnitude and angle
           Vmag  = ones(NUM_OF_BUSES, 1) * Vmag1;
           delta = ones(NUM_OF_BUSES, 1) * delta1;
           
           % Initialize number of iterations
           iterations = 0;
           while true
           
           % Update number of iterations
           iterations = iterations + 1;
           
           % Calculated real and reactive power at each node
           P_calc = zeros(NUM_OF_BUSES,1);
           Q_calc = zeros(NUM_OF_BUSES,1);
           for i = 2:NUM_OF_BUSES
           P_calc(i) = sum( Vmag(i) * Vmag .* Ymag(:,i) .* cos( theta(:,i) + delta - delta(i) ) );
           Q_calc(i) = -sum( Vmag(i) * Vmag .* Ymag(:,i) .* sin( theta(:,i) + delta - delta(i) ) );
           end
           
           % Scheduled real and reactive power at each node
           P_sched = Pgen(:,n) - Pload(:,n);
           % Q_sched = Qgen(:,n) - Qload(:,n);
           Q_sched = Qgen - Qload(:,n);        
           
           % Mismatch active and reactive power
           deltaP = P_sched(2:NUM_OF_BUSES) - P_calc(2:NUM_OF_BUSES);
           deltaQ = Q_sched(2:NUM_OF_BUSES) - Q_calc(2:NUM_OF_BUSES); 
           
           % Stop iteration if tolerance is reached or if maximum number of
           % iterations is reached
           if abs(sum(deltaP)) < mySettings.tolerance
           %exitFlags(n) = 1;
           break;
           end
           if iterations > mySettings.maxIterations
           %exitFlags(n) = -1;
           break;        
           end
           % Basic submatrix elements
           M = - (Vmag*ones(1,NUM_OF_BUSES)) .* (ones(NUM_OF_BUSES,1)*Vmag') .* Ymag .* sin(theta - delta*ones(1,NUM_OF_BUSES) + ones(NUM_OF_BUSES,1)*delta');
           for i = 2:NUM_OF_BUSES
           M(i,i) = -Q_calc(i) - Vmag(i)^2 * B(i,i);
           end
           N = - (Vmag*ones(1,NUM_OF_BUSES)) .* (ones(NUM_OF_BUSES,1)*Vmag') .* Ymag .* cos(theta - delta*ones(1,NUM_OF_BUSES) + ones(NUM_OF_BUSES,1)*delta');
           for i = 2:NUM_OF_BUSES
           N(i,i) = P_calc(i) - Vmag(i)^2 * G(i,i);
           end
           N2 = -N;
           for i = 2:NUM_OF_BUSES
           N2(i,i) = N(i,i) + 2*Vmag(i)^2 * G(i,i);
           end
           M2 = M;
           for i = 2:NUM_OF_BUSES
           M2(i,i) = -M(i,i) - 2*Vmag(i)^2 * B(i,i);
           end
           
           % Jacobian submatrices
           J11 = M(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J21 = N(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J12 = N2(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J22 = M2(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           
           % Jacobian
           J = [J11 J12; J21 J22];
           
           % Right-hand side
           b = [deltaP; deltaQ];
           
           % Turn ill-conditioned matrix warning temporarily off
           warning('off','MATLAB:illConditionedMatrix')
           warning('off','MATLAB:singularMatrix')
           
           % Solve equation system
           x = mldivide(J,b);
           
           % Turn warning on again
           warning('on','MATLAB:illConditionedMatrix')
           warning('on','MATLAB:singularMatrix')
           
           % Improved estimate changes
           d_delta = x(1:NUM_OF_BUSES-1);
           d_Vmag  = x(NUM_OF_BUSES:length(x));
           
           % Adjust estimates
           delta(2:NUM_OF_BUSES) = delta(2:NUM_OF_BUSES) + d_delta;
           Vmag(2:NUM_OF_BUSES)  = Vmag(2:NUM_OF_BUSES) .* (1 + d_Vmag);
           
           % Slack bus voltage
           delta(1) = delta1;
           Vmag(1)  = Vmag1;   
           
           end % This concludes the power flow analysis.
           
           % Store voltage values
           busVoltage(:,n) = complex(Vmag.*cos(delta), Vmag.*sin(delta) );
           currentV = complex(Vmag.*cos(delta), Vmag.*sin(delta) );
           currentV = currentV(idx,1);
           ") # The last line means: only take the voltage at DG nodes
  V <- getVariable(matlab, "currentV"); V <- Re(V$currentV) # This is now a vector of real numbers
  V <- V/4160 
  return(V)
}
# Run the following when n>1
myEdaRun <- function (eda, f, lower, upper, n, prevPop, prevPopEval) { # Define my own so I can get the copula
  gen <- 0
  terminate <- FALSE
  fEvals <- 0; fWrap <- function (...) { fEvals <<- fEvals + 1; f(...) }
  # bestEval <- NA # Commented this on 2020-05-13
  bestEval <- matrix(data = NA, nrow = eda@parameters$maxGen) # Added this on 2020-05-13 15 because 15 generations
  bestSol <- NA
  startTime <- proc.time()
  
  while (!terminate) {
    gen <- gen + 1
    
    if (gen == 1 && n == 1) {
      model <- NULL
      
      pop <- edaSeed(eda, lower, upper)
      popEval <- sapply(seq(length = nrow(pop)),
                        function (i) fWrap(pop[i, ]))
      
      r <- edaOptimize(eda, gen, pop, popEval, fWrap, lower, upper)
      pop <- r$pop
      popEval <- r$popEval
      # i <- which.min(popEval) # Added this on 2020-05-13
      # bestEval[gen,1] <- popEval[i]
    } else if(gen == 1 && n > 1){
      # Select the top 30%/2 from the entire previous population
      pop <- prevPop # Extract pop and popEval from the previous run
      popEval <- prevPopEval
      s <- edaSelect(eda, gen, pop, popEval)
      # s <- edaSelectWithinBounds(eda, gen, pop, popEval, lower, upper)
      selectedPop <- pop[s[1:(length(s)/2)],]
      selectedEval <- popEval[s[1:(length(s)/2)]]
      # Randomly seed a fresh population
      pop <- edaSeed(eda, lower, upper)
      popEval <- sapply(seq(length = nrow(pop)),
                        function (i) fWrap(pop[i, ]))
      # And rbind them with the previous population
      my.sample <- sample(nrow(pop),(length(s)/2))
      selectedPop <- rbind(selectedPop, pop[my.sample,])
      selectedEval <- rbind(selectedEval,popEval[my.sample])
      # # Check if each member of the population is within bounds:
      # for(i in 1:nrow(selectedPop)){ 
      #   x=selectedPop[i,]
      #   x=ifelse(x<lower,lower,x) # assure x within
      #   x=ifelse(x>upper,upper,x) # bounds
      #   selectedPop[i,]=x
      #   if(sum(x==selectedPop[i,])==length(lower)){
      #     selectedEval[i] <- selectedEval[i]
      #   } else {
      #     selectedEval[i]=f(x) # If I changed a population member, I should re-evaluate the function
      #   }
      # }
      # And continue as normal
      model <- edaLearn(eda, gen, model,
                        selectedPop, selectedEval, lower, upper)
      
      sampledPop <- edaSample(eda, gen, model, lower, upper)
      sampledEval <- sapply(seq(length = nrow(sampledPop)),
                            function (i) fWrap(sampledPop[i, ]))
      
      r <- edaOptimize(eda, gen, sampledPop, sampledEval,
                       fWrap, lower, upper)
      sampledPop <- r$pop
      sampledEval <- r$popEval
      
      r <- edaReplace(eda, gen, pop, popEval, sampledPop, sampledEval)
      pop <- r$pop
      popEval <- r$popEval
      # i <- which.min(popEval) # Added this on 2020-05-13
      # bestEval[gen,1] <- popEval[i]
    } else {
      s <- edaSelect(eda, gen, pop, popEval)
      # s <- edaSelectWithinBounds(eda, gen, pop, popEval, lower, upper)
      selectedPop <- pop[s, ]
      selectedEval <- popEval[s]
      # # Check if each member of the population is within bounds:
      # for(i in 1:nrow(selectedPop)){ 
      #   x=selectedPop[i,]
      #   x=ifelse(x<lower,lower,x) # assure x within
      #   x=ifelse(x>upper,upper,x) # bounds
      #   selectedPop[i,]=x
      #   if(sum(x==selectedPop[i,])==length(lower)){
      #     selectedEval[i] <- selectedEval[i]
      #   } else {
      #     selectedEval[i]=f(x) # If a population member is changed, should re-evaluate the function
      #   }
      # }
      
      model <- edaLearn(eda, gen, model,
                        selectedPop, selectedEval, lower, upper)
      
      sampledPop <- edaSample(eda, gen, model, lower, upper)
      sampledEval <- sapply(seq(length = nrow(sampledPop)),
                            function (i) fWrap(sampledPop[i, ]))
      
      r <- edaOptimize(eda, gen, sampledPop, sampledEval,
                       fWrap, lower, upper)
      sampledPop <- r$pop
      sampledEval <- r$popEval
      
      r <- edaReplace(eda, gen, pop, popEval, sampledPop, sampledEval)
      pop <- r$pop
      popEval <- r$popEval
    }
    
    edaReport(eda, gen, fEvals, model, pop, popEval)
    
    terminate <- edaTerminate(eda, gen, fEvals, pop, popEval)
    
    if (is.na(bestEval) || min(popEval) < bestEval) { 
      i <- which.min(popEval)
      # bestEval <- popEval[i]
      bestSol <- pop[i, ]
    }
    i <- which.min(popEval) # Added this on 2020-05-13
    bestEval[gen,1] <- popEval[i]
  }
  
  elapsedTime <- proc.time() - startTime
  result <- list(eda = eda,
                 f = f,
                 lower = lower,
                 upper = upper,
                 numGens = gen,
                 fEvals = fEvals,
                 bestEval = bestEval,
                 bestSol = bestSol,
                 pop = pop,
                 popEval = popEval,
                 cpuTime = sum(elapsedTime, na.rm = TRUE) - elapsedTime[3])
  
  result
}
# Run the following when n==1 because there's no previous population
initialEdaRun <- function (eda, f, lower, upper) {
  gen <- 0
  terminate <- FALSE
  fEvals <- 0; fWrap <- function (...) { fEvals <<- fEvals + 1; f(...) }
  # bestEval <- NA
  bestEval <- matrix(data = NA, nrow = eda@parameters$maxGen) # Added this on 2020-05-13
  bestSol <- NA
  startTime <- proc.time()
  
  while (!terminate) {
    gen <- gen + 1
    
    if (gen == 1) {
      model <- NULL
      
      pop <- edaSeed(eda, lower, upper)
      popEval <- sapply(seq(length = nrow(pop)),
                        function (i) fWrap(pop[i, ]))
      
      r <- edaOptimize(eda, gen, pop, popEval, fWrap, lower, upper)
      pop <- r$pop
      popEval <- r$popEval
    } else {
      s <- edaSelect(eda, gen, pop, popEval)
      # s <- edaSelectWithinBounds(eda, gen, pop, popEval, lower, upper)
      selectedPop <- pop[s, ]
      selectedEval <- popEval[s]
      # # Check if each member of the population is within bounds:
      # for(i in 1:nrow(selectedPop)){ 
      #   x=selectedPop[i,]
      #   x=ifelse(x<lower,lower,x) # assure x within
      #   x=ifelse(x>upper,upper,x) # bounds
      #   selectedPop[i,]=x
      #   if(sum(x==selectedPop[i,])==length(lower)){
      #     selectedEval[i] <- selectedEval[i]
      #   } else {
      #     selectedEval[i]=f(x) # If a population member is changed, re-evaluate the function
      #   }
      # }
      
      model <- edaLearn(eda, gen, model,
                        selectedPop, selectedEval, lower, upper)
      
      sampledPop <- edaSample(eda, gen, model, lower, upper)
      sampledEval <- sapply(seq(length = nrow(sampledPop)),
                            function (i) fWrap(sampledPop[i, ]))
      
      r <- edaOptimize(eda, gen, sampledPop, sampledEval,
                       fWrap, lower, upper)
      sampledPop <- r$pop
      sampledEval <- r$popEval
      
      r <- edaReplace(eda, gen, pop, popEval, sampledPop, sampledEval)
      pop <- r$pop
      popEval <- r$popEval
    }
    
    edaReport(eda, gen, fEvals, model, pop, popEval)
    
    terminate <- edaTerminate(eda, gen, fEvals, pop, popEval)
    
    if (is.na(bestEval) || min(popEval) < bestEval) { # Commented this on 2020-05-13
      i <- which.min(popEval)
      # bestEval <- popEval[i]
      bestSol <- pop[i, ]
    }
    i <- which.min(popEval) # Added this on 2020-05-13
    bestEval[gen,1] <- popEval[i]
  }
  
  elapsedTime <- proc.time() - startTime
  result <- list(eda = eda,
                 f = f,
                 lower = lower,
                 upper = upper,
                 numGens = gen,
                 fEvals = fEvals,
                 bestEval = bestEval,
                 bestSol = bestSol,
                 pop = pop,
                 popEval = popEval,
                 cpuTime = sum(elapsedTime, na.rm = TRUE) - elapsedTime[3])
  
  result
}

# UMDA function
umdaRun <- function (eda, f, lower, upper) {
  gen <- 0
  terminate <- FALSE
  fEvals <- 0; fWrap <- function (...) { fEvals <<- fEvals + 1; f(...) }
  # bestEval <- NA
  bestEval <- matrix(data = NA, nrow = eda@parameters$maxGen) # Added this on 2020-05-13
  bestSol <- NA
  startTime <- proc.time()
  
  while (!terminate) {
    gen <- gen + 1
    
    if (gen == 1) {
      model <- NULL
      
      pop <- edaSeed(eda, lower, upper)
      popEval <- sapply(seq(length = nrow(pop)),
                        function (i) fWrap(pop[i, ]))
      
      r <- edaOptimize(eda, gen, pop, popEval, fWrap, lower, upper)
      pop <- r$pop
      popEval <- r$popEval
    } else {
      s <- edaSelect(eda, gen, pop, popEval)
      selectedPop <- pop[s, ]
      selectedEval <- popEval[s]
      
      model <- edaLearn(eda, gen, model,
                        selectedPop, selectedEval, lower, upper)
      
      sampledPop <- edaSample(eda, gen, model, lower, upper)
      sampledEval <- sapply(seq(length = nrow(sampledPop)),
                            function (i) fWrap(sampledPop[i, ]))
      
      r <- edaOptimize(eda, gen, sampledPop, sampledEval,
                       fWrap, lower, upper)
      sampledPop <- r$pop
      sampledEval <- r$popEval
      
      r <- edaReplace(eda, gen, pop, popEval, sampledPop, sampledEval)
      pop <- r$pop
      popEval <- r$popEval
    }
    
    edaReport(eda, gen, fEvals, model, pop, popEval)
    
    terminate <- edaTerminate(eda, gen, fEvals, pop, popEval)
    
    if (is.na(bestEval) || min(popEval) < bestEval) { # Commented this on 2020-05-13
      i <- which.min(popEval)
      # bestEval <- popEval[i]
      bestSol <- pop[i, ]
    }
    i <- which.min(popEval) # Added this on 2020-05-13
    bestEval[gen,1] <- popEval[i]
  }
  
  elapsedTime <- proc.time() - startTime
  result <- list(eda = eda,
                 f = f,
                 lower = lower,
                 upper = upper,
                 numGens = gen,
                 fEvals = fEvals,
                 bestEval = bestEval,
                 bestSol = bestSol,
                 pop = pop,
                 popEval = popEval,
                 cpuTime = sum(elapsedTime, na.rm = TRUE) - elapsedTime[3])
  return(result)
}

NoControl <- function(){
  source("StartMatlabAndLoadData.R")
  
  evaluate(matlab,
           "
           %%%%%%%%%%%%%%%% Create the bus admittance matrix %%%%%%%%%%%%%%%%%%
           
           % Read line interconnection and impedance data
           NUM_OF_BUSES = myLineData.numOfBuses;
           NUM_OF_LINES = myLineData.numOfConnections;
           start_points = myLineData.connections(:,1);
           end_points   = myLineData.connections(:,2);
           resistance   = myLineData.impedances(:,1);
           reactance    = myLineData.impedances(:,2);
           
           % Line admittance matrix
           Yline = diag(1 ./ complex(resistance, reactance));
           
           % Bus incidence matrix
           Ibus = zeros(NUM_OF_BUSES, NUM_OF_LINES);
           for i = 1:NUM_OF_LINES
           start_ind = start_points(i);
           end_ind   = end_points(i);
           Ibus(start_ind, i) = 1;
           Ibus(end_ind, i) = -1;
           end
           
           % Bus admittance matrix
           Y = Ibus * Yline * Ibus';
           
           % Impedance matrix
           Z = zeros(size(Y));
           for i = 1:size(Y,1)
           for j = 1:size(Y,2)
           if Y(i,j) == 0
           Z(i,j) = 0;
           else
           Z(i,j) = -1/Y(i,j);
           end
           end
           end
           
           % Resistance and reactance matrices
           R = real(Z);
           X = imag(Z);
           
           % Real and imaginary parts of Y
           G = real(Y);
           B = imag(Y);
           
           % Magnitude and angle of Y
           Ymag  = abs(Y);
           theta = angle(Y);
           
           % Load and generation data
           Pload = myPowerData.activeLoad';
           Qload = myPowerData.reactiveLoad';
           Pgen  = myPowerData.activeGeneration';
           Qgen  = myPowerData.reactiveGeneration';
           % Number of time steps
           NUM_OF_TIME_STEPS = size(Pload, 2);
           % Initialize result matrices
           busVoltage = zeros(myLineData.numOfBuses, NUM_OF_TIME_STEPS);
           lineCurrent = zeros(myLineData.numOfConnections, NUM_OF_TIME_STEPS);
           lineLossesActive = zeros(myLineData.numOfConnections, NUM_OF_TIME_STEPS);
           lineLossesReactive = zeros(myLineData.numOfConnections, NUM_OF_TIME_STEPS);
           
           % Index the nodes with DGs
           idx = myPowerData.activeGeneration(144,:) ~= 0;
           lineData = myLineData;
           powerData = myPowerData;
           settings = mySettings;
           ")
  
  evaluate(matlab,
           "
           lineData = myLineData;
           powerData = myPowerData;
           settings = mySettings;
           %function [powerFlowResults, simulationTime, exitFlags, iterations] = runPowerFlow(lineData, powerData, settings)
           %% RUNPOWERFLOW
           %
           % Solves the power flow for the system defined by the input data structs. 
           % The power flow solution is found using the Newton-Raphson method.
           %
           % The INPUT data are:
           % 
           % lineData 
           % --------
           % .numOfBuses: Total number of buses in the grid.
           % .numOfConnections: Total number of lines in the grid.
           % .connections: Matrix [numOfConnections x 2] specifying the two inter-
           %  connected buses (1 to numOfBuses) for each line.
           % .impedances: Matrix [numOfConnections x 2] specifying the resistance and 
           %  reactance of each line [Ohm]. 
           % .maxCurrents: Array [numOfConnections x 1] specifying the maximum current 
           %  for each line [A].
           %
           % powerData
           % ---------
           % .numOfTimeSteps: Total number of time steps in the data.
           % .activeLoad: Matrix [numOfTimeSteps x numOfBuses] containing the active 
           %  power consumption data for each bus [W].
           % .reactiveLoad: As previous but for reactive power consumption [W].
           % .activeGeneration: As previous but for active power generation [W].    
           % .reactiveGeneration: As previous but for reactive power generation [W].
           %
           % settings
           % ---------
           % .slackNodeVoltage: Array [numOfTimeSteps x 1] with fixed slack node 
           %  voltage in each time step [V].
           % .slackNodePhaseAngle: Array [numOfTimeSteps x 1] with fixed slack node 
           %  phase angle in each time step [rad].
           % .tolerance: Size of the total power mismatches below which a solution is 
           %  considered to be found [W]. A reasonable value is 1.0000e-04.
           % .maxIterations: Maximum number of iterations, after which the solver is
           %  stopped with no solution found. A reasonable value is 1000 iterations.
           %
           % The OUTPUT data are:
           % 
           % powerFlowResults
           % ----------------
           % .busVoltage: Matrix [numOfBuses x numOfTimeSteps] containing the 
           %  resulting complex bus voltages [V].
           % .lineCurrent: Matrix [numOfConnections x numOfTimeSteps] containing the 
           %  resulting complex line currents [A].
           % .lineLossesActive: As previous but for active line losses [W].
           % .lineLossesReactive: As previous but for reactive line losses [W].
           %
           % simulationTime
           % --------------
           % Time for solving the total power flow over all time steps [sec].
           %
           % exitFlags
           % ---------
           % Array [1 x numOfTimeSteps] with exit flags for the power flow solution
           % method. 1 if the tolerance level is reached (solution found). -1 if the 
           % maximum number of iterations is reached (solution not found).
           %
           
           %% Create bus admittance matrix
           
           % Number of buses
           NUM_OF_BUSES = lineData.numOfBuses; %size(Y,1);
           
           
           %% Read load and generation data
           
           % Load and generation data
           Pload = powerData.activeLoad';
           Qload = powerData.reactiveLoad';
           Pgen  = powerData.activeGeneration';
           Qgen  = powerData.reactiveGeneration';
           
           
           %% Solve power flow in each time step
           
           % Start recording simulation time
           % tic
           
           % Number of time steps
           NUM_OF_TIME_STEPS = size(Pload, 2);
           
           % Initialize result matrices
           busVoltage = zeros(NUM_OF_TIME_STEPS, lineData.numOfBuses);
           lineCurrent = zeros(NUM_OF_TIME_STEPS, lineData.numOfConnections);
           lineLossesActive = zeros(NUM_OF_TIME_STEPS, lineData.numOfConnections);
           lineLossesReactive = zeros(NUM_OF_TIME_STEPS, lineData.numOfConnections);
           
           % Loop through all time steps
           for n = 1:NUM_OF_TIME_STEPS
           
           %disp(num2str(n))
           
           % Slack bus voltage magnitude and angle
           Vmag1  = settings.slackNodeVoltage(n);
           delta1 = settings.slackNodePhaseAngle(n);
           
           % Initial values for voltage magnitude and angle
           Vmag  = ones(NUM_OF_BUSES, 1) * Vmag1;
           delta = ones(NUM_OF_BUSES, 1) * delta1;
           
           % Initialize number of iterations
           iterations = 0;
           
           % Loop until mismatch tolerance level is reached
           while true
           
           % Update number of iterations
           iterations = iterations + 1;
           
           % Calculated real and reactive power at each node
           P_calc = zeros(NUM_OF_BUSES,1);
           Q_calc = zeros(NUM_OF_BUSES,1);
           for i = 2:NUM_OF_BUSES
           P_calc(i) = sum( Vmag(i) * Vmag .* Ymag(:,i) .* cos( theta(:,i) + delta - delta(i) ) );
           Q_calc(i) = -sum( Vmag(i) * Vmag .* Ymag(:,i) .* sin( theta(:,i) + delta - delta(i) ) );
           end
           
           % Scheduled real and reactive power at each node
           P_sched = Pgen(:,n) - Pload(:,n);
           Q_sched = Qgen(:,n) - Qload(:,n);        
           
           % Mismatch active and reactive power
           deltaP = P_sched(2:NUM_OF_BUSES) - P_calc(2:NUM_OF_BUSES);
           deltaQ = Q_sched(2:NUM_OF_BUSES) - Q_calc(2:NUM_OF_BUSES);
           
           % Stop iteration if tolerance is reached or if maximum number of
           % iterations is reached
           if abs(sum(deltaP)) < settings.tolerance
           exitFlags(n,1) = 1;
           break;
           end
           if iterations > settings.maxIterations
           exitFlags(n,1) = -1;
           break;        
           end
           
           % Basic submatrix elements
           M = - (Vmag*ones(1,NUM_OF_BUSES)) .* (ones(NUM_OF_BUSES,1)*Vmag') .* Ymag .* sin(theta - delta*ones(1,NUM_OF_BUSES) + ones(NUM_OF_BUSES,1)*delta');
           for i = 2:NUM_OF_BUSES
           M(i,i) = -Q_calc(i) - Vmag(i)^2 * B(i,i);
           end
           N = - (Vmag*ones(1,NUM_OF_BUSES)) .* (ones(NUM_OF_BUSES,1)*Vmag') .* Ymag .* cos(theta - delta*ones(1,NUM_OF_BUSES) + ones(NUM_OF_BUSES,1)*delta');
           for i = 2:NUM_OF_BUSES
           N(i,i) = P_calc(i) - Vmag(i)^2 * G(i,i);
           end
           N2 = -N;
           for i = 2:NUM_OF_BUSES
           N2(i,i) = N(i,i) + 2*Vmag(i)^2 * G(i,i);
           end
           M2 = M;
           for i = 2:NUM_OF_BUSES
           M2(i,i) = -M(i,i) - 2*Vmag(i)^2 * B(i,i);
           end
           
           % Jacobian submatrices
           J11 = M(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J21 = N(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J12 = N2(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           J22 = M2(2:NUM_OF_BUSES, 2:NUM_OF_BUSES);
           
           % Jacobian
           J = [J11 J12; J21 J22];
           
           % Right-hand side
           b = [deltaP; deltaQ];
           
           % Turn ill-conditioned matrix warning temporarily off
           warning('off','MATLAB:illConditionedMatrix')
           warning('off','MATLAB:singularMatrix')
           
           % Solve equation system
           x = mldivide(J,b);
           
           % Turn warning on again
           warning('on','MATLAB:illConditionedMatrix')
           warning('on','MATLAB:singularMatrix')
           
           % Improved estimate changes
           d_delta = x(1:NUM_OF_BUSES-1);
           d_Vmag  = x(NUM_OF_BUSES:length(x));
           
           % Adjust estimates
           delta(2:NUM_OF_BUSES) = delta(2:NUM_OF_BUSES) + d_delta;
           Vmag(2:NUM_OF_BUSES)  = Vmag(2:NUM_OF_BUSES) .* (1 + d_Vmag);
           
           % Slack bus voltage
           delta(1) = delta1;
           Vmag(1)  = Vmag1;   
           
           end
           
           % Store voltage values
           busVoltage(n,:) = complex(Vmag.*cos(delta), Vmag.*sin(delta) );
           
           % Calculate currents and losses in all lines
           for i = 1:size(lineData.connections,1)
           
           % Admittance and voltages for the connection
           bus1 = lineData.connections(i,1);
           bus2 = lineData.connections(i,2);
           Yij = abs(-Y(bus1,bus2));
           tij = angle(-Y(bus1,bus2));
           Vi = Vmag(bus1)/sqrt(3); Vj = Vmag(bus2)/sqrt(3); % Line-to-neutral voltages!
           di = delta(bus1); dj = delta(bus2);
           
           % Line current
           lineCurrent(n,i) = complex( Yij*Vi*cos(tij+di), Yij*Vi*sin(tij+di) ) - complex( Yij*Vj*cos(tij+dj), Yij*Vj*sin(tij+dj) );
           
           % Line losses
           lineLossesActive(n,i) = 3 * abs(lineCurrent(n,i)).^2 .* R(bus1,bus2);
           lineLossesReactive(n,i) = 3 * abs(lineCurrent(n,i)).^2 .* X(bus1,bus2);
           
           end
           
           end
           
           % Simulation time
           % simulationTime = toc;
           
           
           %% Set result struct
           
           %powerFlowResults.busVoltage = busVoltage;
           %powerFlowResults.lineCurrent = lineCurrent;
           %powerFlowResults.lineLossesActive = lineLossesActive;
           %powerFlowResults.lineLossesReactive = lineLossesReactive;
  
           myRes = abs(busVoltage) ./ slackNodeVoltage;
           ")
  
  # evaluate(matlab,'myRes')
  
  PF_res <- getVariable(matlab,'myRes'); PF_res <- PF_res$myRes
  PF_res <- data.frame(PF_res,time=seq(1,nrow(PF_res),1))
  
  return(PF_res)
}

