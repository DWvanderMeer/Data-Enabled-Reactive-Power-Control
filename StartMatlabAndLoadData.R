# Author: Dennis van der Meer and Joakim Widén
# E-mail: dennis.van_der_meer[at]minesparis.psl.eu

# This script starts the Matlab server and loads the required data. 

library(R.matlab)
setwd(file.path(WORKING_DIR))
####################################### Start the Matlab server #########################################
options(matlab="/Applications/MATLAB_R2021b.app/bin/matlab") 
Matlab$startServer()
print(system.file("externals", package = "R.matlab"))
Matlab$startServer(port = 9999)
matlab <- Matlab()
print(matlab)
isOpen <- open(matlab)
if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")
print(matlab)
#########################################################################################################

####################### Load the necessary data and create bus admittance matrix ########################
evaluate(matlab,
         "warning('off','MATLAB:table:ModifiedAndSavedVarnames') % Suppress warning
         %%%%%%%%%%%%%%%%%%%%%%%%%%%% Line data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         % The configurations are almost always the same, just in different orders.
         % They are listed first:
         a = [0.4576 1.0780];
         b = [0.4615 1.0651];
         c = [0.4666 1.0482];
         d = [0 0];
         e = [1.3292 1.3475];
         
         % The single line equivalent of the three phase system then becomes:
         Config1 = (1/3) * sum([a; c; b],1);
         Config2 = (1/3) * sum([c; b; a],1);
         Config3 = (1/3) * sum([b; a; c],1);
         
         Config4 = (1/3) * sum([b; c; a],1);
         Config5 = (1/3) * sum([c; a; b],1);
         Config6 = (1/3) * sum([a; b; c],1);
         
         Config7 = (1/3) * sum([a; d; b],1);
         Config8 = (1/3) * sum([a; b; d],1);
         Config9 = (1/3) * sum([e; d; d],1);
         
         Config10 = (1/3) * sum([d; e; d],1);
         Config11 = (1/3) * sum([d; d; e],1);
         Config12 = (1/3) * sum([1.5209 0.7521; 1.5329 0.7162; 1.5209 0.7521],1);
         
         Configurations = {Config1,Config2,Config3,Config4,Config5,Config6,Config7,Config8,Config9,Config10,Config11,Config12,Config12};
         
         % Read the line data from the .xls file:
         T = readtable('myLineData.xls', 'Range', 'A3:D121');
         connections = [T.NodeA T.NodeB];
         impedances = (1/5280)*(T.Length_ft__ .* reshape(cell2mat(Configurations(T.Config_)),[2,118])');
         
         numOfConnections = size(connections,1);
         numOfBuses = max(max(connections));
         maxCurrents = repmat(275,numOfConnections,1);
         
         myLineData.numOfBuses = numOfBuses;
         myLineData.numOfConnections = numOfConnections;
         myLineData.connections = connections;
         myLineData.impedances = impedances;
         myLineData.maxCurrents = maxCurrents;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%% Power data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Instead of one time step, I want to run a day with the load and 
         % generation profiles given. The nodes from below still apply.
         PV_CS  = csvread('cleardayprofile.csv'); PV_CS = 0.9*PV_CS;
         Load_P = csvread('peakloadnormalized.csv');
         Load_M = csvread('minloadnormalized.csv');

         % Average PV to 15 minute resolution:
         % x = (1:length(PV_CS))';
         % xi = (15:15:1440)';
         % PV_CS = interp1q(x,PV_CS,xi);
         % PV_OC = interp1q(x,PV_OC,xi);
         % Average PV to 5 minute resolution:
         x = (1:length(PV_CS))';
         xi = (5:5:1440)';
         PV_CS = interp1q(x,PV_CS,xi);
         load('GHIM.mat')
         xx = (1:length(GHIM))';
         xixi = (20:20:size(GHIM,1))';
         GHIM_int = zeros(size(PV_CS,1),size(GHIM,2));
         for i=1:size(GHIM,2)
             GHIM_int(:,i) = 0.85*normalize(interp1q(xx,GHIM(:,i),xixi),'range');
         end
         
         % Interpolate Load to 5 minute resolution:
         x = (1:length(Load_P))';
         xi = (1:.33:96)';
         Load_P = interp1q(x,Load_P,xi);
         Load_M = interp1q(x,Load_M,xi);

         AL = readtable('spotloadsdata.xls', 'Range', 'A4:D89');
         RL = readtable('spotloadsdata.xls', 'Range', 'F4:I89');
         CG = readtable('capdata.xls', 'Range', 'A4:D8');
         numOfTimeSteps = length(PV_CS);
         
         activeLoad = zeros(numOfTimeSteps,numOfBuses);
         for i=1:numOfTimeSteps
           for j=1:height(AL)
            activeLoad(i,AL.Node(j)) = 1000*mean(AL{j,2:4},2) * Load_P(i,1);
           end
         end
         
         reactiveLoad = zeros(numOfTimeSteps,numOfBuses);
         for i=1:numOfTimeSteps
           for j=1:height(RL)
            reactiveLoad(i,RL.Node(j)) = 1000*mean(RL{j,2:4},2) * Load_P(i,1);
           end
         end
         
         % PV power generation:
         activeGeneration = zeros(numOfTimeSteps,numOfBuses);
         lstOfPVnodes = [19 49 57 67 80 84 96 116 118 119];
         lstOfPVsizes = 8405*[70 140 200 50 90 160 200 140 60 80]; % in W, leads to 200% PV penetration
         
         % Clear sky:
         %for i=1:length(lstOfPVnodes)
         %   activeGeneration(:,lstOfPVnodes(i)) = PV_CS*lstOfPVsizes(i);
         %end
         
         % Mixed sky:
         for i=1:length(lstOfPVnodes)
            activeGeneration(:,lstOfPVnodes(i)) = GHIM_int(:,i)*lstOfPVsizes(i);
         end
         
         reactiveGeneration = zeros(numOfTimeSteps,numOfBuses);
         for j=1:height(CG)
          reactiveGeneration(:,CG.Node(j)) = 1000*mean(CG{j,2:4},2);
         end
         
         maxReactiveGeneration = zeros(numOfTimeSteps,numOfBuses);
         maxReactiveGeneration(:,lstOfPVnodes) = sqrt(lstOfPVsizes.^2 - activeGeneration(:,lstOfPVnodes).^2);
         
         myPowerData.numOfTimeSteps = numOfTimeSteps;
         myPowerData.activeLoad = activeLoad;
         myPowerData.reactiveLoad = reactiveLoad;
         myPowerData.activeGeneration = activeGeneration;
         myPowerData.reactiveGeneration = reactiveGeneration;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         slackNodeVoltage = repmat(4160,length(activeGeneration),1); % Volt
         slackNodePhaseAngle = zeros(length(activeGeneration),1);
         tolerance = 1.000e-04;
         maxIterations = 1000;
         
         mySettings.slackNodeVoltage = slackNodeVoltage;
         mySettings.slackNodePhaseAngle = slackNodePhaseAngle;
         mySettings.tolerance = tolerance;
         mySettings.maxIterations = maxIterations;

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
         ")

