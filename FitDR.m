%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit dosage-response curves to the data
% Dosage-response curves are solutions of the cubic equation that describes growth inhibition 
% upon exposure to ribosome-targeting drugs for a given organism as described in 
% Mol Syst Biol 2015 Mar;11(3):796. doi: 10.15252/msb.20145949

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS (data is available from in Supplementary Table 1):  

% (1) drug values in mg/mL 
% (2) growth rate scaled in units of the wild-type growth
%      data is fit to the mean over 3 replicates
% (3) standard deviation of the growth rate across replicates; one value per drug concentration 

% OUTPUTS 

% Dose-response parameters: 
% (1) rate scale lambda*
% (2) concentration scale d* 
% (3) drug-free growth rate lambda0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generics:

% Vectors starting with 's3P..' are used to store dose-response parameters. 
% Entries: (lambda*/2, 1/d*, 1/lambda0)

% ErrorInferredPars stores the error score for models with different
% assumptions on the parameters. It contains 3 entries corresponding to:
% (1) Model with 3 free parameters
% (2) Model with 2 free parameters; lambda* is fixed
% (3) Model with 2 free parameters; lambda*d* is fixed

% Constrained models are used to generate initial conditions for the
% unconstrained model (3 parameters model). Constrained models can also be
% used for model comparison.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all

% Module 1 ----------------------------------------------------------------------

% Input here the data of your favorite organism (available from Supplementary Table 1)
drugData = [];
growthData = [];
stdGrowthAux = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: dosage-response curve of the wild-type - uncomment to run
% drugData = [0     1     2     4     8    10]; % Drug values
% growthData = [1.00    0.93    0.95    0.81    0.62   0 ];  % Growth values
% stdGrowthAux = [0.03    0.05    0.06    0.04    0.03    0.03]; % Standard deviation of growth measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Minimum growth rate and detection threshold
minStd = 0.03;
stdGrowthData = max(stdGrowthAux, minStd);
weights = 1./stdGrowthData.^2; % fitting weights


% Simplify notation to send it into the fitting procedure
x = drugData;
y = growthData;

% Estimate initial condition for the parameter d*
ic50min = x(find(y > y(1)/2, 1, 'last'));
ic50max = x(find(y < 0.1, 1, 'first'));

if isempty(ic50min)
    ic50min = 0;
end

if isempty(ic50max)
    ic50max = x(end);
end


%%% Optimization settings nonlinear least squares fit..   
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt', 'Display', 'off');
options.MaxFunctionEvaluations = 3000;
options.MaxIterations = 1000;
options.TolX = 1.000000e-25;
options.TolFun = 1.000000e-25;
options.MaxFunEvals = 500*3;
OptimalityTolerance =  1.000000e-25;


% Define vectors to store parameters
s3PConstLsIC50 = zeros(3, 1);
s3PConstLs = zeros(3, 1);
s3P = zeros(3, 1);
errorInferredPars = zeros(3, 1);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIT MODEL WITH CONSTANT lambda*d* = 1/const
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The value of lambdaStarHalf is defined before the loop and is used as initial condition for different fits. We use
% the value inferred for the wild-type with the 3-parameters model

% The value of const in lambda*d* = 1/const used in the constrained model gives a low cummulative error score 
% accross the set of mutants. This model is only used to generate initial
% conditions to be fed into the fitting routine of the other models. Fitting alternative models is not necessary but 
% gives more robust initial conditions and good check points for the analysis

lambdaStarHalf = 0.0918; % wild-type Lambda*/2
const = 0.2728;  % constant for constrained model with lambda*d* = 1/const
detectionThreshold = 0.03; % growth detection threshold
minLDGrowthSTD =  0.03; % minimum growth std


% Fit model constrained model with non-linear least squares routine; 
minFconstLSicR = @(pars) sseval2ParsLsIC50(pars, const, x, y, weights, detectionThreshold);
[bestxNLinCr, resnormNLinCr] =lsqnonlin(minFconstLSicR, [lambdaStarHalf, 1/(0.5*(ic50min + ic50max)), 1/y(1)], [], [], options);
bestx = bestxNLinCr;

bestx(2) = const*(2*bestx(1)); 
bestxError = resnormNLinCr;

% Store parameters of constrained model
s3PConstLsIC50(1) = bestx(1);
s3PConstLsIC50(2) = bestx(2);
s3PConstLsIC50(3) = bestx(3);    
errorInferredPars(3) = bestxError;

%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIT MODEL WITH CONSTANT lambda* 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fun = @(pars) sseval2ParsLs(pars, lambdaStarHalf, x, y, weights, detectionThreshold);
optionsFminsearch = optimset('TolX', 10^-25, 'TolFun', 10^-25, 'MaxIter', 500*2, 'DiffMinChange', 10^-10, 'MaxFunEvals',  500*2, 'Display', 'off'); %, 'Display','iter');

bestx(1) =   lambdaStarHalf;
bestx(2:3) = fminsearch(fun, [s3PConstLsIC50(2), 1/y(1)], optionsFminsearch);
bestxError =  sseval3Pars(bestx, x, y, weights, detectionThreshold);


% Store parameters of constrained model
s3PConstLs(1) = bestx(1);
s3PConstLs(2) = bestx(2);
s3PConstLs(3) = bestx(3);    
errorInferredPars(2) = bestxError;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIT 3-PARAMETERS MODEL TO DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial condition
initialStep = s3PConstLsIC50;

fun = @(pars) sseval3Pars(pars, x, y, weights, detectionThreshold);
optionsFminsearch = optimset('TolX', 10^-25, 'TolFun', 10^-25, 'MaxIter', 500*3, 'DiffMinChange', 10^-10, 'MaxFunEvals',  500*3,  'Display', 'off');

bestx = fminsearch(fun, initialStep, optionsFminsearch);
bestxError =  sseval3Pars(bestx, x, y, weights, detectionThreshold);


% Store parameters of the unconstrained model
s3P(1) = bestx(1);
s3P(2) = bestx(2);
s3P(3) = bestx(3);    
errorInferredPars(1) = bestxError;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots dosage-response curve of a given organism
% Data + dosage-response curves

% Define drug range for plotting..

dDrug = 0.02; % Size of step
dMax = 64; % Maximum growth data. Should be set as required by tested drug range
d = 0:dDrug:dMax; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values used in the paper; here to illustrate plotting procedure

dScale = 8.66; % d50wt; use dScale = 1 for plots in mug/mL; adjust figure label accordingly
gScale = 0.982; % Lambda0 of the wild-type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Define colors for data and model dosage-response curve
dataColor = [0 0 0];
modelColor = [0.8 0.2 0.2];


figure
hold on
errorbar(drugData/dScale, growthData/gScale, stdGrowthData, '.', 'Markersize', 10, 'Color', dataColor, 'Linewidth', 1.2)
plot(d/dScale, makeG(s3P, d), 'Linewidth', 1.2, 'Color', modelColor)

ylim([0 1.05])
xlabel('d/d50_{wt}')
ylabel('scaled growth, \lambda')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Module 2 ----------------------------------------------------------------------
% Ensemble of synthetic growth inhibition curves for error analysis

nReps = 1000; % Set number of replicates as required. Paper analysis uses 1000 fits tht were tested for convergence (exitFlag = 1)

if nReps

    stdForError = stdGrowthData; % Std will be use to define distribution to sample synthetic dosage-response curves
    nDataPoints = length(growthData);
    initialGuess = s3P; % Use optimal parameters as initial condition
    nMeasurements = length(growthData); 
    
    % Generate nReps synthetic dosage-response curves
    for n = 1:nReps
       
        clear syntheticDR
        for i = 1:nMeasurements
            if growthData(i) == 0
               % If no growth is reported, sample from a uniform distribution between (0, detectionThreshold) 
               syntheticDR(i) = detectionThreshold*rand;
            else
                
               % If there is growth, sample from a Gaussian with mean and standard deviation set
               % by three replicates
               syntheticDR(i) = normrnd(growthData(i), stdForError(i)); 
            end
        end
      
        % Adjust synthetic data if we have two subsequent points with no
        % growth -- if needed, i.e., sampled growth rate at the largest drug concentration is
        % higher
        
       if sum(growthData(end-1:end) == 0 & syntheticDR(end-1) < syntheticDR(end)) 
     
           auxGrowth = syntheticDR(end);
           syntheticDR(end) = syntheticDR(end-1);
           syntheticDR(end-1) = auxGrowth;

  
       end
   
        saveGrowth(n, :) = max(syntheticDR, 0);

    end

    % Use parfor to speed up the fits; define vectors to save parametes, error scores and exit flags
    auxParForA = zeros(nReps, 1); auxParForB = zeros(nReps, 1); auxParForC = zeros(nReps, 1); auxParForE = zeros(nReps, 1); auxForExitFlag = zeros(nReps, 1);
   
    parfor n = 1:nReps

        fun = @(pars) sseval3Pars(pars, x, saveGrowth(n, :), weights, detectionThreshold);
        optionsFminsearch = optimset('TolX', 10^-6, 'TolFun', 10^-6, 'MaxIter', 600, 'MaxFunEvals',  600, 'Display', 'off'); %, 'Display','iter');

        % Store exit flag to inspect algo convergence
        [bestxFmin, errorNotToSave, exitFlag] = fminsearch(fun, initialGuess, optionsFminsearch)
        errorBestxFmin =  sseval3Pars(bestxFmin, x, saveGrowth(n, :), weights, detectionThreshold);
    

        bestx = bestxFmin;
        bestxError = errorBestxFmin;


        auxParForA(n) = bestx(1); % Save Lambda*/2
        auxParForB(n) = bestx(2); % Save 1/d*
        auxParForC(n) = bestx(3); % Save Lambda0
        auxParForE(n) = bestxError; % Save error score
        auxForExitFlag(n) = exitFlag; % Save exit flag; we look for exitflag = 1, means convergence


    end

    % Store dose-response parameters, error scores and convergence info
    savePars(:, 1) = auxParForA; savePars(:, 2) = auxParForB; savePars(:, 3) = auxParForC; saveError = auxParForE;  saveFlag= auxForExitFlag; 

end



% Store/save savePars as a nRepsx3 matrix to input in getErrorBars.m (for error analysis and parameters bounds)















