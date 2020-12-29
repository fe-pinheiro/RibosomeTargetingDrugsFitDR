% Load WT and mutant data as matrices (from ensemble generated by FitDR.m)
% saveParsWT and saveParsMut are nReps*3 matrices, generated by FitDR.m and
% have the same output as savePars. 

% Explicitly, we input nReps rows with 3 columns and the following parameters inferred for wild-type and mutants:
% (a = Lambda/2; b = 1/d*, c = 1/Lambda0)

% The example given below illustrates the procedure. It requires input data for the wild-type 
% that is available from files ensembleParsWT.mat and for a mutant, here
% membrane mutant 1, in ensembleParsMM1.mat. Uncomment to run. .mat files
% for all other mutants are also available in the repository. The numbers
% correspon

% Fill here your the output from FitDR for wild-type and mutant of
% interest or uncomment below to run the example.
saveParsWT = [];
saveParsMut = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example:
% A = load('ensembleParsWT.mat'); 
% saveParsWT = A.saveParsWT;
% 
% B = load('ensembleParsMM1.mat'); 
% saveParsMut = B.saveParsMM1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nReps = length(saveParsWT); % We use the same number of replicates for wild-type and mutants
nReps2 = nReps^2;

% Scaled parameters and inference of model-dependent quantities (with propagation of error)
% Propagation of error for the scaled variables; Define vectors:
scaledLambda0 = zeros(nReps2,1); scaledResistance = zeros(nReps2,1);
saveScaledPin = zeros(nReps2,1); saveScaledPout = zeros(nReps2,1);
saveScaledLambdaStar = zeros(nReps2,1); saveScaledDStar = zeros(nReps2,1);

L0wt = zeros(nReps,1); D50wt = zeros(nReps,1); Ls_wt = zeros(nReps,1);
d50s_wt = zeros(nReps,1); pin_wt = zeros(nReps,1); pout_wt = zeros(nReps,1);

wildTypeID = zeros(nReps2, 1);
caseID = zeros(nReps2, 1);

z = 0;
for n1 = 1:nReps

    % Wild-type quantities..
    L0wt(n1) = (1/saveParsWT(n1, 3));
    D50wt(n1) = (0.5./(saveParsWT(n1, 2))).*( 2*saveParsWT(n1, 1).*saveParsWT(n1, 3) + 1./( 2*saveParsWT(n1, 1).*saveParsWT(n1, 3)));
    Ls_wt(n1) = 2*saveParsWT(n1, 1);
    d50s_wt(n1) = 1/saveParsWT(n1, 2);

    pin_wt(n1) = Ls_wt(n1)/d50s_wt(n1); 
    pout_wt(n1) = Ls_wt(n1)^2; 


    % Mutant quantities..
    for n2 = 1:nReps
        z = z+1;

        scaledLambda0(z) = (1/saveParsMut(n2, 3))/L0wt(n1);
        scaledResistance(z) = ((0.5./(saveParsMut(n2, 2))).*( 2*saveParsMut(n2, 1).*saveParsMut(n2, 3) + 1./( 2*saveParsMut(n2,1).*saveParsMut(n2, 3))))/D50wt(n1);   
        saveScaledPin(z) = 2*saveParsMut(n2, 1)*saveParsMut(n2, 2)/pin_wt(n1);
        
        saveScaledPout(z) = (2*saveParsMut(n2, 1))^2/pout_wt(n1);  % For Extended Data Table 2, column 7  
        saveScaledLambdaStar(z) = (2*saveParsMut(n2, 1))/(Ls_wt(n1)); % For Extended Data Table 2, column 5 
        saveScaledDStar(z) = (1/saveParsMut(n2, 2))/d50s_wt(n1); % For Extended Data Table 2, column 4 
        
  
        % In case we would like to track the wild-type
        wildTypeID(z) =  n1;
        caseID(z) = n2;

    end 
end

    
%% Error bars for different observables used in Fig. 3
% Observables/parameter bounds are given as 90-10 percentiles; 
% The following observables correspond to columns 3, 4 and 6 of Extended Data Table 2

ub = floor(0.9*nReps^2); % Upper bound
lb = floor(0.1*nReps^2); % Lower bound

% Bounds of scaled resistance
sortScaledResistance = sort(scaledResistance);
lowerScaledResistance = sortScaledResistance(lb);
upperScaledResistance = sortScaledResistance(ub);
fprintf('\n Resistance bounds:\n Lower: %f\n Upper %f\n', lowerScaledResistance, upperScaledResistance)


% Bounds of scaled drug-free growth rate
sortScaledLambda0 = sort(scaledLambda0);
lowerLambda0 = sortScaledLambda0(lb);
upperLambda0 = sortScaledLambda0(ub);
fprintf('\n Drug-free growth rate bounds:\n Lower: %f\n Upper %f\n', lowerLambda0, upperLambda0)


% Bounds of Pin
sortScaledPin = sort(saveScaledPin);
lowerPin = sortScaledPin(lb);
upperPin = sortScaledPin(ub);
fprintf('\n Pin bounds:\n Lower: %f\n Upper %f\n', lowerPin, upperPin)












