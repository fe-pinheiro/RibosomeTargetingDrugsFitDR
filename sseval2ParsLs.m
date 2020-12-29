function sseWithSolver = sseval2ParsLs(pars, a, xData, yData, weights, detectionThreshold)

% Function to fit the metabolic model with constrained Lambda* 
% Inputs: 
% (1) a - takes Lambda*/2
% (2) xData - takes the drug values
% (3) yData - takes the measured growth at the drug values given in xData
% (4) weigths - takes the weights used in the fitting procedure; we use inverse variance of measurements across 3 replicates
% (5) detectionThreshold - takes the value of the detection threshold

b = pars(1); % 1/d*
c = pars(2);  % 1/Lambda0

% Load function definitions..
loadFunctionsGreulichModel;

w_i = ones(length(xData), 1);

if ~isempty(weights)
w_i = weights;
end


% Compute values to of the growth function and deal with branching points..

Y = nan(length(xData), 1);
kk = 0;
for X = xData
           
           kk = kk+1;
           
            f = 24*sqrt(3)*sqrt(-(a^2)*(b^12)*(c^2) + 81*(a^4)*(b^12)*(c^4) - 2187*(a^6)*(b^12)*(c^6) +19683*(a^8)*(b^12)*(c^8));

            g  = b^6 - 540*(a^2)*(b^6)*(c^2) - 5832*(a^4)*(b^6)*(c^4);

           drugSecondBranch = -((-b^2 + 12*(a^2)*(b^2)*(c^2))/(12*a*(b^3)*c)) - (-729*(b^4) - 157464*(a^2)*(b^4)*(c^2))/(8748*a*(b^3)*c*(g + f)^(1/3)) + (g + f)^(1/3)/(12*a*(b^3)*c);

           
           drugThirdBranch = (1 - 3*(a^2)*(c^2))/(3*a*b*c);
           
           
           if X> drugSecondBranch & X< drugThirdBranch
                Y(kk) = real(g3(a, b, c , X));
           else
               Y(kk) = real(g1(a, b, c , X));
           end
                     
end


% Make error score

auxSSE = 0;
thetaError = 1;
absError = 0;

for i = 1:length(yData)
   
    absError = yData(i) - Y(i);
    if  yData(i) == 0
        absError = detectionThreshold/2 - Y(i);
    end

    
    if yData(i) < detectionThreshold && Y(i) < detectionThreshold
       thetaError = 0;
    end
    
    auxSSE = auxSSE + thetaError*w_i(i)*(absError)^2;
  
end

sseWithSolver = auxSSE +  heaviside(a - 1) + heaviside(-a) + heaviside(-b) +  heaviside(b - 0.25) + heaviside(-c);

% Heaviside functions are added to avoid unphysical solutions
