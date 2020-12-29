function sseWithThreshold = sseval3Pars(pars, xData, yData, weights, detectionThreshold)

% Function to fit the metabolic model with 3 parameters
% Inputs: 
% (1) xData - takes the drug values
% (2) yData - takes the measured growth at the drug values given in xData
% (3) weigths - takes the weights used in the fitting procedure; we use inverse variance of measurements across 3 replicates
% (4) detectionThreshold - takes the value of the detection threshold

a = pars(1); % Lambda*/2
b = pars(2); % 1/d*
c = pars(3); % 1/Lambda0


% Add option of not using weights..
w_i = ones(length(xData), 1);
if ~isempty(weights)
    w_i = weights;
end

Y = nan(length(xData), 1); % Vector allocated to store growth values

% Load functions of metabolic model and deal with branch changes
loadFunctionsGreulichModel;

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
        absError = detectionThreshold - Y(i);
    end
    
    
    if yData(i) < detectionThreshold && Y(i) < detectionThreshold
       thetaError = 0;
    end
    auxSSE = auxSSE + abs(thetaError*(w_i(i))*(absError).^2);

end


sseWithThreshold = (auxSSE)  + 100*(heaviside(a-1) + heaviside(-a) + heaviside(-b)); 

 % The term 100*(..) is added to penalize eventual non-physical parameter values that might be found by the fitting procedure in pathological cases 
