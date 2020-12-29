function minimizeF = sseval2ParsLsIC50(pars, const, xData, yData, weights, detectionThreshold)

% Function to fit the metabolic model with the constraint: IC50*Lambda* = 1/const
% Inputs: 
% (1) const
% (2) xData - takes the drug values
% (3) yData - takes the measured growth at the drug values given in xData
% (4) weigths - takes the weights used in the fitting procedure; we use inverse variance of measurements across 3 repliicates
% (5) detectionThreshold - takes the value of the detection threshold

a = pars(1); % Lambda*/2
b = const*(2*pars(1));  % 1/d*
c = pars(3); % 1/Lambda0

Y = nan(length(xData), 1); % Vector allocated to store growth values

w_i = ones(length(xData), 1)';
if ~isempty(weights)
w_i = weights;
end


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


% Compute error score

nDataPoints = length(yData);
thetaError = ones(1, nDataPoints);

for ii = 1:nDataPoints
    
    if (yData(ii) <detectionThreshold && Y(ii) < detectionThreshold)
       thetaError(ii) = 0;
    end
end


 minimizeF  = thetaError.*abs(yData - Y').*sqrt(w_i) + 100*(heaviside(-c) + heaviside(-a) + heaviside(a - 1) +  heaviside(b - 0.25) + heaviside(-b)); 

 % The term 100*(..) is added to penalize eventual non-physical parameter values that might be found by the fitting procedure in pathological cases 
 
