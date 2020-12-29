function [ Y ] = makeG(pars, xData)

% This function computes the model dosage-response curve and is used for
% plotting
%
% Input the parameters for the metabolic drug response as described in
% loadFunctionsGreulichModel

a = pars(1); % Lambda*/2
b = pars(2); % 1/d*
c = pars(3);  %1/L0

% Load function definitions..
loadFunctionsGreulichModel;


 
%  The next step deals with jumps between branches 
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

end

