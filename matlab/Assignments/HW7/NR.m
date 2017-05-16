classdef NR
    %NR is a static class, containing functions that are 
    %useful for solving Newton Rapson Equations for problem
    %4. particularly, calculating the jacobian and G matix
  
    methods(Static)
    
function RHS = G(h,Y_nm1, Y_n, alpha, beta)
    %input:
    %   h is the step size (scalar)
    %   Y_nm1 = 2x1 vector which is known.
    %   Y_n   = 2x1 vector which is being solved for
    %   alpha is a constant 
    %   beta is a constant 
    %output:
    
    %break out vectors into scalar values:
    x_nm1   = Y_nm1(1);
    y_nm1   = Y_nm1(2);
    x_n   = Y_n(1);
    y_n   = Y_n(2);
    
    %   G 2x1 vector representing RHS of newton rapson equation
    RHS = [x_n*(1+h) + 4*h*x_n/(1 + x_n^2) - x_nm1 - h*alpha , ...
          -h*beta*x_n + y_n + h*beta*x_n*y_n/(1+x_n^2) - y_nm1 ]';

end

function LHS = J(h,Y_n, beta)
    %input:
    %   Y_n = 2x1 vector which is being solved for
    %   beta is a constant 
    %output:
    %   J is the jacobian, which is part of the LHS of 
    %   the newton rapson matrix formlation
    
    %break out vectors into scalar values:
    x_n   = Y_n(1);
    y_n   = Y_n(2);
    den = 1 + x_n^2;  
    
    %form jacobian:
    LHS = [1 + h + 4*h*y_n*(1-x_n^2)/(den)^2 , 4*h*x_n/(den) ; ...
        -h*beta + h*beta*(1-x_n^2)/(den)^2 ,  1 + beta*h*x_n/(den)];
        
        
    end
    
end
end

