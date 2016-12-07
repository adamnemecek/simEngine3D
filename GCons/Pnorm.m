classdef Pnorm
    %Pnorm is a constraint using r-p formulation. Pnorm is the Euler
    %Parameter (p) normalization constraint which statest that the dot
    %product of P with itself (p'* p) must be equal to 1. this is a
    %consiquence of the way that P is defined, and the trig identity,
    %sinx^2 + cosx^2 = 1. 
    
    %consult Haug pg. 340 , or lecture notes S13 , 9.21 for derivation
  
    
    properties
        DOFremoved = 1;  %DOF removed
        bodyi; %instance of a body object, i 
        
 
    end
    
    %calculated values 
    properties(Dependent)
        phi;   %constraint equation             [1x1] 
        nu;    %-dPhi/ dt                       [1x1]
        gamma; %RHS of acceleration equation    [1x1]
        
        phi_ri;%derivative of phi wrt ri        [1x3]  or [ ] if  i is ground
        phi_pi;%derivative of phi wrt ri        [1x4]  or [ ] if  i is ground
 

    end
    
    methods
        %constructor
        function pnorm = Pnorm(bodyi)
            pnorm.bodyi = bodyi;
        end
        
%--------------------getters for dependent properties-------------------------
        
        function phi = get.phi(pnorm) %S24 9.28
            %phi [1x1] constraint equation
            phi = pnorm.bodyi.p'*pnorm.bodyi.p - 1;
           
        end
        
        function nu = get.nu(pnorm) %S25 9.28
           %nu [1x1] velocity equation
           nu = 0; 
        end
        
        function gamma = get.gamma(pnorm)%s25 9.28
            %gamma [1x1] is the RHS of the acceleration equation
            gamma = -2*pnorm.bodyi.pdot'*pnorm.bodyi.pdot;
        end
       
      
        
        function phi_ri = get.phi_ri(pnorm) %haug pg. 383
            %phi_ri [1x3] derivative of phi wrt location parameters
            phi_ri = 0; %this agrees with our intuition
          
        end
       
        function phi_pi = get.phi_pi(pnorm) %haug pg. 383
            %phi_pi [1x4] derivative of phi wrt orientation parameters
            phi_pi = 2*pnorm.bodyi.p';
           
          
        end
    end
end
    

