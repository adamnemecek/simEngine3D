classdef DP1
    %DP1 constraint using r-p. dp acts between two markers, from two 
    %separate ridgid bodies
    %   Detailed explanation goes here
    
    properties
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        Qi; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        Pi; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        Qj; %index of a [3x1] marker contained within bodyj, representing head of ajBar
        Pj; %index of a [3x1] marker contained within bodyj, representing head of ajBar
        f; %prescribed value that the dot product should assume
        fdot; %derivative of above
        fddot;%derivative of above       
    end
    
    %calculated values 
    properties(Dependent)
        aiBar; % vector in local reference frame of body i
        ajBar; % vector in local reference frame of body j
        phi; %jacobian
        nu;  % dPhi_dp1/ dt
        gamma; % RHS of acceleration equation
        phiR;
        phiP;
        
        
        
    end
    
    methods
        %constructor
        function dp1 = DP1(bodyi,bodyj,Qi,Pi,Qj,Pj,f,fdot,fddot)
            dp1.bodyi = bodyi;
            dp1.bodyj = bodyj;
            dp1.Qi = Qi;
            dp1.Pi = Pi;
            dp1.Qj = Qj;
            dp1.Pj = Pj;
            dp1.f = f;
            dp1.fdot = fdot;
            dp1.fddot = fddot;
            
        
        %getters for dependent properties
        function aiBar = get.aiBar() %9.26 slide 10
            aiBar = dp1.bodyi.markers{Qi}-dp1.bodyi.markers{Pi};
        end
        
        function ajBar = get.ajBar() %9.26 slide 10
            ajBar = dp1.bodyj.markers{Qj}-dp1.bodyi.markers{Pj};
        end
        
        function phi = get.phi() %9.26 slide 11
            %phi = aiBar'* Ai'* Aj * ajBar
            phi = (aiBar' * dp1.bodyi.A')*(dp1.bodyj.A * ajBar) - dp1.f;
        end
        
        function nu = get.nu() %9.26 slide 12
            nu = dp1.fdot;
        end
        
        function gamma = get.gamma()% 9.26 slide 12
            gamma = 
            
        end
       
        function phiR = phiR()
        end
        
        function phiP = phiP()
        end
            
    end
    
end

