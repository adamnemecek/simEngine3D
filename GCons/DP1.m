classdef DP1
    %DP1 constraint using r-p formulation. The DP1 constraint exists
    %between a vector aiBar on body i, and a vector ajBar, on body j. DP1
    %constraint specifies that the dot product of the two vectors aiBar
    %and ajBar assume a certian value, F(t), which is frequently a constant
    
    properties
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        
        %these point definitions are taken from S10 lecture 9.26
        Qi; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        Pi; %index of a [3x1] marker contained within bodyi, representing tail of aiBar
        Qj; %index of a [3x1] marker contained within bodyj, representing head of ajBar
        Pj; %index of a [3x1] marker contained within bodyj, representing tail of ajBar
        f; %prescribed value that the dot product should assume
        fdot; %derivative of above
        fddot;%derivative of above 
        t;  
    end
    
    %calculated values 
    properties(Dependent)
        aiBar; % vector in local reference frame of body i
        ajBar; % vector in local reference frame of body j
        phi; %jacobian
        nu;  % dPhi_dp1/ dt
        gamma; %RHS of acceleration equation
        phi_ri;%derivative of phi wrt ri
        phi_rj;%derivative of phi wrt rj
        phi_pi;%derivative of phi wrt ri
        phi_pj;%derivative of phi wrt ri
        
        
        
    end
    
    methods
        %constructor
        function dp1 = DP1(bodyi,bodyj,Qi,Pi,Qj,Pj,f,fdot,fddot,t)
            dp1.bodyi = bodyi;
            dp1.bodyj = bodyj;
            dp1.Qi = Qi;
            dp1.Pi = Pi;
            dp1.Qj = Qj;
            dp1.Pj = Pj;
            dp1.f = f;
            dp1.fdot = fdot;
            dp1.fddot = fddot;
            dp1.t = t;
            
        end
        %getters for dependent properties
        
        function aiBar = get.aiBar(dp1) %9.26 slide 10
            %ai vector in the local reference frame
            aiBar = dp1.bodyi.markers{dp1.Qi}-dp1.bodyi.markers{dp1.Pi};
        end
        
        function ajBar = get.ajBar(dp1) %9.26 slide 10
            ajBar = dp1.bodyj.markers{dp1.Qj}-dp1.bodyi.markers{dp1.Pj};
        end
        
        function phi = get.phi(dp1) %9.26 slide 11
            %phi = aiBar'* Ai'* Aj * ajBar
            phi = (dp1.aiBar' * dp1.bodyi.A')*(dp1.bodyj.A * dp1.ajBar) - dp1.f;
        end
        
        function nu = get.nu(dp1) %9.26 slide 12
            nu = dp1.fdot;
        end
        
        function gamma = get.gamma(dp1)% 10.7 slide 8
            %gamma [1x1] is the RHS of the acceleration equation
            ai = dp1.bodyi.A*dp1.aiBar; % a = A*abar
            bpdotj = Matop.B(dp1.bodyj.pdot, dp1.ajBar);
            pdotj = dp1.bodyj.pdot;
            aj = dp1.bodyj.A*dp1.ajBar;
            bpdoti = Matop.B(dp1.bodyi.pdot, dp1.aiBar);
            pdoti = dp1.bodyi.pdot;
            adoti = Matop.B(dp1.bodyi.p, dp1.aiBar)*dp1.bodyi.pdot; %aidot = B(p,abar)*pdot
            adotj = Matop.B(dp1.bodyj.p, dp1.ajBar)*dp1.bodyj.pdot;
            
            
            gamma = -ai'*bpdotj*pdotj - aj'*bpdoti*pdoti - 2*adoti'*adotj + d1.fddot;
            
        end
       
        function phi_ri = get.phi_ri(dp1) % s13 lecture 9.28.2016
            %phi_ri [1x3] is the derivative of phi with respect to the 
            %positional coordinates of bodyi. If I is the ground, then
            %phi_ri is [0x0]
            if dp1.bodyi.isground
                phi_ri = [];
            else
                phi_ri = [0 0 0]
            end
        end
        
        function phi_rj = get.phi_rj(dp1) % s13 lecture 9.28.2016
            %phi_ri [1x3] is the derivative of phi with respect to the 
            %positional coordinates of bodyi. If I is the ground, then
            %phi_rj is [0x0]
            if dp1.bodyj.isground
                phi_rj = [];
            else
                phi_rj = [0 0 0];
            end
        end
       
        
        function phi_pi = get.phi_pi(dp1) % s13 lecture 9.28.2016
            %phi_pi [1x4] is the derivative of phi with respect to the 
            %orientational coordinates of bodyi. If bodyi is the ground
            %then phi_pi is [0x0]
            if dp1.bodyi.isground
                phi_pi = []
            else
                phi_pi = (dp1.bodyj.A*dp1.ajBar)'*MatOp.B(dp1.bodyi.p, dp1.aiBar);
            end
        end
        
        function phi_pj = get.phi_pj(dp1) % s13 lecture 9.28.2016
            %phi_pi [1x4] is the derivative of phi with respect to the 
            %orientational coordinates of bodyi. If bodyi is the ground
            %then phi_pi is [0x0]
            if dp1.bodyj.isground
                phi_pj = [];
            else
                phi_pj =(dp1.bodyi.A*dp1.aiBar)'*MatOp.B(dp1.bodyj.p, dp1.ajBar);
            end
            
            end
        end
           
            
    
end
