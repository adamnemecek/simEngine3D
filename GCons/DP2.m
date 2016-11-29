classdef DP2
    %dp2 constraint using r-p formulation. The dp2 constraint exists
    %between a vector aiBar on body i, and a vector dij which is the distance 
    %between point Pi and Point Qj, represented in the G-RF. dp2
    %constraint specifies that the dot product of these two vectors assume a
    %certian value, F(t), which is frequently a constant
    
    properties
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        
        %these point definitions are taken from S13 lecture 9.26
        Pi; %index of a [3x1] marker contained within bodyi, representing head of SiBar
        Qj; %index of a [3x1] marker contained within bodyj, representing head of SjBar
        ai; %index of a [3x1] marker contained within bodyi, representing the head of AiBar
        f; %prescribed value that the dot product should assume
        fdot; %derivative of above
        fddot;%derivative of above 
        t;  
    end
    
    %calculated values 
    properties(Dependent)
        aiBar; %vector AiBar the first vector in the DP2 constraint
        dij;   %distance in G-RF from point Pi to point Qj, the second vector in the DP2 constraint
        dijdot;%derivative of dij
        phi; %jacobian
        nu;  % dPhi_dp2/ dt
        gamma; %RHS of acceleration equation
        phi_ri;%derivative of phi wrt ri
        phi_rj;%derivative of phi wrt rj
        phi_pi;%derivative of phi wrt ri
        phi_pj;%derivative of phi wrt ri
        
        
        
    end
    
    methods
        %constructor
        function dp2 = DP2(bodyi,bodyj,Pi,Qj,ai,f,fdot,fddot,t)
            dp2.bodyi = bodyi;
            dp2.bodyj = bodyj;
            dp2.Pi = Pi;
            dp2.Qj = Qj;
            dp2.ai = ai;
            dp2.f = f;
            dp2.fdot = fdot;
            dp2.fddot = fddot;
            dp2.t = t;
            
        end
        %getters for dependent properties
        
        %consider revision that makes dij an abstract method
          function dij = get.dij(dp2)
            %returns the distance vector [3x1] in G-RF that connects point
            %Pi on body i to point Qj on body j
            dij = dp2.bodyj.r + dp2.bodyj.A*dp2.bodyj.markers{dp2.Qj} - ...
                 (dp2.bodyi.r - dp2.bodyi.A*dp2.bodyi.markers{dp2.Pi});
        end
        
        function dijdot = get.dijdot(dp2)
            %returns the derivative of the distance vector [3x1]
            dijdot = dp2.bodyj.rdot + MatOp.B(dp2.bodyj.p, dp2.bodyj.marker{dp2.Qj})*dp2.bodyj.pdot ...
                    -dp2.bodyi.rdot + MatOp.B(dp2.bodyi.p, dp2.bodyi.marker{dp2.Pi})*dp2.bodyi.pdot;
        end  
        
        function aiBar = get.aiBar(dp2)
            % [3x1] local vector
            aiBar = dp2.bodyi.markers{dp2.ai};
        end
       
        function phi = get.phi(dp2) %9.26 slide 14
            %DP2 [1x1] constraint equation
            phi = (dp2.aiBar' * dp2.bodyi.A') * dp2.dij - dp2.f;
        end
        
        function nu = get.nu(dp2) %9.26 slide 12
            nu = dp2.fdot;
        end
        
        function gamma = get.gamma(dp2)% 10.7 slide 8
            %gamma [1x1] is the RHS of the acceleration equation
            Ai      = dp2.bodyi.A*dp2.aiBar; % a = A*abar
            bpdotj  = MatOp.B(dp2.bodyj.pdot, dp2.bodyj.markers{dp2.Qj});
            pdotj   = dp2.bodyj.pdot;
            bpdoti  = MatOp.B(dp2.bodyi.pdot, dp2.bodyi.markers{dp2.Pi});
            pdoti   = dp2.bodyi.pdot;
            bpdotai = MatOp.B(dp2.bodyi.pdot, dp2.aiBar);
            adoti   = MatOp.B(dp2.bodyi.p, dp2.aiBar)*dp2.bodyi.pdot; %aidot = B(p,abar)*pdot
           
            gamma = -Ai'*bpdotj*pdotj + Ai'*bpdoti*pdoi - dp2.dij'*bpdotai*pdoti - 2*adoti'*dp2.dijdot + dp2.fddot;
            
        end
       
        function phi_ri = get.phi_ri(dp2) % s13 lecture 9.28.2016
            %phi_ri [1x3] is the derivative of phi with respect to the 
            %positional coordinates of bodyi. If I is the ground, then
            %phi_ri is [0x0]
            if dp2.bodyi.isground
                phi_ri = [];
            else
                phi_ri = -dp2.aiBar';
            end
        end
        
        function phi_rj = get.phi_rj(dp2) % s13 lecture 9.28.2016
            %phi_ri [1x3] is the derivative of phi with respect to the 
            %positional coordinates of bodyi. If I is the ground, then
            %phi_rj is [0x0]
            if dp2.bodyj.isground
                phi_rj = [];
            else
                phi_rj = dp2.aiBar';
            end
        end
       
        
        function phi_pi = get.phi_pi(dp2) % s13 lecture 9.28.2016
            %phi_pi [1x4] is the derivative of phi with respect to the 
            %orientational coordinates of bodyi. If bodyi is the ground
            %then phi_pi is [0x0]
            if dp2.bodyi.isground
                phi_pi = [];
            else
                phi_pi = dp2.dij*MatOp.B(dp2.bodyi.p,dp2.aiBar) - ...
                         dp2.aiBar'*MatOp.B(dp2.bodyi.p,dp2.bodyi.markers{dp2.Pi}) ;
            end
        end
        
        function phi_pj = get.phi_pj(dp2) % s13 lecture 9.28.2016
            %phi_pi [1x4] is the derivative of phi with respect to the 
            %orientational coordinates of bodyi. If bodyi is the ground
            %then phi_pi is [0x0]
            if dp2.bodyj.isground
                phi_pj = [];
            else
                 phi_pj = dp2.aiBar'*MatOp.B(dp2.bodyj.p,dp2.bodyj.markers{dp1.Qj});
            end
            end
            
    end
end
           
            
    

