classdef D
    %D or distance constraint using r-p formulation. The Distance constraint 
    %exists between a between a vector siBar(defined by Pi) on body i, and
    %a vector sjBar (defined by point Qj)on body j (see slide 16, lecture
    %9.26.2016 for a graphical representation). The constraint specifies 
    %that the distance between the heads of vectors aiBar and Ajbar maintain 
    %assume a certian value, F(t) which is frequently a constant,
    %rather than an explicit function of time
    
    properties
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        %these point definitions are taken from S10 lecture 9.2
        Pi; %integer index of a [3x1] marker contained within bodyi, representing head of siBar
        Qj; %integer index of a [3x1] marker contained within bodyj, representing head of sjBar
        f; %prescribed value that the dot product should assume
        fdot; %derivative of above
        fddot;%derivative of above
        t;
 
    end
    
    %calculated values 
    properties(Dependent)
        dij %vector in G-RF representing the distance from Pi to Qj
        dijdot;%derivative of above
        phi; %jacobian
        nu;  %  - d/dt(Phi)
        gamma; %RHS of acceleration equation
        phi_ri;%derivative of phi wrt ri
        phi_rj;%derivative of phi wrt rj
        phi_pi;%derivative of phi wrt pi
        phi_pj;%derivative of phi wrt pj   
    end
    
    
    methods
        %constructor
        function d = D(bodyi,bodyj,Pi,Qj,f,fdot,fddot,t)
            d.bodyi = bodyi;
            d.bodyj = bodyj;
            d.Pi = Pi;
            d.Qj = Qj;
            d.f = f;
            d.fdot = fdot;
            d.fddot = fddot;
            d.t = t;
            
        end
        
        
        %-----------getters for dependent properties-----------------------
        function dij = get.dij(d)
            %returns the distance vector [3x1] in G-RF that connects point
            %Pi on body i to point Qj on body j
            dij = d.bodyj.r + d.bodyj.A*d.bodyj.markers{d.Qj} - ...
                  (d.bodyi.r - d.bodyi.A*d.bodyi.markers{d.Pi});
        end
        
        function dijdot = get.dijdot(d)
            %returns the derivative of the distance vector [3x1]
            dijdot = d.bodyj.rdot + MatOp.B(d.bodyj.p, d.bodyj.marker{d.Qj})*d.bodyj.pdot ...
                    -d.bodyi.rdot + MatOp.B(d.bodyi.p, d.bodyi.marker{d.Pi})*d.bodyi.pdot;
        end
        
        function phi = get.phi(d) %9.26 slide 17
            %phi is the constraint equation [1x1] which states that 
            %the magnitude of the distance between Pi and Qj must be 
            %equal to ft
            dij = d.dij; %only calculate once
            phi =  dij'*dij - d.f;
        end
        
        function nu = get.nu(d) %9.26 slide 12
            nu = d.fdot;
        end
        
        function gamma = get.gamma(d)% 10.7 slide 8
            %gamma [1x1] is the RHS of the acceleration equation
            dij = d.dij; %calculate dij once per call
            dijdot = d.dijdot;
            bpdoti = MatOp.B(d.bodyi.p,d.bodyi.markers{d.Pi});
            bpdotj = MatOp.B(d.bodyj.p,d.bodyj.markers{d.Qj});
            
            gamma = -2*dij'*bpdotj*d.bodyj.pdot + 2*dij*bpdoti*d.bodyi.pdot + ...
                    -2*dijdot'*dijdot + d.fddot;   
        end
       
        
        function phi_ri = get.phi_ri(d) % s16 lecture 9.28.2016
            %phi_ri [1x3] is the derivative of phi with respect to the 
            %positional coordinates of bodyi. If I is the ground, then
            %phi_ri is [0x0]
            if d.bodyi.isground
                phi_ri = [];
            else
                phi_ri = -2*d.dij';
            end
        end
        
        function phi_rj = get.phi_rj(d) % s16 lecture 9.28.2016
            %phi_ri [1x3] is the derivative of phi with respect to the 
            %positional coordinates of bodyi. If I is the ground, then
            %phi_rj is [0x0]
            if d.bodyj.isground
                phi_rj = [];
            else
                phi_rj = 2*d.dij';
            end
        end
       
        
        function phi_pi = get.phi_pi(d) % s16 lecture 9.28.2016
            %phi_pi [1x4] is the derivative of phi with respect to the 
            %orientational coordinates of bodyi. If bodyi is the ground
            %then phi_pi is [0x0]
            if d.bodyi.isground
                phi_pi = [];
            else
                phi_pi = 2*d.dij'*-MatOp.B(d.bodyi.p,d.bodyi.markers{d.Pi});
            end
        end
        
        function phi_pj = get.phi_pj(d) % s16 lecture 9.28.2016
            %phi_pi [1x4] is the derivative of phi with respect to the 
            %orientational coordinates of bodyi. If bodyi is the ground
            %then phi_pi is [0x0]
            if d.bodyj.isground
                phi_pj = [];
            else
                phi_pj = 2*d.dij'*MatOp.B(d.bodyj.p,d.bodyj.markers{d.Qj});
            
            end
        end
           
            
    
    end
end
