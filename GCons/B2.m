classdef B2
    %B2 constraint using r-p formulation. The B2 is an intermediate level
    %constraint which exists between a plane, formed by 2 vectors aiBar and
    %biBar on body i, and a vector dij from point Pi to point Qj
    %(note: see slide 24 from 9/26 for details and notation)
    %the constraint is made using 2 dp2 constraints, and removes 2 degrees
    %of freedom by enforcing 2 ACE's
  
    
    properties
        DOFremoved = 2;  %DOF removed
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        dp2_a;%instance of dp2 constraint between aiBar and dij
        dp2_b;%instance of dp2 constraint between biBar and dij
        
        %these point definitions are taken from S24 lecture 9.26
        aiHead; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        aiTail; %index of a [3x1] marker contained within bodyj, representing tail of aiBar
        biHead; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        biTail; %index of a [3x1] marker contained within bodyj, representing tail of aiBar 
        Pi;     %index of a [3x1] marker contained within bodyi, representing point Pi
        Qj;     %index of a [3x1] marker contained within bodyj, representing Qj
          
        
    end
    
    %calculated values 
    properties(Dependent)
        phi;   %constraint equation             [2x1] 
        nu;    %-dPhi/ dt                       [2x1]
        gamma; %RHS of acceleration equation    [2x1]
        
        phi_ri;%derivative of phi wrt ri        [2x3]  or [ ] if  i is ground
        phi_rj;%derivative of phi wrt rj        [2x3]  or [ ] if  j is ground
        phi_pi;%derivative of phi wrt ri        [2x4]  or [ ] if  i is ground
        phi_pj;%derivative of phi wrt ri        [2x4]  or [ ] if  j is ground

    end
    
    methods
        %constructor
        function b2 = B2(bodyi,bodyj,aiHead,aiTail,biHead,biTail,Pi,Qj)
            
            %bodies
            b2.bodyi = bodyi;
            b2.bodyj = bodyj;
            
            %vector head and tail indices
            b2.aiHead = aiHead;
            b2.aiTail = aiTail;
            b2.biHead = biHead;
            b2.biTail = biTail;
            b2.Pi  = Pi;
            b2.Qj = Qj;
           
            %form basic Gcons which make up B2 
            b2.dp2_a = DP2(b2.bodyi,b2.bodyj,b2.Pi,b2.Qj, b2.aiHead, b2.aiTail);
            b2.dp2_b = DP2(b2.bodyi,b2.bodyj,b2.Pi,b2.Qj, b2.biHead, b2.biTail);
            
            
        end
        
%--------------------getters for dependent properties-------------------------
        
        function phi = get.phi(b2)
            %phi [2x1] constraint equation
            phi = [b2.dp2_a.phi , b2.dp2_b.phi]';
        end
        
        function nu = get.nu(b2) 
           %nu [2x1] velocity equation [nu(1) , nu(2)]';
           nu = [b2.dp2_a.nu , b2.dp2_b.nu]';
        end
        
        function gamma = get.gamma(b2)
            %gamma [2x1] is the RHS of the acceleration equation
            gamma = [b2.dp2_a.gamma , b2.dp2_b.gamma]';
        end
       
      
        
        function phi_ri = get.phi_ri(b2)
            %phi_ri [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_ri is []
            phi_ri = [b2.dp2_a.phi_ri, b2.dp2_b.phi_ri]';
          
        end
        
         function phi_rj = get.phi_rj(b2)
            %phi_rj [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_rj is []
            phi_rj = [b2.dp2_a.phi_rj, b2.dp2_b.phi_rj]';
          
         end
        
         function phi_pi = get.phi_pi(b2)
            %phi_pi [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_pi is []
            phi_pi = [b2.dp2_a.phi_pi, b2.dp2_b.phi_pi]';
          
         end
        
          function phi_pj = get.phi_pj(b2)
            %phi_pj [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_pj is []
            phi_pj = [b2.dp2_a.phi_pj, b2.dp2_b.phi_pj]';
          
          end
        
       
    end
end
    

