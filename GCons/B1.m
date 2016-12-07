classdef B1
    %B1 constraint using r-p formulation. The B1 is an intermediate level
    %constraint which exists between a plane, formed by 2 vectors aiBar and
    %biBar on body i, and a vector cjBar on body j.
    %(note: see slide 22 from 9/26 for details and notation)
    %the constraint is made using 2 dp1 constraints, and removes 2 degrees
    %of freedom by enforcing 2 ACE's
  
    
    properties
        DOFremoved = 2;  %DOF removed
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        dp1_ac;%instance of dp1 constraint between aiBar cjBar
        dp1_bc;%instance of dp1 constraint between biBar cjBar 
        
        %these point definitions are taken from S22 lecture 9.26
        aiHead; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        aiTail; %index of a [3x1] marker contained within bodyj, representing tail of aiBar
        biHead; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        biTail; %index of a [3x1] marker contained within bodyj, representing tail of aiBar 
        cjHead; %index of a [3x1] marker contained within bodyi, representing head of aiBar
        cjTail; %index of a [3x1] marker contained within bodyj, representing tail of aiBar 
 
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
        function b1 = B1(bodyi,bodyj,aiHead,aiTail,biHead,biTail,cjHead,cjTail)
            
            %bodies
            b1.bodyi = bodyi;
            b1.bodyj = bodyj;
            
            %vector head and tail indices
            b1.aiHead = aiHead;
            b1.aiTail = aiTail;
            b1.biHead = biHead;
            b1.biTail = biTail;
            b1.cjHead = cjHead;
            b1.cjTail = cjTail;
            
            %form basic Gcons which make up B1 Slide 
            b1.dp1_ac = DP1(b1.bodyi,b1.bodyj,b1.aiHead,b1.aiTail, b1.cjHead, b1.cjTail);
            b1.dp1_bc = DP1(b1.bodyi,b1.bodyj,b1.biHead,b1.biTail, b1.cjHead, b1.cjTail);
            
            
        end
        
%--------------------getters for dependent properties-------------------------
        
        function phi = get.phi(b1)
            %phi [2x1] constraint equation
            phi = [b1.dp1_ac.phi , b1.dp1_bc.phi]';
        end
        
        function nu = get.nu(b1) 
           %nu [2x1] velocity equation [nu(1) , nu(2)]';
           nu = [b1.dp1_ac.nu , b1.dp1_bc.nu]';
        end
        
        function gamma = get.gamma(b1)
            %gamma [2x1] is the RHS of the acceleration equation
            gamma = [b1.dp1_ac.gamma , b1.dp1_bc.gamma]';
        end
       
      
        
        function phi_ri = get.phi_ri(b1)
            %phi_ri [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_ri is []
            phi_ri = [b1.dp1_ac.phi_ri, b1.dp1_bc.phi_ri]';
          
        end
        
         function phi_rj = get.phi_rj(b1)
            %phi_rj [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_rj is []
            phi_rj = [b1.dp1_ac.phi_rj, b1.dp1_bc.phi_rj]';
          
         end
        
         function phi_pi = get.phi_pi(b1)
            %phi_pi [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_pi is []
            phi_pi = [b1.dp1_ac.phi_pi, b1.dp1_bc.phi_pi]';
          
         end
        
          function phi_pj = get.phi_pj(b1)
            %phi_pj [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_pj is []
            phi_pj = [b1.dp1_ac.phi_pj, b1.dp1_bc.phi_pj]';
          
          end
        
       
    end
end
    

