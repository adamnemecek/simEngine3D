classdef SJ
    %SJ (spherical)constraint using r-p formulation. The SJ is an High level
    %constraint which specifies that a point located on body i, and a point
    %located on body j remain coincident at all times. 
    
    %(note: see slide 26 from 9/26 for details and notation)
    
    %the constraint is made using 3 CD constraints, and removes 3 degrees
    %of freedom by enforcing 3 ACE's
  
    
    properties
        DOFremoved = 3;  %DOF removed
        bodyi; %instance of a body object, i 
        bodyj; %instance of a body object, j
        
        cdi;%instance of the CD constraint, specifying coincidence of x constraint in G-RF 
        cdj;%instance of the CD constraint, specifying coincidence of y constraint in G-RF 
        cdk;%instance of the CD constraint, specifying coincidence of z constraint in G-RF 
        
        %these point definitions are taken from S22 lecture 9.26
        Pi; %index of a [3x1] marker contained within bodyi, representing head of sIbarP
        Qj; %index of a [3x1] marker contained within bodyj, representing head of sJbarP
       
    end
    
    %calculated values 
    properties(Dependent)
        phi;   %constraint equation             [3x1] 
        nu;    %-dPhi/ dt                       [3x1]
        gamma; %RHS of acceleration equation    [3x1]
        
        phi_ri;%derivative of phi wrt ri        [3x3]  or [ ] if  i is ground
        phi_rj;%derivative of phi wrt rj        [3x3]  or [ ] if  j is ground
        phi_pi;%derivative of phi wrt ri        [3x4]  or [ ] if  i is ground
        phi_pj;%derivative of phi wrt ri        [3x4]  or [ ] if  j is ground

    end
    
    methods
        %constructor
        function sj = SJ(bodyi,bodyj,PiIndex,QjIndex)
            
            %bodies
            sj.bodyi = bodyi;
            sj.bodyj = bodyj;
            
            %vector head and tail indices
            sj.Pi = PiIndex;
            sj.Qj = QjIndex;
            
            %form basic Gcons which make up SJ Slide 
            sj.cdi = CD(sj.bodyi,sj.bodyj,sj.Pi, sj.Qj,'x');
            sj.cdj = CD(sj.bodyi,sj.bodyj,sj.Pi, sj.Qj,'y');
            sj.cdk = CD(sj.bodyi,sj.bodyj,sj.Pi, sj.Qj,'z');
            
            
        end
        
%--------------------getters for dependent properties-------------------------
        
        function phi = get.phi(sj)
            %phi [3x1] constraint equation
            phi = [sj.cdi.phi , sj.cdj.phi , sj.cdk.phi]';
        end
        
        function nu = get.nu(sj) 
           %nu [3x1] velocity equation [nu(1) , nu(2)]';
           nu = [sj.cdi.nu , sj.cdj.nu, sj.cdk.nu]';
        end
        
        function gamma = get.gamma(sj)
            %gamma [3x1] is the RHS of the acceleration equation
            gamma = [sj.cdi.gamma , sj.cdj.gamma, sj.cdk.gamma]';
        end
       
      
        
        function phi_ri = get.phi_ri(sj)
            %phi_ri [3x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_ri is []
            phi_ri = [sj.cdi.phi_ri, sj.cdj.phi_ri, sj.cdk.phi_ri]';
          
        end
        
         function phi_rj = get.phi_rj(sj)
            %phi_rj [3x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_rj is []
             phi_rj = [sj.cdi.phi_rj, sj.cdj.phi_rj, sj.cdk.phi_rj]';
          
         end
        
         function phi_pi = get.phi_pi(sj)
            %phi_pi [3x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_pi is []
            phi_pi = [sj.cdi.phi_pi, sj.cdj.phi_pi, sj.cdk.phi_pi]';
          
         end
        
          function phi_pj = get.phi_pj(sj)
            %phi_pj [2x3] is the derivative of both phi equations, with respect
            %to the positional coordinates of body i, if body i is ground,
            %phi_pj is []
            phi_pj = [sj.cdi.phi_pj, sj.cdj.phi_pj, sj.cdk.phi_pj]';
          
          end
        
       
    end
end
    

