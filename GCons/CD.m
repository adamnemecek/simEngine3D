classdef CD
    %CD implements the CD constraint, using the r-p formulation
        %a constraint exists between two bodies, i and j
        
    
    properties
        
    end
    
    properties(Dependent)
        %these represent the outputs that a calling function
        %may request of the CD constraint (problem 2 prompt)
        %the request will be made with a getter
        
       
        nu;    %RHS of the velocity equation
        gamma; %RHS of the acceleration equation
        phi;   %constraint 
        phi_r; %partial derivative wrt r
        phi_p; %partial derivative wrt p
        
    end
    
    
    methods
        %constructor
        function CD = CDconst()
    end
    
end

