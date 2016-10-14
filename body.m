classdef Body
    %Body class defines a rigid body and its associated makers
   
    
    properties
        ID; %identifying number for the body
        
        %kinematics
        r; % [x,y,z] vector defining the location of the body in the global RF
        p; % [e0,e1,e2,e3] vector defining the oriation
        pdot; %time derivative of p
        
        %mass properties
        m; %[1x1] defining the mass of the ridgid body in kg
        I; %[3x3] inertia tensor of the ridgid body
        
        %markers - used to define constraints
        markers; % cell array of markers, accessed by their index 
                 % ex) body.markers{1} -> [0,0,0]
                 %first marker should always be the origin
        
    end
    
    methods
        %default constructor
        function body = constBody()
        end
        
    end
    
end

