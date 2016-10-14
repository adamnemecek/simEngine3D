classdef Body < handle
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
        markers{1} = [0,0,0]'; % cell array of markers, accessed by their index 
                 % ex) body.markers{1} -> [0,0,0]
                 %first marker should always be the origin
        
    end
    
    methods
        %default constructor
        function body = makeBody(ID,r,p,pdot,m,I,markers)
            body.ID = ID;
            body.r = r;
            body.p = p;
            body.pdot = pdot;
            body.m = m;
            body.I = I;
            body.markers = markers;
        end
            
        %add marker to list of markers
        function appendMarker(body, newMarker)
            %input: body object, marker [3x1]
            body.markers{length(markers) + 1} = newMarker; 
        end
            
    end
    
end
