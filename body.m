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
        
        %getters
        function A = get.A()
            A = MatOp.A2P(body.p);
        end
        
        
        methods(Access = private)
            
        %body static type and dimension checking
        function dimCheck(body)
            if ~isequal(size(body.ID),[1 1])   || ~isnumeric(body.ID)
                error('ID is not a number')
            end
            if ~isequal(size(body.r),[3 1])    || ~isnumeric(body.r)
                error('r is not a [3x1] matrix')
            end
            if ~isequal(size(body.p),[4 1])    || ~isnumeric(body.p)
                error('p is not a [4x1]')
            end
            if ~isequal(size(body.rdot),[4 1]) || ~isnumeric(body.rdot)
                error('rdot is not a [4x1]')
            end
            if ~isequal(size(body.pdot),[4 1]) || ~isnumeric(body.pdot)
                error('pdot is not a [4x1]')
            end
            if ~isequal(size(body.rddot),[4 1])|| ~isnumeric(body.rddot)
                error('rddot is not a [4x1]')
            end
            if ~isequal(size(body.pddot),[4 1])|| ~isnumeric(body.pddot)
                error('pddot is not a [4x1]')
            end
            
            
    end
    
end
