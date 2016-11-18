classdef Body < handle
    %Body class defines a rigid body and its associated makers
    properties
        ID;   %identifying number for the body
        
        %kinematics
        r;      %[x,y,z]' vector defining the location of the body in the G-RF
        rdot;   %d/dt[x,y,z]' 
        rddot;  %d2t/dt2[x,y,z]'
        p;      %[e0,e1,e2,e3]' vector defining the orientation of the body relative to G-RF
        pdot;   %time derivative of p
        pddot;  %2nd time derivative of p
        
        %mass properties
        m;       %[1x1] defining the mass of the ridgid body in kg
        I;       %[3x3] inertia tensor of the ridgid body
        
        %markers - used to define constraints
        bodyType;%string indicating if the body is 'ground' or 'free'
        markers = {};
    end
    
    properties (Dependent)
        A; %rotation matrix A of the body
        nMarkers; %total # of markers on the body
        isGround;
    end
    
    methods
        function body = Body(ID,bodyType) %default constructor
            body.ID = ID;
            body.bodyType = bodyType;
            body.markers{1} = [0,0,0]; %first marker is origin
            
            r = [0 0 0]';
            p = [0 0 0 0]';
            body.setKinematics(r,p);
        end
        
        function setKinematics(body,r,p,rdot,pdot,rddot,pddot)
            %handle underspecified kinematics gracefully
            if ~exist('rdot','var')
                rdot = [0 0 0]';
            end
            if ~exist('pdot','var')
                pdot =[0 0 0 0]';
            end
            if ~exist('rddot','var')
                rddot=[0 0 0]';
            end
            if ~exist('pddot','var')
                pddot = [0 0 0 0]';
            end
            
            %specified properties to body
            body.r = r;
            body.p = p;
            body.rdot = rdot;
            body.pdot = pdot;
            body.rddot =rddot;
            body.pddot = pddot;
            
            %perform dimension check
            body.dimCheck()
            
        end
        
%         function setMassProperties(m,I);
%            %future
%         end
            
        %add marker to list of markers
        function appendMarker(body, newMarker)
             %add a critical point(maker) that will be used
              %in the definition of contraints
            %input: 
                %body -instance of body class
                %marker [3x1] point in L-RF
            body.markers{length(body.markers) + 1} = newMarker; 
        end
        
        %getters
        function A = get.A(body)
            %returns [3x3] A matrix, calculated from the euler parameters P
            %[4x1] of the body
            A = MatOp.P2A(body.p);
        end
        
        function nMarkers = get.nMarkers(body)
            nMarkers = length(body.markers);
        end
        
        function B = Bp(body,markerInd)
            % input:
            %   markerInd - integer index of the marker number
            %   p (implicit) - [4x1]
            % output: 
            %   [3x4] B matrix
            B = MatOp.calcB(body.p, body.markers(markerInd));
        end
            
            
        function B = Bpdot(body,markerInd)
            % input:
            %   markerInd - integer index of the marker number
            %   pdot (implicit) - [4x1]
            % output: 
            %   [3x4] B matrix
            B = MatOp.calcB(body.pdot, body.markers(markerInd));
        end
            
        
        function isGround = ground(body)
            %by convension, the ground will be labeled with ID 0
            if body.ID == 0
                isGround = true;
            else
                isGround = false;
            end
        end
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
            if ~isequal(size(body.rdot),[3 1]) || ~isnumeric(body.rdot)
                error('rdot is not a [3x1]')
            end
            if ~isequal(size(body.pdot),[4 1]) || ~isnumeric(body.pdot)
                error('pdot is not a [4x1]')
            end
            if ~isequal(size(body.rddot),[3 1])|| ~isnumeric(body.rddot)
                error('rddot is not a [3x1]')
            end
            if ~isequal(size(body.pddot),[4 1])|| ~isnumeric(body.pddot)
                error('pddot is not a [4x1]')
            end
        end
    end
end