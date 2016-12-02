classdef System
    %System contains arrays of all bodies and constraints, as well as 
    % the global reference frame
    properties
        %setup global frame
        rGlobal = [0 0 0]'; %[3x1] global origin
        PGlobal = [1 0 0 0]'; %[4x1] global Euler parameters
        
        %cell arrays containing all bodies and constraints within system
        bodies = {}; %init with no bodies
        const = {};  
    end
    
    methods
        %constructor
        function sys = SysConst()
           sys.bodies = bodies;
           sys.const = const;
        end
        
    end
    
end

