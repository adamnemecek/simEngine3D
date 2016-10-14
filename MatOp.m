classdef MatOp
    % MatOps collection of useful matrix helper functions
    %   Detailed explanation goes here

    methods
        function a_tilde = tilde(a) 
            %takes a 1x3 input vector and returns a 3x3
            %matrix, capable of performing the cross product
            a_tilde = [ 0 -a(3) a(2);
                      a(3)  0  -a(1);
                     -a(2) a(1)  0];
        end
        
        %conversion functions from A2P and P2A
        function  P = A2P(A)
        end
        
        function A = P2A(P)
        end
        end
        
    end
    
end

