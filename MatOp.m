classdef MatOp
    % MatOps collection of useful matrix helper functions
    %   Detailed explanation goes here

    methods(Static) %don't require instantiation
        function a_tilde = tilde(a) 
            %takes a 1x3 input vector and returns a 3x3
            %matrix, capable of performing the cross product
            a_tilde = [ 0 -a(3) a(2);
                      a(3)  0  -a(1);
                     -a(2) a(1)  0];
        end
        
        %conversion functions from A2P and P2A
        function  P = A2P(A)
            %to be implemented
        end
        
        function A = P2A(P)
            e0 = p(1);
            e = p(2:4);
            %from key kinematic formulas
            E = [-e,  MatOp.tilde(e) +e0*eye(3)];
            G = [-e, -MatOp.tilde(e) +e0*eye(3)];
            A = EG';
        end
        
        function B = calcB(p,aBar) %9.28 slide 12
            %input: p [4x1] orientation vector
                   %abar [3x1] position in L-RF
            e = p(2:4);
            etil = MatOp.tilde(e);
            %B = [[3x3]*[3x1] , [3x1][1x3] - [3x3]] = [3x4] 
            B = 2*[(p(0)*eye(3) + etil)*aBar , e*aBar' ...
                -(p(0)*eye(3) + etil) *MatOp.tilde(aBar)];
        end
        end
        
    end
    

