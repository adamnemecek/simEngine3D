%Author: Alex Dawson-Elli
%Assignment: ME 751 HW 7

%reference configuration
r_ref    =[ 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ]; % [3x4] [r1ref; r2ref; r3ref ; r4ref]
drdz_ref =[ 0 0 1 ; 0 0 1 ; 0 0 1 ; 0 0 1 ]; % [3x4] [dr1/dz; dr2/dz; dr3/dz; dr4/dz]

%current configuration
r    = [0 0 0; 1, -.1, 0; 1.1 1.1 .1; 0 1 0];
drdz = [0 0 1; -.1407125, .1407125, .98; -.1407125, -.1407125, .98; 0 0 1];

%element parameters:
E = 0;    
N = sqrt(3/5);  
zeta   = sqrt(3/5); 

%define shape functions:
s1 = 1/4*(1-E)*(1-N); s2 = 1/4*(1+E)*(1-N);
s3 = 1/4*(1+E)*(1+N); s4 = 1/4*(1-E)*(1+N);
Sm = [s1*eye(3), s2*eye(3), s3*eye(3), s4*eye(3)]; % [3x12]


%obtain the Green-Lagrance strain tensor
F = 
