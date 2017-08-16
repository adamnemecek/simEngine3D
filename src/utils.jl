
#the utilities module provides utility functionality useful throughout
#Simbody3D
"""
B matrix = [A(p)sbar] p is useful in calculating partial derivative of ϕ wrt GC's
inputs: Pi = euler params [4x1], siBar = point location in LRF [3x1]
output: B = [3x4]
"""

function B(Pi::Array,siBar::Array) #9.28.2016 slide 12

  e0 = Pi[1]
  e  = Pi[2:4]
  b = 2*[(e0*eye(3) + tilde(e))*siBar   e*siBar' - (e0*eye(3) + tilde(e))*tilde(siBar) ]
end

"""[3x4] matrix used for converting between EP and ω orientation representations"""
function G(p::Array)
  e0 = p[1] ; e = p[2:4]
  G = [-e -tilde(e) + e0*eye(3)]
end

"""[3x4] matrix used for converting between EP and ω orientation representations"""
function E(p::Array)
  e0 = p[1] ; e = p[2:4]
  E = [-e tilde(e) + e0*eye(3)]
end

""" K matrix is useful for calculating derivatives of ϕλ_(rr,rp,pr,pp) 10.19 slide 9"""
function K(abar, b)
  k = 2*[abar'*b          abar'*tilde(b);
         tilde(abar)*b    abar*b'+ b*abar'- (abar'*b)[1]*eye(3)];

end



"""takes a 4x1 array of euler parameters and returns a 3x3 rotation matrix"""
P2A(p::Array) = E(p)*(G(p)')

"""orientational representation converters"""
pdot2ωbar(p::Array, pdot::Array) = 2*G(p)*pdot
pddot2ωdotbar(p::Array, pddot::Array) = 2*G(p)*pddot

"""takes a 3x3 rotation matrix and converts it to a 4x1 array of euler parameters"""
function A2P(A::Array) #9.21 slide 20
  e0 = sqrt((trace(A) + 1)/4)
  if e0 != 0
    e1 = (A[3,2] - A[2,3])/(4*e0)
    e2 = (A[1,3] - A[3,1])/(4*e0)
    e3 = (A[2,1] - A[1,2])/(4*e0)
  end

  if e0 == 0  #implies Χ = π
    #figure out which e terms are non-zero
    e1flag = false; e2flag = false; e3flag = false;
    if A[1,1] + 1 != 0  e1flag = true end
    if A[2,2] + 1 != 0  e2flag = true end
    if A[3,3] + 1 != 0  e3flag = true end

    if e1flag
      e1 = sqrt((A[1,1] + 1)/2)
      e2 = (A[2,1] + A[1,2])/(4*e1)
      e3 = (A[3,1] + A[1,3])/(4*e1)
    elseif e2flag
      e2 = sqrt((A[2,2] + 1)/2)
      e1 = (A[2,1] + A[1,2])/(4*e2)
      e3 = (A[3,2] + A[2,3])/(4*e2)
    elseif e3flag
      e3 = sqrt((A[3,3] + 1)/2)
      e1 = (A[3,1] + A[1,3])/(4*e3)
      e2 = (A[3,2] + A[2,3])/(4*e3)
    else
      warn("something is wrong with A value  $A ")
    end
  end

  p = [e0 e1 e2 e3]'
end

"""tilde takes a 3x1 vector and makes it the cross product operator matrix ~ """
function tilde(a::Array)  #kinematic key formulas
  Atil = [ 0  -a[3]  a[2];
         a[3]    0  -a[1];
        -a[2]  a[1]    0  ]
  return Atil
end

"""determines the distance in the GRF between point Pi and Qj"""
function dij(bi::Body,bj::Body,PiBar::Array,QjBar::Array) #9.26.2017 slide 14
  return r(bj) + A(bj)*QjBar -(r(bi) + A(bi)*PiBar )

end

"""time derivative of dij"""
function dijdot(bi::Body,bj::Body,PiBar::Array,QjBar::Array) #10.7.2017 slide 7
  return rdot(bj)+ B(p(bj),QjBar)*pdot(bj) -rdot(bi) - B(p(bi),PiBar)*pdot(bi)
end


"""
insert a smaller matrix h, into a bigger matrix A by overwriting elements of A
with h. The upper left hand corners of both matricies at ind (x,y)
input:
  A - larger, outer matrix where a part is being overwritten
  h - smaller matrix to be inserted into a
  ind - (x,y) location
output:
  !indicates A is motified inplace , nothing is returned
"""
function insertUL!(A, h , ind)
  row, col = size(h)
  x , y = ind[1] , ind[2]
  A[x:x+row-1 ,y:y+col-1] = h
end

"""sum each element of a row and return the resulting column vector"""
function rowSum(arr::Array)
  result = zeros(size(arr)[1],1) #column vector
  for row in 1:size(arr)[1]
    for col in 1:size(arr)[2]
      result[row,1] += arr[row,col]
    end
  end
  return result
end

#---------------------setup functions-------------------------------------------
"""principle rotations are useful when setting up problems"""
Rx(Θ) = [1 0 0 ; 0 cos(Θ) -sin(Θ) ; 0 sin(Θ) cos(Θ)]
Ry(Θ) = [ cos(Θ) 0 sin(Θ) ; 0 1 0 ; -sin(Θ) 0 cos(Θ)]
Rz(Θ) = [ cos(Θ) -sin(Θ) 0 ; sin(Θ) cos(Θ) 0  ; 0 0 1]
