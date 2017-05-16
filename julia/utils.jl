
#this may be moved in the future.
push!(LOAD_PATH , "C:\\Users\\Alex\\Documents\\gitRepos\\simEngine3D")

module Utils
"""
the utilities module provides utility functionality useful throughout
Simbody3D
"""

# function A2P (A::Array{float64})
#   #takes a 3x3 rotation matrix and converts it to a 4x1 array of euler parameters
#   #e0sqrd = (trace(A) + 1)/4
#
#
#   return P
# end

function B(p::Array{float64},sBar::Array{float64}) #9.28.2016 slide 12
  """
  B matrix = [A(p)sbar]_p is useful in calculating partial derivative of ϕ wrt GC's
  inputs: p = euler params [4x1], sBar = point location in LRF [3x1]
  output: B = [3x4]
  """
  e0 = p[1]
  e  = p[2:4]
  b = 2*[(e0*eye(3) + tilde(e))*sBar   e*sBar' - (e0*eye(3) + tilde(e))*tilde(sBar) ]
end

function P2A(P::Array{Float64})   #kinematic key formulas
  """takes a 4x1 array of euler parameters and returns a 3x3 rotation matrix"""
  e0 = p[1]
  e = p[2:4]
  E = [-e , tilde(e) , e0*eye(3)]
  G = [-e , tilde(e) , e0*eye(3)]
  A = EG'
  return A
end

function tilde(a::Array{Float64})  #kinematic key formulas
  """tilde takes a 3x1 vector and makes it the cross product operator matrix ~ """
  Atil = [ 0  -a[1]  a[2];
         a[3]   0   -a[1];
        -a[2]  a[1]    0  ]
  return Atil
end

function insertUL!(A, h , ind)
  """
  insert a smaller matrix h, into a bigger matrix A by overwriting elements of a
  with h. The upper left hand corners of both matricies at ind (x,y)
  input:
    A - larger, outer matrix where a part is being overwritten
    h - smaller matrix to be inserted into a
    ind - (x,y) location
  output:
    !indicates A is motified inplace , nothing is returned
  """
  row, col = size(h)
  x , y = ind[1] , ind[2]
  A[x:x+row-1 ,y:y+col-1] = h
end

#export only needed to have non-dot access
export P2A, tilde, insertUL!


end
