
#this may be moved in the future.
#push!(LOAD_PATH , "C:\\Users\\Alex\\Documents\\gitRepos\\simEngine3D")
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

function B(Pi::Array{Float64},siBar::Array{Float64}) #9.28.2016 slide 12
  """
  B matrix = [A(p)sbar]_p is useful in calculating partial derivative of ϕ wrt GC's
  inputs: Pi = euler params [4x1], siBar = point location in LRF [3x1]
  output: B = [3x4]
  """
  e0 = Pi[1]
  e  = Pi[2:4]
  b = 2*[(e0*eye(3) + tilde(e))*siBar   e*siBar' - (e0*eye(3) + tilde(e))*tilde(siBar) ]
end

function G(Pi::Array{Float64})
  e0 = Pi[1] ; e = Pi[2:4]
  G = [-e , tilde(e) , e0*eye(3)]
end


function P2A(Pi::Array{Float64})   #kinematic key formulas
  """takes a 4x1 array of euler parameters and returns a 3x3 rotation matrix"""
  e0 = Pi[1]
  e = Pi[2:4]
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

function dij(bi::Body,bj::Body,Si::Array,Sj::Array) #9.26.2017 slide 14
  """determines the distance in the GRF between point Pi and Qj"""
  return r(bj) + A(bj)*pt(bj,Sj) -(r(bi) + A(bi)*pt(bi,Si) )
end

function dijdot(bi::Body,bj::Body,Si::Array,Sj::Array) #10.7.2017 slide 7
  """time derivative of dij"""
  return rdot(bj)+ B(p(bj),Sj)*pdot(bj) -rdot(bi) - B(p(bi),Si)*pdot(bi)
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

function flatten(arr::Array)
  """flattens an array of arrays into a array"""
  [arr[1,1]  arr[1,2] ; arr[2,1]  arr[2,2]]
  ## in the future, upgrade this function to work for systems that are not simply 2x2
end
