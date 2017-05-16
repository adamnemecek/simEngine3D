#GCons are Geometic Constraints
module GCons

type dp1
  """
  The dp1 or dot product 1 constraint reflects the fact that the dot product of
  between a vector on body i and a second vector on body j assume a specified value
  d1 is one of 4 basic GCons that each remove one DOF
  """
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  Pi::Int       #index of point p on body i (tail)
  Qi::Int       #index of point q on body i (head)
  Pj::Int       #index of point P on body j (tail)
  Qj::Int       #index of point Q on body j (head)
  rDOF::Int     #number of degrees of freedom removed by constraint
  f             #lambda function describing how constraint changes with time
  fdot          #lambda function describing how constraint changes with time
  fdot          #lambda function describing how constraint changes with time



  #constructor function
  function dp1(bi::Body,bj::Body,Qi::Int,Qj::Int,Pi=1,Pj=1,f = t->0 , fdot = t->0, fddot = t->0)
    rDOF = 1; #dp1 removes one DOF
    new(bi,bj,Pi,Qi,Pj,Qj,rDOF) #make a new dp1
  end
end

#----------------begin functions associated with dp1----------------------------
#pseudo - getter methods.
aBari(con::dp1) = con.bodyi.pts[:,con.Qi] - con.bodyi.pts[:,con.Pi] #[3x1]
aBarj(con::dp1) = con.bodyj.pts[:,con.Qi] - con.bodyi.pts[:,con.Pi] #[3x1]


function ϕ(con::dp1,sys::system)   #9.26.2016 - slide
  """
  constraint equation ϕ
  output: [1 x 1] evaluation of constraint equation value
  """
  phi = aBari(con)'*A(con.bodyi)'*A(bodyj)*aBarj(con) - con.f(sys.t)
end

function ν(con::dp1,sys::system)  #9.26.2016 - slide 12
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = con.fdot(sys.time)
end

function γ(con::dp1,sys::system)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [1 x 1] evaluation ν
"""
ai = aBari(con) ; aj = aBarj(con)
aiBar = A(con.bodyi)*ai ; ajBar = A(cons.bodyj)*aj
̇pᵢ =body.pdot(con.bodyi) ; pdotj =body.pdot(con.bodyj)

gamma = -ai'B(pdotj,ajBar)*pdotj - aj'B(̇pᵢ,aiBar)*̇pᵢ - fddot(system.time)
end

function ϕ_r(con::dp1,sys::system)
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([1,3],[1,3]) vectors of
  """

  return phi_ri , ph_rj
end

function ϕ_p(r)
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
"""
end

end
