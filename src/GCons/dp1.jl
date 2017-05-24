#GCons are Geometic Constraints
type dp1
  """
  The dp1 or dot product 1 constraint reflects the fact that the dot product of
  between a vector on body i and a second vector on body j assume a specified value
  d1 is one of 4 basic GCons that each remove one DOF
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  Pi::Int       #index of point P on body i (tail)
  Qi::Int       #index of point Q on body i (head)
  Pj::Int       #index of point P on body j (tail)
  Qj::Int       #index of point Q on body j (head)
  rDOF::Int     #number of degrees of freedom removed by constraint
  f             #lambda function describing how constraint changes with time
  fdot          #lambda function describing how constraint changes with time
  fddot         #lambda function describing how constraint changes with time



  #constructor function
  function dp1(sim::Sim,bi::Body,bj::Body,Qi::Int,Qj::Int,Pi=1,Pj=1,f = t->0 , fdot = t->0, fddot = t->0)
    rDOF = 1; #dp1 removes one DOF
    new(sim,bi,bj,Pi,Qi,Pj,Qj,rDOF,f,fdot,fddot)
  end
end

#----------------begin functions associated with dp1----------------------------
#pseudo - getter methods.
aBari(con::dp1) = pt(con.bodyi,con.Qi) - pt(con.bodyi,con.Pi) #[3x1]
aBarj(con::dp1) = pt(con.bodyi,con.Qj) - pt(con.bodyi,con.Pj) #[3x1]


function ϕ(con::dp1)   #9.26.2016 - slide 11
  """
  constraint equation ϕ
  output: [1 x 1] evaluation of constraint equation value
  """
  phi = aBari(con)'*A(con.bodyi)'*A(con.bodyj)*aBarj(con) - con.f(con.sim.t)
end

function ν(con::dp1)  #9.26.2016 - slide 12
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = con.fdot(con.sim.t)
end

function 	𝛾(con::dp1)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [1 x 1] evaluation ν
"""
aiBar = aBari(con) ; ajBar = aBarj(con)
ai = A(con.bodyi)*aiBar ; aj = A(con.bodyj)*ajBar
aidot =B(p(con.bodyi) ,aiBar)*pdot(con.bodyi)
ajdot =B(p(con.bodyj) ,ajBar)*pdot(con.bodyj)
pdoti = pdot(con.bodyi) ; pdotj = pdot(con.bodyj)

gamma = -ai'B(pdotj,ajBar)*pdotj - aj'B(pdoti,aiBar)*pdoti - 2*aidot'ajdot + con.fddot(con.sim.t)
end

function ϕ_r(con::dp1)  #9.28.2016 slide 13
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([1x3],[1x3])
  """
  phi_ri = zeros(1,3) ; phi_rj = zeros(1,3)
  return phi_ri , phi_rj
end

function ϕ_p(con::dp1)  # #9.28.2016 slide 13
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
output:([1x4],[1x4])
"""
phi_pi = aBarj(con)'*A(con.bodyj)'*B(p(con.bodyi),aBari(con))
phi_pj = aBari(con)'*A(con.bodyi)'*B(p(con.bodyj),aBarj(con))
return phi_pi , phi_pj
end
