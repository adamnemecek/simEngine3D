#GCons are Geometic Constraints
type cd
  """
  The cd or coordinate difference constraint requires that the absolute distance between a point
  Pi (on body i) and a second point Qj (on body j) maintain a particular magnitude greater than zero
  cd is one of 4 basic GCons that each remove one DOF
  """
  sim::Sim      #reference to the simulation data structures
  bi::Body      #bodyi
  bj::Body      #bodyj
  Pi::Int       #index of point p on body i (head)
  Qj::Int       #index of point Q on body j (head)
  c::Array      #column vect specifying coord of interest [1 0 0]' for x
  rDOF::Int     #number of degrees of freedom removed by constraint
  f             #lambda function describing how constraint changes with time
  fdot          #lambda function describing how constraint changes with time
  fddot         #lambda function describing how constraint changes with time


  #constructor function
  function cd(sim::Sim,bi::Body,bj::Body,Pi::Int,Qj::Int,c::Array, f = t->0 , fdot = t->0, fddot = t->0)
    rDOF = 1; #cd 1 removes one DOF
    new(sim,bi,bj,Pi,Qj,c,rDOF,f,fdot,fddot)
  end
end

#----------------begin functions associated with dp1----------------------------
#pseudo - getter methods.
PiQj(con::cd)  = dij(con.bi,cons.bj,pt(con.bi,con.Pi),pt(con.bj,con.Qj))

function ϕ(con::cd)   #9.26.2016 - slide 20
  """
  constraint equation ϕ
  output: [1 x 1] evaluation of constraint equation value
  """
  phi = con.c'*PiQj(con) - con.f(con.sim.t)
end

function ν(con::cd)  #9.26.2016 - slide 21
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = con.fdot(con.sim.t)
end

function 	𝛾(con::cd)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [1 x 1] evaluation 𝛾
"""
SiBar = pt(con.bi,con.Pi) ; SjBar = pt(con.bj,con.Qj)
pdoti = pdot(con.bi) ; pdotj = pdot(con.bj)

gamma = con.c'*B(pdoti,siBar)*pdoti - con.c'*B(pdotj,sjBar)*pdotj + con.fddot(con.sim.t)
end

function ϕ_r(con::cd)  #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([1x3],[1x3])
  """
  return -con.c' , con.c'
end

function ϕ_p(con::cd)  # #9.28.2016 slide 17
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
output:([1x4],[1x4])
"""
SjBar = pt(con.bodyj,con.Qj) ; SiBar = pt(con.bodyi,con.Pi)
Pj = p(con.bodyj) ; Pi = p(con.bodyi)

phi_pi = -con.c'*B(Pi,SiBar)
phi_pj =  con.c'*B(Pj,SjBar)

return phi_pi , phi_pj
end
