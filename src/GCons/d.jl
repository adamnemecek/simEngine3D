#GCons are Geometic Constraints
type d
  """
  The d or distance constraint requires that the absolute distance between a point
  Pi (on body i) and a second point Qj (on body j) maintain a particular magnitude greater than zero
  d is one of 4 basic GCons that each remove one DOF
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body      #bodyi
  bodyj::Body      #bodyj
  Pi::Int       #index of point p on body i (head)
  Qj::Int       #index of point Q on body j (head)
  rDOF::Int     #number of degrees of freedom removed by constraint
  f             #lambda function describing how constraint changes with time
  fdot          #lambda function describing how constraint changes with time
  fddot         #lambda function describing how constraint changes with time


  #constructor function
  function d(sim::Sim,bodyi::Body,bodyj::Body,Pi::Int,Qj::Int, f = t->0 , fdot = t->0, fddot = t->0)
    rDOF = 1; #d removes one DOF
    new(sim,bodyi,bodyj,Pi,Qj,rDOF,f,fdot,fddot)
  end
end

#----------------begin functions associated with dp1----------------------------
#pseudo - getter methods.
PiQj(con::d)  = dij(con.bodyi,cons.bodyj,pt(con.bodyi,con.Pi),pt(con.bodyj,con.Qj))

function ϕ(con::d)   #9.26.2016 - slide 17
  """
  constraint equation ϕ
  output: [1 x 1] evaluation of constraint equation value
  """
  phi = PiQj(con)'PiQj(con) - con.f(con.sim.t)
end

function ν(con::d)  #9.26.2016 - slide 18
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = con.fdot(con.sim.t)
end

function 	𝛾(con::d)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [1 x 1] evaluation 𝛾
"""
SiBar = pt(con.bodyi,con.Pi) ; SjBar = pt(con.bodyj,con.Qj)
pdoti = pdot(con.bodyi) ; pdotj = pdot(con.bodyj)
d_ijdot = dijdot(con.bodyi,con.bodyj,Si,Sj)

gamma = -2*PiQj(con)'B(pdotj,SjBar)*pdotj + 2*PiQj(con)'B(pdoti,SiBar)pdoti - 2*d_ijdot'd_ijdot + con.fddot(con.sim.t)
end

function ϕ_r(con::d)  #9.28.2016 slide 15
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([1x3],[1x3])
  """
  phi_ri = -2*PiQj(con)'
  phi_rj = -phi_ri
  return phi_ri,phi_rj
end

function ϕ_p(con::d)  # #9.28.2016 slide 15
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
output:([1x4],[1x4])
"""
SjBar = pt(con.bodyj,con.Qj) ; SiBar = pt(con.bodyi,con.Pi)
Pj = p(con.bodyj) ; Pi = p(con.bodyi)

phi_pi = -2*PiQj(con)'*B(Pi,SiBar)
phi_pj =  2*PiQj(con)'*B(Pj,SjBar)

return phi_pi , phi_pj
end
