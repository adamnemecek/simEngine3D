#GCons are Geometic Constraints
type dp2
  """
  The dp2 or dot product 2 constraint reflects the fact that the dot product of
  between a vector ai and a second vector PiQj assume a specified value.
  d2 is one of 4 basic GCons that each remove one DOF
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  Pi::Int       #index of point p on body i (head)
  Qj::Int       #index of point Q on body j (head)
  ai_head::Int  #index of the point at the head of vector ai
  ai_tail::Int  #index of the point at the tial of vector ai
  rDOF::Int     #number of degrees of freedom removed by constraint
  f             #lambda function describing how constraint changes with time
  fdot          #lambda function describing how constraint changes with time
  fddot         #lambda function describing how constraint changes with time


  #constructor function
  function dp2(sim::Sim,bi::Body,bj::Body,Pi::Int,Qj::Int,ai_head::Int, ai_tail = 1, f = t->0 , fdot = t->0, fddot = t->0)
    rDOF = 1; #dp2 removes one DOF
    new(sim,bi,bj,Pi,Qj,ai_head,ai_tail,rDOF,f,fdot,fddot)
  end
end

#----------------begin functions associated with dp1----------------------------
#pseudo - getter methods.
aBari(con::dp2) = pt(con.bodyi,con.ai_head) - pt(con.bodyi,con.ai_tail) #[3x1]
ai(con::dp2) = A(con.bodyi)*aBari(con)
PiQj(con::dp2)  = dij(con.bodyi,con.bodyj,pt(con.bodyi,con.Pi),pt(con.bodyj,con.Qj))
siBar(con::dp2) = pt(con.bodyi,con.Pi)
sjBar(con::dp2) = pt(con.bodyj,con.Qj)

function ϕ(con::dp2)   #9.26.2016 - slide 14
  """
  constraint equation ϕ
  output: [1 x 1] evaluation of constraint equation value
  """
  phi = aBari(con)'*A(con.bodyi)'PiQj(con) - con.f(con.sim.t)
end

function ν(con::dp2)  #9.26.2016 - slide 15
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = con.fdot(con.sim.t)
end

function 	𝛾(con::dp2)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [1 x 1] evaluation ν
"""
aidot =B(p(con.bodyi) , aBari(con))*pdot(con.bodyi)
pdoti = pdot(con.bodyi) ; pdotj = pdot(con.bodyj)
d_ijdot = dijdot(con.bodyi,con.bodyj,siBar(con),siBar(con))

gamma = -ai(con)'*B(pdotj,sjBar(con))*pdotj +
         ai(con)'*B(pdoti,siBar(con))*pdoti -
         PiQj(con)'*B(pdoti,aBari(con))*pdoti -
         2*aidot'*d_ijdot + con.fddot(con.sim.t)
end

function ϕ_r(con::dp2)  #9.28.2016 slide 15
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([1x3],[1x3])
  """
  phi_ri = -aBari(con)'*A(con.bodyi)'
  phi_rj = -phi_ri;
  return phi_ri , phi_rj
end

function ϕ_p(con::dp2)  # #9.28.2016 slide 15
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
output:([1x4],[1x4])
"""
Pj = p(con.bodyj) ; Pi = p(con.bodyi)

phi_pi = PiQj(con)'*B(Pi,aBari(con)) - ai(con)'*B(Pi,siBar(con))
phi_pj = ai(con)'*B(Pj,sjBar(con))

return phi_pi , phi_pj
end


#---------------------ɸλ_qq values for ψFull------------------------------------
#ϕ_rr
ϕ_riri(con::dp2) = zeros(3,3)
ϕ_rirj(con::dp2) = zeros(3,3)
ϕ_rjrj(con::dp2) = zeros(3,3)

#ϕ_rp
ϕ_ripi(con::dp2) = -B(p(con.bodyi),aBari(con))
ϕ_ripj(con::dp2) = zeros(3,4)
ϕ_rjpi(con::dp2) = -1*ϕ_ripi(con)
ϕ_rjpj(con::dp2) = zeros(3,4)

#ϕ_pp
ϕ_pipi(con::dp2) =  K(aBari(con),PiQj(con)) - K(siBar(con), A(con.bodyi)*aBari(con)) -
                    B(p(con.bodyi),aBari(con))'*B(p(con.bodyi),siBar(con)) -
                    B(p(con.bodyi),siBar(con))'*B(p(con.bodyi),aBari(con))
ϕ_pipj(con::dp2) =  B(p(con.bodyi),aBari(con))'*B(p(con.bodyj),sjBar(con))
ϕ_pjpj(con::dp2) =  K(sjBar(con),A(con.bodyi)*aBari(con))
