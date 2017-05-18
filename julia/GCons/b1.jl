#GCons are Geometic Constraints
include("..\\body.jl") #this needs work !!
type b1
  """
  The b1 constraint specifies that a vector cj is perpendicular to a plane
  defined by two vectors on bodyi , ai and bi
  b1 is one of two intermediate constraints
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  ai_head::Int  #index of ai_head
  ai_tail::Int  #index of a1_tail
  bi_head::Int  #index of bi_head
  bi_tail::Int  #index of bi_tail
  cj_head::Int  #index of cj_head
  cj_tail::Int  #index of cj_tail
  rDOF::Int     #number of degrees of freedom removed by constraint
  subGCs::Array #consider using inhepitance here for speedup?


  #constructor function
  function b1(sim::Sim,bodyi::Body,bodyj::Body,ai_head,bi_head,cj_head, ai_tail = 1, bi_tail = 1,cj_tail = 1)
    rDOF = 2; #cd 1 removes one DOF

    subGCs = Array(Any,rDOF)
    subGCs[1] = dp1(sim,bodyi,bodyj,ai_head,cj_head,ai_tail,cj_tail)
    subGCs[2] = dp1(sim,bodyi,bodyj,bi_head,cj_head,bi_tail,cj_tail)

    new(sim,bodyi,bodyj,ai_head,ai_tail,bi_head,bi_tail,cj_head,cj_tail,rDOF,subGCs)
  end
end

#----------------begin functions associated with b1----------------------------

function ϕ(con::b1)   #9.26.2016 - slide 20
  """
  constraint equation ϕ
  output: [2 x 1] evaluation of constraint equation value
  """
  phi = [ ϕ(subGCs[1]) ; ϕ(subGCs[2])]
end

function ν(con::b1)  #9.26.2016 - slide 21
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = [ ν(subGCs[1]) ; ν(subGCs[2])]
end

function 	γ(con::b1)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [2 x 1] evaluation ν
"""
gamma = [ γ(subGCs[1]) ; γ(subGCs[2])]
end

function ϕ_r(con::b1)  #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([2x3],[2x3])
  """
  phi_r = Array(Float64,2,2)
  phi_r[1,1], phi_r[1,2] = ϕ_r(subGCs[1])
  phi_r[2,1], phi_r[2,2] = ϕ_r(subGCs[2])
  return phi_r[:,1] , phi_r[:,2]
end

function ϕ_p(con::b1)  # #9.28.2016 slide 17
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
output:([2x4],[2x4])
"""
phi_pi_1, phi_pj_1 = ϕ_p(subGCs[1])
phi_pi_2, phi_pj_2 = ϕ_p(subGCs[2])
phi_pi = [phi_pi_1 ; phi_pi_2] ; phi_pj = [phi_pj_1 ; phi_pj_2]
return phi_pi , phi_pj
end
