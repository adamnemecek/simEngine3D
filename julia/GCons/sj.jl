#GCons are Geometic Constraints
include("..\\body.jl") #this needs work !!
type sj
  """
  The sj constraint specifies that two points, Pi and Qj remain coincident
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  Pi::Int       #index of head of vector Pi
  Qj::Int       #index of head of vector Qj
  rDOF::Int     #number of degrees of freedom removed by constraint
  subGCs::Array #consider using inhepitance here for speedup?


  #constructor function
  function sj(sim::Sim,bodyi::Body,bodyj::Body,Pi,Qj)
    rDOF = 3;

    subGCs = Array(Any,3)
    subGCs[1] = cd(sim,bodyi,bodyj,Pi,Qj,[1 0 0])
    subGCs[2] = cd(sim,bodyi,bodyj,Pi,Qj,[0 1 0])
    subGCs[3] = cd(sim,bodyi,bodyj,Pi,Qj,[0 0 1])

    new(sim,bodyi,bodyj,Pi,Qj, rDOF,subGCs)
  end
end

#----------------begin functions associated with sj----------------------------

function ϕ(con::sj)   #9.26.2016 - slide 24
  """
  constraint equation ϕ
  output: [3 x 1] evaluation of constraint equation value
  """
  phi = [ ϕ(subGCs[1]) ; ϕ(subGCs[2]) ; ϕ(subGCs[3])]
end

function ν(con::sj)
  """
  RHS of vel equation
  output: [3 x 1] evaluation ν
  """
  nu = [ ν(subGCs[1]) ; ν(subGCs[2]);  ν(subGCs[3])]
end

function 	γ(con::b2)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [3 x 1] evaluation ν
"""
gamma = [ γ(subGCs[1]) ; γ(subGCs[2]); γ(subGCs[3])]
end

function ϕ_r(con::sj)  #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([3x3],[3x3])
  """
  phi_r = Array(Float64,con.rDOF,2)
  phi_r[1,1], phi_r[1,2] = ϕ_r(subGCs[1])
  phi_r[2,1], phi_r[2,2] = ϕ_r(subGCs[2])
  phi_r[3,1], phi_r[3,2] = ϕ_r(subGCs[3])
  return phi_r[:,1] , phi_r[:,2]
end

function ϕ_p(con::sj)  # #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
  output:([3x4],[3x4])
  """
  phi_p = Array(Float64,con.rDOF,2)
  phi_p[1,1], phi_p[1,2] = ϕ_p(subGCs[1])
  phi_p[2,1], phi_p[2,2] = ϕ_p(subGCs[2])
  phi_p[3,1], phi_p[3,2] = ϕ_p(subGCs[3])
  return phi_p[:,1] , phi_p[:,2]
end
