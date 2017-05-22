#GCons are Geometic Constraints
type uj
  """
  The uj or Universal Joint combines a spherical joint with a dp1 to for for a
  total rDOF of 4
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  Pi::Int       #index of head of vector Pi
  Qj::Int       #index of head of vector Qj
  ai_head::Int  #index of ai_head
  ai_tail::Int  #index of a1_tail
  aj_head::Int  #index of ai_head
  aj_tail::Int  #index of aj_tail
  rDOF::Int     #number of degrees of freedom removed by constraint
  subGCs::Array #consider using inhepitance here for speedup?

  #constructor function
  function uj(sim::Sim,bodyi::Body,bodyj::Body,Pi,Qj,ai_head,aj_head,ai_tail = 1,aj_tail = 1)
    rDOF = 4;

    subGCs = Array(Any,2)
    subGCs[1] =  sj(sim,bodyi,bodyj,Pi,Qj)
    subGCs[2] = dp1(sim,bodyi,bodyj,ai_head,aj_head,a,ai_tail,aj_tail)

    new(sim,bodyi,bodyj,Pi,Qj,ai_head,ai_tail,bi_head,bi_tail,cj_head,cj_tail,rDOF,subGCs)
  end
end

#----------------begin functions associated with sj-----------------------------

function ϕ(con::uj)   #9.26.2016 - slide 24
  """
  constraint equation ϕ
  output: [4 x 1] evaluation of constraint equation value
  """
  phi = [ ϕ(con.subGCs[1]) ; ϕ(con.subGCs[2]) ]
end

function ν(con::uj)
  """
  RHS of vel equation
  output: [4 x 1] evaluation ν
  """
  nu = [ ν(con.subGCs[1]) ; ν(con.subGCs[2]) ]
end

function 	𝛾(con::uj)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [4 x 1] evaluation ν
"""
gamma = [ 𝛾(con.subGCs[1]) ; 𝛾(con.subGCs[2])]
end

function ϕ_r(con::uj)  #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([4x3],[4x3])
  """
  phi_r = Array(Array{Float64},2,2)
  phi_r[1,1], phi_r[1,2] = ϕ_r(con.subGCs[1])
  phi_r[2,1], phi_r[2,2] = ϕ_r(con.subGCs[2])
  phi_r = flatten(phi_r)
  return phi_r[:,1] , phi_r[:,2]
end

function ϕ_p(con::uj)  # #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
  output:([4x4],[4x4])
  """
  phi_p = Array(Array{Float64},2,2)
  phi_p[1,1], phi_p[1,2] = ϕ_p(con.subGCs[1])
  phi_p[2,1], phi_p[2,2] = ϕ_p(con.subGCs[2])
  phi_p = flatten(phi_p)
  return phi_p[:,1] , phi_p[:,2]
end
