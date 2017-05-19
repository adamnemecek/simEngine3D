#GCons are Geometic Constraints
type cj
  """
  The cj or cylindrical joint  allows for translation along an axis, and rotation about
  the same axis, it is constructed using 2 intermediate constraints.
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  bodyj::Body   #bodyj
  Pi::Int       #index of head of vector Pi
  Qj::Int       #index of head of vector Qj
  ai_head::Int  #index of ai_head
  ai_tail::Int  #index of a1_tail
  bi_head::Int  #index of bi_head
  bi_tail::Int  #index of bi_tail
  cj_head::Int  #index of cj_head
  cj_tail::Int  #index of cj_tail
  rDOF::Int     #number of degrees of freedom removed by constraint
  subGCs::Array #consider using inhepitance here for speedup?

  #constructor function
  function cj(sim::Sim,bodyi::Body,bodyj::Body,Pi,Qj,ai_head,bi_head,cj_head,ai_tail = 1,bi_tail = 1,cj_tail = 1)
    rDOF = 5;

    subGCs = Array(Any,2)
    subGCs[1] = b1(sim,bodyi,bodyj,ai_head,bi_head,cj_head,ai_tail,bi_tail,cj_tail)
    subGCs[2] = b2(sim,bodyi,bodyj,Pi,Qj,ai_head,bi_head,ai_tail,bi_tail)

    new(sim,bodyi,bodyj,Pi,Qj,ai_head,ai_tail,bi_head,bi_tail,cj_head,cj_tail,rDOF,subGCs)
  end
end

#----------------begin functions associated with sj-----------------------------

function ϕ(con::cj)   #9.26.2016 - slide 24
  """
  constraint equation ϕ
  output: [5 x 1] evaluation of constraint equation value
  """
  phi = [ ϕ(subGCs[1]) ; ϕ(subGCs[2]) ]
end

function ν(con::cj)
  """
  RHS of vel equation
  output: [5 x 1] evaluation ν
  """
  nu = [ ν(subGCs[1]) ; ν(subGCs[2]) ]
end

function 	𝛾(con::cj)  #10.7.2016 - slide 8
"""
RHS of accel equation
output: [5 x 1] evaluation 𝛾
"""
gamma = [ 𝛾(subGCs[1]) ; 𝛾(subGCs[2])]
end

function ϕ_r(con::cj)  #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([5x3],[5x3])
  """
  phi_r = Array(Array{Float64},2,2)
  phi_r[1,1], phi_r[1,2] = ϕ_r(subGCs[1])
  phi_r[2,1], phi_r[2,2] = ϕ_r(subGCs[2])
  phi_r = flatten(phi_r)
  return phi_r[:,1] , phi_r[:,2]
end

function ϕ_p(con::cj)  # #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
  output:([5x4],[5x4])
  """
  phi_p = Array(Array{Float64},2,2)
  phi_p[1,1], phi_p[1,2] = ϕ_p(subGCs[1])
  phi_p[2,1], phi_p[2,2] = ϕ_p(subGCs[2])
  phi_p = flatten(phi_p)
  return phi_p[:,1] , phi_p[:,2]
end
