#GCons are Geometic Constraints
type b2
  """
  The b2 constraint specifies that a vector PiQj is perpendicular to a plane
  defined by two vectors on bodyi , ai and bi
  b1 is one of two intermediate constraints
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
  rDOF::Int     #number of degrees of freedom removed by constraint
  subGCs::Array #consider using inhepitance here for speedup?


  #constructor function
  function b2(sim::Sim,bodyi::Body,bodyj::Body,Pi,Qj, ai_head, bi_head, ai_tail = 1, bi_tail = 1)
    rDOF = 2;

    subGCs = Array(Any,2)
    subGCs[1] = dp2(sim,bodyi,bodyj,Pi,Qj,ai_head,ai_tail)
    subGCs[2] = dp2(sim,bodyi,bodyj,Pi,Qj,bi_head,bi_tail)

    new(sim,bodyi,bodyj,Pi,Qj,ai_head,ai_tail,bi_head,bi_tail,rDOF,subGCs)
  end
end

#----------------begin functions associated with b2----------------------------

function ϕ(con::b2)   #9.26.2016 - slide 24
  """
  constraint equation ϕ
  output: [2 x 1] evaluation of constraint equation value
  """
  phi = [ ϕ(subGCs[1]) ; ϕ(subGCs[2])]
end

function ν(con::b2)
  """
  RHS of vel equation
  output: [2 x 1] evaluation ν
  """
  nu = [ ν(subGCs[1]) ; ν(subGCs[2])]
end

function 	𝛾(con::b2)  #10.7.2016 - slide 8
  """
  RHS of accel equation
  output: [2 x 1] evaluation ν
  """
  gamma = [ 𝛾(subGCs[1]) ; 𝛾(subGCs[2])]
end

function ϕ_r(con::b2)  #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([2x3],[2x3])
  """
  phi_r = Array(Float64,con.rDOF,2)
  phi_r[1,1], phi_r[1,2] = ϕ_r(subGCs[1])
  phi_r[2,1], phi_r[2,2] = ϕ_r(subGCs[2])
  return phi_r[:,1] , phi_r[:,2]
end

function ϕ_p(con::b2)  # #9.28.2016 slide 17
  """
  partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
  output:([2x4],[2x4])
  """
  phi_p = Array(Float64,con.rDOF,2)
  phi_p[1,1], phi_p[1,2] = ϕ_p(subGCs[1])
  phi_p[2,1], phi_p[2,2] = ϕ_p(subGCs[2])
  return phi_p[:,1] , phi_p[:,2]
end
