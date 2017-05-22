#GCons are Geometic Constraints
type ep
  """
  The p constraint is the euler normalization constraint, which specifies that
  p'p = 1
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body   #bodyi
  rDOF::Int     #number of degrees of freedom removed by constraint


  #constructor function
  function ep(sim::Sim,bodyi::Body)
    rDOF = 1;
    new(sim,bodyi,rDOF)
  end
end

#----------------begin functions associated with b1----------------------------

function ϕ(con::ep)   #9.26.2016 - slide 23
  """
  constraint equation ϕ
  output: [1 x 1] evaluation of constraint equation value
  """
  phi = .5*p(con.bodyi)'*p(con.bodyi) - .5
end

function ν(con::ep)
  """
  RHS of vel equation
  output: [1 x 1] evaluation ν
  """
  nu = 0
end

function 	𝛾(con::ep)
  """
  RHS of accel equation
  output: [1 x 1] evaluation ν
  """
  gamma = -2*pdot(con.bodyi)'*pdot(con.bodyi)
end

function ϕ_r(con::ep)
  """
  partial derivative of ϕ WRT position position GC's of bodyi
  output: ([1x3])
  """
  return zeros(1,3)
end

function ϕ_p(con::ep)
  """
  partial derivative of ϕ WRT position orientation GC's of bodyi
  output:([1x4])
  """
  return phi_p = 2*p(con.bodyi)'

end
