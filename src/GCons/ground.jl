#GCons are Geometic Constraints
type ground
  """
  The absolute ground constraint ties the single supplied bodyi to the ground frame
  """
  sim::Sim      #reference to the simulation data structures
  bodyi::Body      #bodyi
  Pi::Int       #index of point p on body i (head)
  rDOF::Int     #number of degrees of freedom removed by constraint

  #constructor function
  function ground(sim::Sim,bodyi::Body,Pi::Int)
    rDOF = 6; #remove all dofs from the body
    new(sim,bodyi,Pi,rDOF)
  end
end

#----------------begin functions associated with ground----------------------------

function ϕ(con::ground)   #Haug 358
  """
  constraint equation ϕ
  output: [6 x 1] evaluation of constraint equation
  """
  e = p(con.bodyi); e = e[2:4,1:1]
  phi = [r(con.bodyi) ; e ]
end

function ν(con::ground)  #
  """
  RHS of vel equation
  output: [6 x 1] evaluation ν
  """
  nu = zeros(6,1) #this should work as there is no explicit fxn of time
end

function 	𝛾(con::ground)  #Haug 386
"""
RHS of accel equation
output: [6 x 1] evaluation ν
"""
gamma = zeros(6,1) #all qdot terms are zero, so this should be zero
end

function ϕ_r(con::ground)  #Haug 360
  """
  partial derivative of ϕ WRT position position GC's bodyi
  output: ([6x3])
  """
  [eye(3); zeros(3,3)]
end

function ϕ_p(con::ground)  #Haug 360
"""
partial derivative of ϕ WRT position orientation GC's bodyi
output:([6x4])
"""
P_i = p(con.bodyi); e0 = P_i[1] ; e = P_i[2:4,1:1] ; siBar = pt(con.bodyi,con.Pi)
phi_pi =[-2*A(con.bodyi)*tilde(siBar)'*G(P_i) ; (tilde(e) + e0*eye(3))*G(P_i)]
return phi_pi
end
