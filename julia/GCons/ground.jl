#GCons are Geometic Constraints
include("..\\body.jl") #this needs work !!
type ground
  """
  The absolute ground constraint ties the single supplied bodyi to the ground frame
  """
  sim::Sim      #reference to the simulation data structures
  bi::Body      #bodyi
  Pi::Int       #index of point p on body i (head)
  rDOF::Int     #number of degrees of freedom removed by constraint

  #constructor function
  function ground(sim::Sim,bi::Body,Pi::Int)
    rDOF = 6; #remove all dofs from the body
    new(sim,bi,Pi,rDOF)
  end
end

#----------------begin functions associated with ground----------------------------

function ϕ(con::ground)   #9.26.2016 - slide 20
  """
  constraint equation ϕ
  output: [6 x 1] evaluation of constraint equation
  """
  e = p(con.bi); e = e[2:4,1:1]
  phi = [r(con.bi) ; e ]
end

function ν(con::ground)  #9.26.2016 - slide 21
  """
  RHS of vel equation
  output: [6 x 1] evaluation ν
  """
  nu = zeros(6,1) #this should work as there is no explicit fxn of time
end

function 	γ(con::ground)  #Haug 386
"""
RHS of accel equation
output: [6 x 1] evaluation ν
"""
gamma = zeros(6,1) #all qdot terms are zero, so this should be zero
end

function ϕ_r(con::ground)  #Haug 360
  """
  partial derivative of ϕ WRT position position GC's of both bodyi and bodyj
  output: ([6x3])
  """
  [2*eye(3)*G(p(con.bi)); zeros(3,3)] , 0
end

function ϕ_p(con::ground)  #Haug 360
"""
partial derivative of ϕ WRT position orientation GC's of both bodyi and bodyj
output:([6x4])
"""
 Pi = p(con.bi); e0 = Pi[1] ; e = Pi[2:4,1:1]
 phi_pi = (tilde(e) + e0*eye(3))*G(p(con.bi))
 return phi_pi , 0

end
