"""
the body class contains information about the mass properties, location and orientation
and critical points used to form constraints and forces acting on the body.
"""
type Body
  sim::Sim              #This is a reference / pointer to sim level variables (same memory location for each body)
  ID::Int               #unique identifier of the order in which bodies where added
  pts::Array            #a [3 x n] array where n is the number of points. the first column is always zeros

  #dynamics
  m::Real              #mass of the body
  J::Array             #inertia tensor

  #constructor function
  function Body(sim::Sim, ID::Int,m = 1, J = eye(3))
    pts = [0.0 0.0 0.0]'
    new(sim,ID,pts,m,J)
  end
end

#---------------------index accesser methods-----------------------------------
"""returns a range object that specifies where the coordinates of a body are located in the system matrix"""
rr(body::Body) = 3*(body.ID-1)+1:3*body.ID
pr(body::Body) = 3*body.sim.nb+4(body.ID-1)+1:3*body.sim.nb+4*body.ID

#-----------------pseudo getter methods----------------------------------------
"""retrieve a point from the point body point matrix pts"""
pt(bd::Body,ind::Int) =  bd.pts[:,ind:ind]

"""get translations kinematics quantities"""
r(bd::Body) = bd.sim.q[rr(bd) , 1:1]
rdot(bd::Body) =  bd.sim.qdot[rr(bd) , 1:1]
rddot(bd::Body) = bd.sim.qddot[rr(bd), 1:1]

"""get orientation kinematics quantities"""
p(bd::Body) =  bd.sim.q[pr(bd) , 1:1]
pdot(bd::Body) = bd.sim.qdot[pr(bd), 1:1]
pddot(bd::Body) = bd.sim.qddot[pr(bd) , 1:1]

"""returns rotation matrix of the body, calculated from p"""
A(bd::Body) = P2A(p(bd))

"""returns G matrix of a body, calculated from p"""
G(bd::Body) =  G(p(bd))

"""returns G matrix of a body, calculated from p"""
Gdot(bd::Body) = G(pdot(bd))

"""take a 3x1 ωbar and a body and returns global pdot (for setting up IC's)"""
ωbar2pdot(bd::Body, ωbar::Array) = .5*G(bd)'*ωbar

"""returns G matrix of a body, calculated from p"""
E(bd::Body) =  E(p(bd))


#----------------------pseudo setter methods------------------------------------
"""set translations kinematics quantities"""
function set_r!(bd::Body, r::Array)
   bd.sim.q[rr(bd) , 1:1] = r
end
function set_rdot!(bd::Body, rdot::Array)
 bd.sim.qdot[rr(bd) , 1:1] = rdot
end
function set_rddot!(bd::Body , rddot::Array)
 bd.sim.qddot[rr(bd), 1:1] = rddot
end

"""set translations kinematics quantities"""
function set_p!(bd::Body , p::Array)
 bd.sim.q[pr(bd) , 1:1] = p
end
function set_pdot!(bd::Body, pdot::Array)
 bd.sim.qdot[pr(bd) , 1:1] = pdot
end
function set_pddot!(bd::Body , pddot::Array)
 bd.sim.qddot[pr(bd), 1:1] = pddot
end

"""adds a new point to a bodies' point matrix"""
function addPoint(bd::Body, point::Array)
 bd.pts = [bd.pts point]
end
