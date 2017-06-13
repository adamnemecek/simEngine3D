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

#-----------------pseudo getter methods----------------------------------------
function pt(bd::Body,ind::Int)
  """retrieve a point from the point body point matrix pts"""
  return bd.pts[:,ind:ind]
end

function r(bd::Body)
  """return r of body by accessing sim q"""
  return bd.sim.q[3*(bd.ID-1)+1:3*bd.ID , 1:1]
end
function rdot(bd::Body)
  """return rdot of body by accessing sim q"""
    return bd.sim.qdot[3*(bd.ID-1)+1:3*bd.ID , 1:1]
end
function rddot(bd::Body)
  """return rddot of body by accessing sim q"""
  return bd.sim.qddot[3*(bd.ID-1)+1:3*bd.ID , 1:1]
end

function p(bd::Body)
  """return p of body by accessing sim q"""
  return bd.sim.q[3*bd.sim.nb+4(bd.ID-1)+1:3*bd.sim.nb+4*bd.ID , 1:1]
end
function pdot(bd::Body)
  """return p of body by accessing sim q"""
  return bd.sim.qdot[3*bd.sim.nb+4(bd.ID-1)+1:3*bd.sim.nb+4*bd.ID , 1:1]
end
function pddot(bd::Body)
  """return p of body by accessing sim q"""
  return bd.sim.qddot[3*bd.sim.nb+4(bd.ID-1)+1:3*bd.sim.nb+4*bd.ID , 1:1]
end

function A(bd::Body)
  """returns rotation matrix of the body, calculated from p"""
  return P2A(p(bd))
end
function G(bd::Body)
  """returns G matrix of a body, calculated from p"""
  return G(p(bd))
end
function Gdot(bd::Body)
  """returns G matrix of a body, calculated from p"""
  return G(pdot(bd))
end

"""take a 3x1 ωbar and a body and returns global pdot (for setting up IC's)"""
function ωbar2pdot(bd::Body, ω::Array)
  pdot = .5*G(bd)'*ω
end


#----------------------pseudo setter methods------------------------------------
function set_r!(bd::Body, r::Array)
 bd.sim.q[3*(bd.ID-1)+1:3*bd.ID , 1:1] = r
end
function set_rdot!(bd::Body, rdot::Array)
 bd.sim.qdot[3*(bd.ID-1)+1:3*bd.ID , 1:1] = rdot
end
function set_rddot!(bd::Body , r::Array)
 bd.sim.qddot[3*(bd.ID-1)+1:3*bd.ID , 1:1] = rddot
end

function set_p!(bd::Body , p::Array)
 bd.sim.q[3*bd.sim.nb+4(bd.ID-1)+1:3*bd.sim.nb+4*bd.ID , 1:1] = p
end
function set_pdot!(bd::Body, r::Array)
 bd.sim.qdot[3*bd.sim.nb+4(bd.ID-1)+1:3*bd.sim.nb+4*bd.ID , 1:1] = pdot
end
function set_pddot!(bd::Body , r::Array)
 bd.sim.qddot[3*bd.sim.nb+4(bd.ID-1)+1:3*bd.sim.nb+4*bd.ID , 1:1] = pddot
end

"""handles a request to add a point, and return the ID of a specified point"""
function _addPoint(bd::Body, point::Array)
  #check if point is already in bd.points, and return it's index
  for pointID in 1:size(bd.pts)[2] #number of points in point matrix
    if pt(bd,pointID) == point #point already exists in list
      return pointID
    end
  end
  #add point to array, return point index
  addNewPoint(bd,point)
  return size(bd.points)[2] #new point stored at end of array

end


"""adds a new point to a bodies' point matrix"""
function addNewPoint(bd::Body, point::Array)
 bd.pts = [bd.pts point]
end

"""adds a new point to a bodies' point matrix"""
function addPoint(bd::Body, point::Array)
 bd.pts = [bd.pts point]
end
