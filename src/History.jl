type History
  #kinematics of interest
  q::Array{Float64}       #[7nb x t]array of system generalized coordinates = [r;p]
  qdot::Array{Float64}    #[7nb x t]array of system generalized coordinates = [rdot;pdot]
  qddot::Array{Float64}   #[7nb x t]array of system generalized coordinates = [rdot;pdot]
  nb::Int64               #number of bodies in the system
  #dynamics of interest
  rForces::Array{Float64}  #[3nb x nc_k x t] vector of system reaction forces
  rTorques::Array{Float64} #[3nb x nc_k x t]  vector of system reaction torques
  νerror::Array{Float64}   #[1 x t]     vector of velocity constraint violations

  tgrid::FloatRange{Float64}

  #constructor
  function History(sim::Sim,tgrid::FloatRange{Float64})
    #pass in sim and range of t over which function is evaluated
    q =     zeros(length(sim.q),length(tgrid))
    qdot =  zeros(length(sim.qdot),length(tgrid))
    qddot = zeros(length(sim.qddot),length(tgrid))
    nb = sim.nb

    rForces  = zeros(3*sim.nb, sim.nc_k, length(tgrid))
    rTorques = zeros(3*sim.nb, sim.nc_k, length(tgrid))

    νerror = zeros(1,length(tgrid))

    new(q,qdot,qddot,nb,rForces,rTorques,νerror,tgrid)
  end
end

"""store the system state from a single instant to the system history"""
function snapShot(sim::Sim,hist::History,tInd::Int64)
    hist.q[:,tInd] = sim.q
    hist.qdot[:,tInd] = sim.qdot
    hist.qddot[:,tInd] = sim.qddot

    hist.rForces[:,:,tInd] = sim.rForces
    hist.rTorques[:,:,tInd] = sim.rTorques

    hist.νerror[1,tInd] = νerror(sim)
end

#-----------------------------extractor functions-------------------------------------
#iso-time extractor (used in dynamics for history)
r(hist::History,tInd::Int64)     =     hist.q[1:3*hist.nb, tInd:tInd]
rdot(hist::History,tInd::Int64)  = hist.qdot[1:3*hist.nb, tInd:tInd]
rddot(hist::History,tInd::Int64) = hist.qddot[1:3*hist.nb, tInd:tInd]

p(hist::History,tInd::Int64)     =    hist.q[3*hist.nb+1:end, tInd:tInd]
pdot(hist::History,tInd::Int64)  = hist.qdot[3*hist.nb+1:end, tInd:tInd]
pddot(hist::History,tInd::Int64) = hist.qddot[3*hist.nb+1:end, tInd:tInd]

#---------------for plotting-------------------
#in all the below functions, force or torque goes along the rows, and time, the columns
"""returns the time history of net reactions forces on each body"""
function Fʳ(hist::History)
  Fr = zeros(3*hist.nb,length(hist.t))
  for instant in 1:length(hist.tgrid)
    Fr[:,instant] = rowSum(hist.rForces[:,:,instant])
  end
  return Fr
end

"""returns the time history of net reactions forces on each body"""
function nbarʳ(hist::History)
  nbarʳ = zeros(3*hist.nb,length(hist.tgrid))
  for instant in 1:length(hist.tgrid)
    nbarʳ[:,instant] = rowSum(hist.rTorques[:,:,instant])
  end
  return nbarʳ
end

"""returns the time history a reaction torque caused by λ_i"""
function rTorque(hist::History, bdID::Int64 , λID::Int64)
  rTorque = zeros(3,length(hist.tgrid))
  for instant in 1:length(hist.tgrid)
    rTorque[:,instant] = hist.rTorques[3*(bdID-1)+1:3*bdID, λID:λID, instant]
  end
  return rTorque
end

"""returns the time history a reaction torque caused by λ_i"""
function rForce(hist::History, bdID::Int64 , λID::Int64)
  rForce = zeros(3,length(hist.tgrid))
  for instant in 1:length(hist.tgrid)
    rForce[:,instant] = hist.rForces[3*(bdID-1)+1:3*bdID, λID:λID, instant]
  end
  return rForce
end

#iso-body extractor (used in plotting functions)
# r(hist::History,tInd::Int64)         hist.q[1:3*sim.nb, tInd:tInd]
# rdot(hist::History,tInd::Int64)   hist.qdot[1:3*sim.nb, tInd:tInd]
# rddot(hist::History,tInd::Int64) hist.qddot[1:3*sim.nb, tInd:tInd]
#
# p(sim::Sim, hist::History,tInd::Int64)         hist.q[3*sim.nb+1:end, tInd:tInd]
# pdot(sim::Sim, hist::History,tInd::Int64)   hist.qdot[3*sim.nb+1:end, tInd:tInd]
# pdot(sim::Sim, hist::History,tInd::Int64)  hist.qddot[3*sim.nb+1:end, tInd:tInd]


#---------------------------calculator functions---------------------------------------
#to calculate values to be sotred in history from simulation state
function νerror(sim::Sim)
  rdot = sim.qdot[1:3*sim.nb,1:1] ; pdot = sim.qdot[3*sim.nb + 1:end, 1:1]
  νi = sim.ɸk_r*rdot + sim.ɸk_p*pdot
  Verrors = νi - sim.νk
  return sqrt(norm(Verrors))
end
