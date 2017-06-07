type History
  #kinematics of interest
  q::Array{Float64}       #[7nb x t]array of system generalized coordinates = [r;p]
  qdot::Array{Float64}    #[7nb x t]array of system generalized coordinates = [rdot;pdot]
  qddot::Array{Float64}   #[7nb x t]array of system generalized coordinates = [rdot;pdot]
  #dynamics of interest
  Fʳ::Array{Float64}       #[3nb x  t] vector of system reaction forces
  nbarʳ::Array{Float64}    #[3nb x  t] vector of system reaction torques
  λk::Array{Float64}       #[nc_k  x t] vector of system lagrange multipliers
  λp::Array{Float64}       #[nc_p  x t] vector of system lagrange multipliers
  νerror::Array{Float64}   #[1 x t]     vector of velocity constraint violations

  t::FloatRange{Float64}
  function History(sim::Sim,tgrid::FloatRange{Float64})
    #pass in sim and range of t over which function is evaluated
    q =     zeros(length(sim.q),length(tgrid))
    qdot =  zeros(length(sim.qdot),length(tgrid))
    qddot = zeros(length(sim.qddot),length(tgrid))

    Fʳ    = zeros(length(sim.Fʳ),length(tgrid))  #might not be initialized at time of hist creation
    nbarʳ = zeros(length(sim.nbarʳ),length(tgrid))
    λk    = zeros(length(sim.λk),length(tgrid))
    λp    = zeros(length(sim.λp),length(tgrid))

    νerror = zeros(1,length(tgrid))

    new(q,qdot,qddot,Fʳ,nbarʳ,λk, λp, νerror, tgrid)
  end
end

"""store the system state from a single instant to the system history"""
function snapShot(sim::Sim,hist::History,tInd::Int64)
    hist.q[:,tInd] = sim.q
    hist.qdot[:,tInd] = sim.qdot
    hist.qddot[:,tInd] = sim.qddot

    hist.Fʳ[:,tInd] = sim.Fʳ
    hist.nbarʳ[:,tInd] = sim.nbarʳ

    hist.νerror[1,tInd] = νerror(sim)
end

#-----------------------------extractor functions-------------------------------------
#iso-time extractor (used in dynamics for history)
r(sim::Sim, hist::History,tInd::Int64)     =     hist.q[1:3*sim.nb, tInd:tInd]
rdot(sim::Sim, hist::History,tInd::Int64)  = hist.qdot[1:3*sim.nb, tInd:tInd]
rddot(sim::Sim, hist::History,tInd::Int64) = hist.qddot[1:3*sim.nb, tInd:tInd]

p(sim::Sim, hist::History,tInd::Int64)     =    hist.q[3*sim.nb+1:end, tInd:tInd]
pdot(sim::Sim, hist::History,tInd::Int64)  = hist.qdot[3*sim.nb+1:end, tInd:tInd]
pddot(sim::Sim, hist::History,tInd::Int64) = hist.qddot[3*sim.nb+1:end, tInd:tInd]

getλk(sim::Sim, hist::History,tInd::Int64) = hist.λk[: , tInd:tInd]
getλp(sim::Sim, hist::History,tInd::Int64) = hist.λP[: , tInd:tInd]

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
  rdot = sim.q[1:3*sim.nb,1:1] ; pdot = sim.q[3*sim.nb + 1:end, 1:1]
  νi = sim.ɸk_r*rdot + sim.ɸk_p*pdot
  Verrors = νi - sim.νk
  return norm(Verrors)
end
