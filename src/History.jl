type History
  #kinematics of interest
  q::Array{Float64}       #[7nb x 1]array of system generalized coordinates = [r;p]
  qdot::Array{Float64}    #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  qddot::Array{Float64}   #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  #dynamics of interest
  Fʳ::Array{Float64}       #[3nb x   1] vector of system reaction forces
  nbarʳ::Array{Float64}    #[3nb x   1] vector of system reaction torques

  t::FloatRange{Float64}
  function History(sim::Sim,tgrid::FloatRange{Float64})
    #pass in sim and range of t over which function is evaluated
    q =     zeros(length(sim.q),length(tgrid))
    qdot =  zeros(length(sim.qdot),length(tgrid))
    qddot = zeros(length(sim.qddot),length(tgrid))
    Fʳ    = zeros(3*sim.nb,length(tgrid))  #might not be initialized at time of hist creation
    nbarʳ = zeros(3*sim.qddot,length(tgrid))
    new(q,qdot,qddot,Fʳ,nbarʳ,tgrid)
  end
end

"""store the system state from a single instant to the system history"""
function snapShot(sim::Sim,hist::History,tInd::Int64)
    hist.q[:,tInd] = sim.q
    hist.qdot[:,tInd] = sim.qdot
    hist.qddot[:,tInd] = sim.qddot

    hist.Fʳ[:,tInd] = sim.Fʳ
    hist.nbarʳ[:,tInd] = sim.nbarʳ

end
