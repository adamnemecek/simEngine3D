type History
  q::Array{Float64}       #[7nb x 1]array of system generalized coordinates = [r;p]
  qdot::Array{Float64}    #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  qddot::Array{Float64}   #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  t::FloatRange{Float64}
  function History(sim::Sim,tgrid::FloatRange{Float64})
    #pass in sim and range of t over which function is evaluated
    q =     zeros(length(sim.q),length(tgrid))
    qdot =  zeros(length(sim.qdot),length(tgrid))
    qddot = zeros(length(sim.qddot),length(tgrid))
    new(q,qdot,qddot,tgrid)
  end
end

"""store the system state from a single instant to the system history"""
function snapShot(sim::Sim,hist::History,tInd::Int64)
    hist.q[:,tInd] = sim.q
    hist.qdot[:,tInd] = sim.qdot
    hist.qddot[:,tInd] = sim.qddot

end
