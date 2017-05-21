#kinematics - all functions required for kinematics analysis.
"""
kinematics analysis - generates the time history (q, qdot, qddot) of a
system between tStart and tStop
Inputs: tStart, tStop, δt
output: hist - an array of history terms
"""
function KinematicsAnalysis(sim::Sim,tStart,tStop,δt = .001)
  if nDOF(sys) > 0
    error("system is underconstrainted, and therefore kinematic analysis cannot continue")
  end
  if Rank(sys.ɸF_q) < nDOF(sys)
    error("simulation is starting in a singular configuration ")
  end

  #tgrid
  tgrid = tStart,δt,tStop # type unit range

  #setup history
  hist = Array{SnapShot}(lenght(tgrid))
  histInd = 1;
  #iterate through grid and solve equations
  for instant in tgrid
    #update time
    sim.t = instant

    #don't do position analysis at t = 0
    if sim.t != tStart
      positionAnalysis(sim)
    end

    #do velocity analysis
    velocityAnalysis(sim)

    #acceleration analysis
    accelerationAnalysis(sim)

    #store simulation history
    hist[histInd] = Snapshot(sim)
    histInd += 1
  end



"""
solve the non-linear equations of constraint using an iterative newton-rapson
approach. results in solution for q at time t
"""
function positionAnalysis(sim::Sim , ϵ =1e-9 . , maxIter = 100)
  initial_q = sim.q ; posErr = 1; counter = 1;
  while posErr > ϵ
    bulidɸ(sim)
    buildɸF_q(sim)
    Δq = sim.ɸF_q \ -sim.ɸF
    sim.q += Δq
    if counter > maxIter
      error('failure to converge')
    end
  end
end

"""
solve a linear system to determine qdot at time t
"""
function velocityAnalysis(sim::Sim)
  buildɸF_q(sim) #most updated version
  buildνF(sim)
  sim.qdot = sim.ɸF_q \ sim.νF
end

"""
solve a linear system to determine qddot a ttime t
"""
function accelerationAnalysis(sim::Sim)
  build𝛾F(sim)
  sim.qdot = sim.ɸF_q \ sim.𝛾F
end
