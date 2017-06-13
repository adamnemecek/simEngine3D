#kinematics - all functions required for kinematics analysis.
"""
Kinematics Analysis - generates the time history (q, qdot, qddot) of a
system between tStart and tStop
Inputs: tStart, tStop, δt
output: hist - an array of history terms
"""
function kinematicsAnalysis(sim::Sim,tStart,tStop,δt = .001)
  if nDOF(sim) > 0
    error("system is underconstrainted, and therefore Kinematic Analysis cannot continue")
  end
  buildɸF_q(sim)
  if rank(sim.ɸF_q) < sim.nc
    error("simulation is starting in a singular configuration ")
  end

  #tgrid
  tgrid = tStart:δt:tStop # type unit range

  #setup history
  hist = History(sim,tgrid)
  tInd = 1
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

    #store simulation state snapshot
    snapShot(sim,hist,tInd)
    tInd += 1
  end

  return hist
end

"""
solve the non-linear equations of constraint using an iterative newton-rapson
approach. results in solution for q at time t
"""
function positionAnalysis(sim::Sim , ϵ = 1e-9 , maxIter = 100) #9.29 S69
  initial_q = sim.q ; ΔqNorm = 1; counter = 1
  while  ΔqNorm  > ϵ
    buildɸF(sim)
    buildɸF_q(sim)
    Δq = sim.ɸF_q \ -sim.ɸF
    sim.q += Δq
     ΔqNorm = norm(Δq)
    counter += 1
    if counter > maxIter
      error("failure to converge")
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
solve a linear system to determine qddot at time t
"""
function accelerationAnalysis(sim::Sim)
  build𝛾F(sim)
  sim.qddot = sim.ɸF_q \ sim.𝛾F
end

"""
used to try and get systems out of a singularity (bricard system)
"""
function posJiggle(sim::Sim , ϵ = 1e-7 , maxIter = 1000) #9.29 S69
  initial_q = sim.q ; ΔqNorm = 1; counter = 1
  while  ΔqNorm  > ϵ
    buildɸF(sim)
    buildɸF_q(sim)

    extraCon = zeros(1,7*sim.nb); extraCon[1,5] = 1;
    RHS = [sim.ɸF ; r(sim.bodies[2])[2] + .5]
    LHS = [sim.ɸF_q ; extraCon]
    #sim.ɸF_q could be non-invertable
    U,Σ,V = svd(LHS)
    Σ = diagm(Σ)
    Jplus = V*(Σ'*(Σ*Σ')^-1)*U'

    Δq = Jplus * - RHS
    sim.q += Δq
     ΔqNorm = norm(Δq)
    counter += 1
    if counter > maxIter
      error("failure to converge")
    end
  end
end
