#inverse dynamics
#kinematics - all functions required for kinematics analysis.
"""
Inverse Dynamics Analysis - for a system with prescribed time history (fully
kinematically constrained system) determine the forces and torques required to
produce that prescribed motion, on the specified time interval
Inputs: tStart, tStop, δt
output: hist - an array of kinematic and kinetic quantities
"""
function InverseDynamicsAnalysis(sim::Sim,tStart,tStop,δt = .001) #10.10 slide 7-10
  if nDOF(sim) > 0
    error("system is underconstrainted, and therefore Inverse Dynamic Analysis cannot continue")
  end
  buildɸF_q(sim)
  if rank(sim.ɸF_q) < sim.nc
    error("simulation is starting in a singular configuration ")
  end

  #construct system constant matricies
  buildM(sim)
  #tgrid
  tgrid = tStart:δt:tStop # type unit range

  #setup history
  hist = History(sim,tgrid)
  tInd = 1;
  #iterate through grid and solve equations
  for instant in tgrid
    #update time
    sim.t = instant
    #-----------------step 1 - solve Kinematics for qddot-----------------------
    #don't do position analysis at t = 0
    if sim.t != tStart
      positionAnalysis(sim)
    end

    #do velocity analysis
    velocityAnalysis(sim)

    #acceleration analysis
    accelerationAnalysis(sim)

    #-----------step 2 - solve for lagrange multipliers using N-E dynamics-----
    buildP(sim)
    build

    #------------------step 3 - calculate reaction forces----------------------
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
function positionAnalysis(sim::Sim , ϵ = 1e-9 , maxIter = 100)
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
