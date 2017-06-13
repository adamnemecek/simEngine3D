#inverse dynamics
#kinematics - all functions required for kinematics analysis.
# """
# Inverse Dynamics Analysis - for a system with prescribed time history (fully
# kinematically constrained system) determine the forces and torques required to
# produce that prescribed motion, on the specified time interval
# Inputs - tStart, tStop, δt
# output - hist - an array of kinematic and kinetic quantities
# """

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
  tInd = 1
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

    #-----------step 2 - solve for lagrange multipliers using N-E dynamics------
    #phi_q already up to date
    buildP(sim)
    buildJᵖ(sim)
    buildF(sim)
    buildτh(sim)

    rddot = sim.qddot[1:3*sim.nb, 1:1]
    pddot = sim.qddot[(3*sim.nb + 1):end, 1:1]

    #from  10.10 slide 12 , in matrix form
    LHS = [sim.ɸk_r'   zeros(3*sim.nb,sim.nb) ;
           sim.ɸk_p'             sim.P'        ]
    RHS = [sim.F - sim.M*rddot ;
          sim.τh -  sim.Jᵖ*pddot]

    sim.λF = LHS \ RHS
    sim.λk = sim.λF[1:sim.nc_k, 1:1]
    sim.λp = sim.λF[sim.nc_k+1:end, 1:1]
    #------------------step 3 - calculate reaction forces-----------------------
    buildrForces(sim)
    buildrTorques(sim)

    #store simulation state snapshot
    snapShot(sim,hist,tInd)
    tInd += 1
  end

  return hist
end
