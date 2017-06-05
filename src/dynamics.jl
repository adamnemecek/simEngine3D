#dynamics
"""
Dynamics Analysis  - for systems that are underconstrainted, determine the kinematics
q,qdot,qddot and the reaction forces given the constraints and the externally applied
forces
"""
function DynamicsAnalysis(sim::Sim,tStart,tStop,δt = .001)  #10.19 slide 17
  if nDOF(sim) = 0
    error("dynamic analysis should happen with degress of freedom present")
  end
  #dynamic analysis iteration steps are layed out on 10.19 slide 17
  #step 0 check the initial conditions to insure that the system starts in a healthy state
  checkInitialConditions(sim)

  #find the initial qddot and lagrange multipliers, based on supplied level  zero and one checkInitialConditions
  findInitialL2conditions(sim)
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

    #-----------step 2 - solve for lagrange multipliers using N-E dynamics------
    #phi_q already up to date
    buildP(sim)
    buildJᵖ(sim)
    buildF(sim)
    buildτh(sim)

    rddot = sim.qddot[1:3*sim.nb, 1:1]
    pddot = sim.qddot[3*sim.nb + 1:end, 1:1]

    #from  10.10 slide 12 , in matrix form
    LHS = [sim.ɸk_r' , zeros(3*sim.nb,sim.nb) ;
           sim.ɸk_p' ,           sim.P'        ]
    RHS = [sim.F - sim.M*rddot ;
          sim.τh -  sim.Jᵖ*pddot]

    sim.λF = LHS / RHS
    sim.λk = sim.λF[1:sim.nc_k, 1:1]
    sim.λp = sim.λF[sim.nc_k+1:end, 1:1]
    #------------------step 3 - calculate reaction forces-----------------------
    buildFʳ(sim)
    buildnbarʳ(sim)

    #store simulation state snapshot
    snapShot(sim,hist,tInd)
    tInd += 1
  end

  return hist
end


"""make sure that the supplied initial conditions are satified at levels zero and on"""
function checkInitialConditions(sim::Sim, ϵ = .0001)
  #based on 10.19 slide 15
  #check position level constraints for consistency
  buildɸF(sim)
  Ind = 0
  for constraintValue in sim.ɸF
    if abs(constraintValue) > ϵ
      error("position constraint $Ind is not consistant")
    end
    Ind += 1
  end

  #check the velocity level constraints are satisfied
  buildνk(sim)
  buildɸk_r(sim)
  buildɸk_p(sim)
  rdot = sim.q[1:3*sim.nb,1:1] ; pdot = sim.q[3*sim.nb + 1, 1:1]
  ν₀ = sim.ɸk_r*rdot + sim.ɸk_p*pdot
  Verrors = ν₀ - sim.νk
  Ind = 0
  for Verror in Verrors
    if abs(Verror) > ϵ
      error("velocity equation of constraint $Ind is not consistant")
    end
    Ind += 1
  end

  #check euler parameter constraints
  buildP(sim)
  pdot = sim.q[3*sim.nb + 1, 1:1]
  Perror = sim.P*pdot
  Ind = 0
  for Perror in Perrors
    if abs(Perror) > ϵ
      error("Euler papameter constraint $Ind is not consistant")
    end
    Ind += 1

end

"""given intial level 0 and 1 conditions, solve for the level 2 conditions"""
function findInitialL2conditions(sim) # 10.19 slide 16
end
