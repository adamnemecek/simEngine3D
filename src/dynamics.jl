#dynamics
"""
Dynamics Analysis  - for systems that are underconstrainted, determine the kinematics
q,qdot,qddot and the reaction forces given the constraints and the externally applied
forces
"""
function DynamicsAnalysis(sim::Sim,tStart,tStop,δt = .001)
  if nDOF(sim) = 0
    error("dynamic analysis should happen with degress of freedom present")
  end

  #make sure level zero and level one constraints are consistent
  checkInitialConditions(sim)

  #find the initial qddot and lagrange multipliers, based on supplied level zero and one checkInitialConditions
  findInitialL2conditions(sim)

  #setup history
  tgrid = tStart:δt:tStop # type unit range
  hist = History(sim,tgrid)
  tInd = 1;

  #iterate through time grid and solve dynamics equations
  for instant in tgrid
    #update time
    sim.t = instant

    if sim.t = tStart #save snapshot @ t = 0
      snapshot(sim)
    end

    if tInd = 2  #use BDF of order to seed future solutions
      BDF(sim,1,δt,tInd,hist)
    end

    if tInd >= 3 #perform full integration with BDF order 2
      BDF(sim,2,δt,tInd,hist)
    end
    #store simulation state snapshot
    snapShot(sim,2,tInd)
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
function findInitialL2conditions(sim::Sim) # 10.19 slide 16
  #construct the required system level state variables
  #build LHS
  buildM(sim)
  buildɸk_r(sim)
  buildɸk_p(sim)
  buildP(sim)
  buildJᵖ(sim)
  z12 = zeros(3*sim.nb,4*sim.nb)
  z13 = zeros(3*sim.nb,sim.nb)
  z21 = z12'
  z31 = z13'
  z34 = zeros(sim.nb,sim.nc)
  z43 = z34'
  z44 = zeros(sim.nc,sim.nc)

  LHS = [sim.M    z12      z13    sim.ɸk_r';
         z21      sim.Jᵖ   sim.P' sim.ɸk_p';
         z31      sim.P    z33    z34      ;
         sim.ɸk_r sim.ɸk_p z43    z44      ]

  #build RHS
  buildF(sim)
  buildτh(sim)
  build𝛾p(sim)
  build𝛾k(sim)

  RHS = [sim.F ; sim.τh ; sim.𝛾p ; sim.𝛾k]

  #solve for the level 2 compontents @ t0
  L2 = LHS \ RHS
  sim.q =  L2[1:7*sim.nb, 1:1]
  sim.λp = L2[(7*sim.nb + 1):(7*sim.nb + 1 + sim.nc_p), 1:1]
  sim.λk = L2[(7*sim.nb + sim.nc_p + 1):end, 1:1]

  #build the reaction force vector for archiving in history  buildFʳ(sim)
  buildFʳ(sim)
  buildnbarʳ(sim)

end

"""
solves for the system state information, at a single time step, by following The
6 steps outlined in 10.19 slide 17, using the quazi newton ψ matrix
"""
BDF(sim::Sim,BDForder::Int64,δt::Float64,tInd::Int64,hist::History ) #10.19 slide 17
  #interation constants
  mxInterations = 100;
  tolerance = 1e-2;

  #step 0: done, system should be incremented and up to derivative
  #step 1: compute the position and velocity using BDF and most recent qddot estimates
    #handle the different order of BDF
  if BDForder == 1  # 10.14 slide 19 table
    β₀ =  1
    α₁ = -1
    Cʳdot = -α₁*rdot(hist,tInd - 1)
    Cᵖdot = -α₁*pdot(hist,tInd - 1)
    Cʳ =    -α₁*r(hist,tInd - 1) + β₀*δt*Cʳdot
    Cᵖ =    -α₁*p(hist,tInd - 1) + β₀*δt*Cᵖdot
  end

  if BDForder == 2
    β₀ =  2/3
    α₁ = -4/3
    α₂ =  1/3
    Cʳdot = -α₁*rdot(hist,tInd - 1)  - α₂*rdot(hist,tInd - 2)
    Cᵖdot = -α₁*pdot(hist,tInd - 1)  - α₂*pdot(hist,tInd - 2)
    Cʳ =    -α₁*r(hist,tInd - 1) - α₂*r(hist,tInd - 2) + β₀*δt*Cʳdot
    Cᵖ =    -α₁*p(hist,tInd - 1) - α₂*p(hist,tInd - 2) + β₀*δt*Cᵖdot

  end
  #step 2: compute the residual in the nonlinear system. i.e. evaluate g(z) at ν
  #step 3: solve the linear system ψ*Δz = -g
  #step 4: Improve the quality of the approximate solution z(ν+1) = z(ν) + Δz
  #step 5: increment ν, check if  satisfactory convergence is reached
  #step 6:
  end
