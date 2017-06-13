#dynamics
"""
Dynamics Analysis  - for systems that are underconstrainted, determine the kinematics
q,qdot,qddot and the reaction forces given the constraints and the externally applied
forces
"""
function DynamicsAnalysis(sim::Sim,tStart,tStop,δt = .001)
  if nDOF(sim) == 0
    warn("dynamic analysis should happen with degress of freedom present")
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

    if tInd == 1 #save snapshot @ t = 0
      snapShot(sim,hist,tInd)
    end

    if tInd == 2  #use BDF of order to seed future solutions
      BDF(sim,1,δt,tInd,hist)
    end

    if tInd >= 3 #perform full integration with BDF order 2
      BDF(sim,2,δt,tInd,hist)
    end
    #store simulation state snapshot
      snapShot(sim,hist,tInd)
    tInd += 1
  end

  return hist
end

"""determine a set of velocities that are consistant with the L1 constraint equations
and a small set of manually specified velocities"""
function setInitialVelocities(sim::Sim)
  #manually determining a consistant set of velocities can be difficult, and there
  #are frequently an infinite set of solutions. this function determines a set of velocities
  #that minimizes the L2 norm of the velocity vector, and satisfies the vel constraints

  #collect bodies with nonzero rdot's and pdots
  rdotICbodies = Array{Any}(0) #set of bodies with non-zero rdot's (by reference)
  pdotICbodies = Array{Any}(0) #set of bodies with non-zero rdot's (by reference)
  for body in sim.bodies
    if rdot(body) != [0 0 0]'
        push!(rdotICbodies, body)
    end
    if pdot(body) != [0 0 0 0]'
      push!(rdotICbodies, body)
    end
  end

  #make simple constraint equations that require each specified velocity to take on it's specified value exactly
  ICconstraints = zeros(0,7*sim.nb) #LHS of the constraint equations
  b = zeros(0,1);  #RHS of constraint equations
  for body in rdotICbodies
    rdotConst = zeros(3,7*sim.nb);
    rdotConst[:, 3*(body.ID - 1) + 1:3*body.ID] = eye(3)
    #update rdot constraint equations
    ICconstraints = [ICconstraints ; rdotConst]
    b = [b ; rdot(body)]
  end

  for body in pdotICbodies
    pdotConst = zeros(4,7*sim.nb);
    pdotConst[:, 4*(body.ID - 1) + 1:4*bodyID] = eye(4)
    #update pdot constraint equations
    ICconstraints = [ICconstraints ; pdotConst]
    b = [b ; pdot(body)]
  end

  #build sim level matricies using initial r and p
  buildɸF(sim)
  buildɸk_r(sim)
  buildɸk_p(sim)
  buildP(sim)
  buildνk(sim)

  #formulate the linear system which is likely underconstrainted
  z = zeros(sim.nb, 3*sim.nb) ; zp = zeros(sim.nb,1)
  LHS = [sim.ɸk_r   sim.ɸk_p ;
            z       sim.P    ;
            ICconstraints    ]
  RHS = [sim.νk ; zp ; b]

  #check at this point that the rank of LHS does not exceed  7nb
  if rank(LHS) > 7*sim.nb
    error("too many velocities prescribed, the system is over constrained")
  end

  #use SVD to find the right pseudoInverse
  U,Σ,V = svd(LHS)
  Σ = diagm(Σ)
  Jplus = V*(Σ'*(Σ*Σ')^-1)*U'

  #solve for one of an infinite number of solutions for velocities
  sim.qdot = Jplus * RHS   #solution minimizes L2 norm of velocity vector

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
  rdot = sim.qdot[1:3*sim.nb,1:1] ; pdot = sim.qdot[3*sim.nb + 1:end, 1:1]
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
  Perrors = sim.P*pdot
  Ind = 0
  for Perror in Perrors
    if abs(Perror) > ϵ
      error("Euler parameter constraint $Ind is not consistant")
    end
    Ind += 1
  end

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
  z21 = zeros(4*sim.nb,3*sim.nb)
  z31 = zeros(sim.nb,3*sim.nb)
  z33 = zeros(sim.nb,sim.nb)
  z44 = zeros(sim.nc_k,sim.nc_k)
  z43 = zeros(sim.nc_k,sim.nb)

  LHS = [sim.M    z21'     z31'   sim.ɸk_r';
         z21      sim.Jᵖ   sim.P' sim.ɸk_p';
         z31      sim.P    z33    z43'     ;
         sim.ɸk_r sim.ɸk_p z43    z44      ]

  #build RHS
  buildF(sim)
  buildτh(sim)
  build𝛾p(sim)
  build𝛾k(sim)

  RHS = [sim.F ; sim.τh ; sim.𝛾p ; sim.𝛾k]

  #solve for the level 2 compontents @ t0
  L2 = LHS \ RHS
  sim.qddot =  L2[1:7*sim.nb, 1:1]
  sim.λp    =  L2[(7*sim.nb + 1):(7*sim.nb + sim.nc_p), 1:1]
  sim.λk    =  L2[(7*sim.nb + sim.nc_p + 1):end, 1:1]

  #build the reaction force vector for archiving in history
  buildFʳ(sim)
  buildnbarʳ(sim)

end

"""
solves for the system state information, at a single time step, by following The
6 steps outlined in 10.19 slide 17, using the quazi newton ψ matrix
"""
function BDF(sim::Sim,BDForder::Int64,δt::Float64,tInd::Int64,hist::History ) #10.19 slide 17
  #interation constants
  mxInterations = 100;
  ϵ = 1e-2;

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

  #step 0: prime new step from previous step
  #sim state variables remain unchanged from previous step
  #initialize constant matricies useful in construction of ψ
  z21 = zeros(4*sim.nb,3*sim.nb)
  z31 = zeros(sim.nb,3*sim.nb)
  z33 = zeros(sim.nb,sim.nb)
  z44 = zeros(sim.nc_k,sim.nc_k)
  z43 = zeros(sim.nc_k,sim.nb)

  ν = 0;
  while ν < mxInterations
    #step 1: compute the position and velocity using BDF and most recent qddot estimates
    rn  = Cʳ + β₀^2*δt^2*rddot(sim) ; rndot = Cʳdot + β₀*δt*rddot(sim)
    pn  = Cᵖ + β₀^2*δt^2*pddot(sim) ; pndot = Cᵖdot + β₀*δt*pddot(sim)
    #update sim level vars with q and qddot estimates
    sim.q = [rn ; pn] ; sim.qdot = [rndot ; pndot]

    #step 2: compute the residual in the nonlinear system. i.e. evaluate g(z) at ν
    #build matricies that depend on new q,qddot estimates
    #buildM(sim) const and already built
    buildɸF(sim)
    buildɸk_r(sim)
    buildɸk_p(sim)
    buildP(sim)
    buildJᵖ(sim)
    buildF(sim)
    buildτh(sim)

    #compute g(z)
    gn = [sim.M*rddot(sim)  + sim.ɸk_r'*sim.λk - sim.F ;
         sim.Jᵖ*pddot(sim)  + sim.ɸk_p'*sim.λk + sim.P'*sim.λp - sim.τh;
         1/(β₀^2*δt^2)*sim.ɸp;
         1/(β₀^2*δt^2)*sim.ɸk                                            ]

    #step 3: solve the linear system ψ*Δz = -g
    #compute ψ
    ψ =   [sim.M    z21'     z31'   sim.ɸk_r';
           z21      sim.Jᵖ   sim.P' sim.ɸk_p';
           z31      sim.P    z33    z43'     ;
           sim.ɸk_r sim.ɸk_p z43    z44      ]

    #solve linear system for correction Δz
    Δz = ψ \ -gn

    #step 4: Improve the quality of the approximate solution z(ν+1) = z(ν) + Δz
    z_old = [sim.qddot ; sim.λp ; sim.λk]
    z_new = z_old + Δz
    #update system variables
    sim.qddot = z_new[1:7*sim.nb, 1:1]
    sim.λp    = z_new[7*sim.nb + 1:8*sim.nb, 1:1]
    sim.λk    = z_new[8*sim.nb + 1:8*sim.nb + sim.nc_k, 1:1]
    #step 5: increment ν, check if  satisfactory convergence is reached
    if norm(Δz) < ϵ
      break
    end
    ν += 1
  end
  if ν == 100 println("max iterations exceeded") end
  #step 6: accept, and perform step 1 to update q and qdot to most recent values
  #step 1: compute the position and velocity using BDF and most recent qddot estimates
  rn  = Cʳ + β₀^2*δt^2*rddot(sim) ; rndot = Cʳdot + β₀*δt*rddot(sim)
  pn  = Cᵖ + β₀^2*δt^2*pddot(sim) ; pndot = Cᵖdot + β₀*δt*pddot(sim)
  #update sim level vars with q and qddot estimates
  sim.q = [rn ; pn] ; sim.qdot = [rndot ; pndot]

  #build the reaction force vector for archiving in history
  buildFʳ(sim)
  buildnbarʳ(sim)
end
