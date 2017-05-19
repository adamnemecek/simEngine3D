#includes

type Sim
  """
  the system contains all the bodies, constraints and forces required to perform
  spacial kinematics and dynamics on a time grid.
  """
  #counters
  nb      #number of bodies in the system
  nc      #number of constraint equations in the system
  nc_k    #number of kinematic and driving contraints in the system
  nc_p    #number of euler parameter normalization constraints == nb
  t       #current time of the system

  #objects
  bodies  #[nb x 1] array of body objects in the system
  cons    #[nCons x 1]array of constraint objects in system

  #state
  q       #[7nb x 1]array of system generalized coordinates = [r;p]
  qdot    #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  qddot   #[7nb x 1]array of system generalized coordinates = [rdot;pdot]

  #constraint matricies
  ɸ       #[nc_k x 1]array of system non-linear equations of constraint
  ɸf      #[nc x 1]  array of system equations of constraint - including euler params
  ɸ_r
  ɸ_p
  ɸ_q

  #dynamics
end
