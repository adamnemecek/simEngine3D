# ME 751 SimEngine3D benchmark - 4 bar linkage
#Alex Dawson-Elli

#simulate the 4 bar linkage, and plot the displacement of B0

#includes
include("../src/SimEngine3D.jl")
using SimEngine3D ; SE3D = SimEngine3D;  #alias

#define the simulation
sim = SE3D.Sim(4)

#----------------------add bodies to system-------------------------------------
#add the ground body in the simulation as body 1
SE3D.addGround!(sim) #g is in -z

#------link 1----------
#initial kinematics
L =.5  #half the length of the rod
p = [1 ; 0; 0; 0]
r = [0 ; 0 ; .5]

#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = mass/12*(2*L)^2
J[2,2] = mass/12*(2*L)^2
J[3,3] = 0

#make a body for the pendulum ,body #2
rod1 = SE3D.Body(sim,2,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, rod1, r, p)

#------link 2---------
#initial kinematics
L =.5  #half the length of rod
p = [1 ; 0; 0; 0]
r = [0; .5 ; 1]

#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = mass/12*(2*L)^2
J[2,2] = 0
J[3,3] = mass/12*(2*L)^2

#make a body for the pendulum ,body #2
rod2 = SE3D.Body(sim,3,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, rod2, r, p)

#------link 3---------
#initial kinematics
L =.5  #half the length of rod
p = [1 ; 0; 0; 0]
r = [0; 1 ; .5 ]

#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = mass/12*(2*L)^2
J[2,2] = mass/12*(2*L)^2
J[3,3] = 0

#make a body for the pendulum ,body #2
rod3 = SE3D.Body(sim,4,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, rod3, r, p)

#---------------------add constraints to the system-----------------------------

#----------joint 1 - rj1 -----------------
#i is ground, j is rod1
#define points
Pi = [0 0 0]'  #this is the default point
Qj = [0 0 -L]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 0 1]'
cj_head = [1 0 0]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
#SE3D.addPoint(sim.bodies[1] , Pi)
SE3D.addPoint(sim.bodies[1] , ai_head) #index 3
SE3D.addPoint(sim.bodies[1] , bi_head) #index 4
SE3D.addPoint(sim.bodies[2] , Qj)      #index 2
SE3D.addPoint(sim.bodies[2] , cj_head) #index 3

#hardcode the indecies
PiID = 1; QjID = 2; ai_headID = 3; bi_headID = 4; cj_headID = 3
#add kinematic constraints
rj1 = SE3D.rj(sim,sim.bodies[1],sim.bodies[2],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rj1)

#----------joint 2 - rj2 -----------------
#i is rod1 , j is rod2
#define points
Pi = [0 0 L]'
Qj = [0 -L 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 0 1]'
cj_head = [1 0 0]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
SE3D.addPoint(sim.bodies[2] , Pi)      #index 4
SE3D.addPoint(sim.bodies[2] , ai_head) #index 5
SE3D.addPoint(sim.bodies[2] , bi_head) #index 6
SE3D.addPoint(sim.bodies[3] , Qj)      #index 2
SE3D.addPoint(sim.bodies[3] , cj_head) #index 3

#hardcode the indecies
PiID = 4; QjID = 2; ai_headID = 5; bi_headID = 6; cj_headID = 3
#add kinematic constraints
rj2 = SE3D.rj(sim,sim.bodies[2],sim.bodies[3],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rj2)

#----------joint 3 - rj3 -----------------
#i rod2 , j is rod3
#define points
Pi = [0 L 0]'
Qj = [0 0 L]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 0 1]'
cj_head = [1 0 0]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
SE3D.addPoint(sim.bodies[3] , Pi)      #index 4
SE3D.addPoint(sim.bodies[3] , ai_head) #index 5
SE3D.addPoint(sim.bodies[3] , bi_head) #index 6
SE3D.addPoint(sim.bodies[4] , Qj)      #index 2
SE3D.addPoint(sim.bodies[4] , cj_head) #index 3

#hardcode the indecies
PiID = 4; QjID = 2; ai_headID = 5; bi_headID = 6; cj_headID = 3
#add kinematic constraints
rj3 = SE3D.rj(sim,sim.bodies[3],sim.bodies[4],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rj3)

#----------joints 4 & 5 - cd_y and cd_z -----------------
#i rod1 , j is rod3
#define points
cy = [0 1 0]'
cz = [0 0 1]'
Pi = [0 0 -L]'
Qj = [0 0 -L]'
f(t) = 1

#add points to bodies
SE3D.addPoint(sim.bodies[2] , Pi)      #index 7
SE3D.addPoint(sim.bodies[4] , Qj)      #index 4

#hardcode the indecies
PiID = 7; QjID = 4

#add kinematic constraints
cd_y = SE3D.cd(sim,sim.bodies[2],sim.bodies[4],PiID,QjID,cy,f)
cd_z = SE3D.cd(sim,sim.bodies[2],sim.bodies[4],PiID,QjID,cz)
SE3D.addConstraint!(sim,cd_y)
SE3D.addConstraint!(sim,cd_z)

#---------------------set system initial velocities-----------------------------
SE3D.set_rdot!(sim.bodies[2], [0,.5,0])  #rod1 moving .5m/s to the right

#-----------------------initialize simulation-----------------------------------
SE3D.initForAnalysis(sim)

#determine remainder of system velocities
SE3D.setInitialVelocities(sim)

#---------------------perform Dynamics Analysis---------------------------------
tstart = 0
tstop = 10
δt = .01

tic()
hist = SE3D.DynamicsAnalysis(sim,tstart,tstop,δt)
toc()

#------------------------------plot---------------------------------------------
penID = 2 #body 2
path = "./unitySim/Assets/data/fourBarLinkage/q_rot.csv"
#SE3D.plotReactionTorque(penID,hist)
#SE3D.exportKinematicsToCSV(hist ,path , SE3D.Rx(-pi/2))
