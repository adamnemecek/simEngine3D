# ME 751 SimEngine3D benchmark - bricard Mechanism
#Alex Dawson-Elli

#simulate bricard's mechanism over time

#includes
include("../src/SimEngine3D.jl")
using SimEngine3D ; SE3D = SimEngine3D;  #alias

#define the simulation
sim = SE3D.Sim(6)

#----------------------add bodies to system-------------------------------------
#add the ground body in the simulation as body 1
SE3D.addGround!(sim) #g is in -z

#------link 1----------
#initial kinematics
L =.5  #half the length of the rod
r = [-.5 ; 0 ; 0]
p = [1 ; 0; 0; 0]


#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = 0
J[2,2] = mass/12*(2*L)^2
J[3,3] = mass/12*(2*L)^2

#add body to system
L1 = SE3D.Body(sim,2,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L1, r, p)

#------link 2----------
#initial kinematics
L =.5  #half the length of the rod
r = [-1 ; 0 ; -.5]
p = [1 ; 0; 0; 0]


#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = mass/12*(2*L)^2
J[2,2] = mass/12*(2*L)^2
J[3,3] = 0

#add body to system
L2 = SE3D.Body(sim,3,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L2, r, p)

#------link 3----------
#initial kinematics
L =.5  #half the length of the rod
r = [-1 ; .5 ; -1]
p = [1 ; 0; 0; 0]

#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = mass/12*(2*L)^2
J[2,2] = 0
J[3,3] = mass/12*(2*L)^2

#add body to system
L3 = SE3D.Body(sim,4,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L3, r, p)

#------link 4----------
#initial kinematics
L =.5  #half the length of the rod
r = [-.5 ; 1 ; -1]
p = [1 ; 0; 0; 0]

#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = 0
J[2,2] = mass/12*(2*L)^2
J[3,3] = mass/12*(2*L)^2

#add body to system
L4 = SE3D.Body(sim,5,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L4, r, p)

#------link 5----------
#initial kinematics
L =.5  #half the length of the rod
r = [0 ; 1 ; -.5]
p = [1 ; 0; 0; 0]

#mass properties - slender rod
mass = 1 #kg
J = zeros(3,3)
J[1,1] = mass/12*(2*L)^2
J[2,2] = mass/12*(2*L)^2
J[3,3] = 0

#add body to system
L5 = SE3D.Body(sim,6,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L5, r, p)

#---------------------add constraints to the system-----------------------------


#----------joint 1 - rj1 -----------------
#i is ground, j is L1
#define points
Pi = [0 0 0]'  #this is the default point
Qj = [.5 0 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [1 0 0]'
cj_head = [0 0 1]'  #vector defining the axis of rotation in the LRF of bodyj

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
#i is L1 , j is L2
#define points
Pi = [-.5 0 0]'
Qj = [0 0 .5]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 0 1]'  #ai and bi vectors define the yz plane of body i
bi_head = [1 0 0]'
cj_head = [0 1 0]'  #vector defining the axis of rotation in the LRF of bodyj

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
#i L2 , j is L3
#define points
Pi = [0 0 -.5]'
Qj = [0 -.5 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 0 1]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 1 0]'
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

#----------joint 4 - rj4 -----------------
#i L3 , j is L4
#define points
Pi = [0 .5 0]'
Qj = [-.5 0 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [1 0 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 1 0]'
cj_head = [0 0 1]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
SE3D.addPoint(sim.bodies[4] , Pi)      #index 4
SE3D.addPoint(sim.bodies[4] , ai_head) #index 5
SE3D.addPoint(sim.bodies[4] , bi_head) #index 6
SE3D.addPoint(sim.bodies[5] , Qj)      #index 2
SE3D.addPoint(sim.bodies[5] , cj_head) #index 3

#hardcode the indecies
PiID = 4; QjID = 2; ai_headID = 5; bi_headID = 6; cj_headID = 3
#add kinematic constraints
rj4 = SE3D.rj(sim,sim.bodies[4],sim.bodies[5],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rj4)

#----------joint 5 - rj5 -----------------
#i L4 , j is L5
#define points
Pi = [.5 0 0]'
Qj = [0 0 -.5]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [1 0 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 0 1]'
cj_head = [0 1 0]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
SE3D.addPoint(sim.bodies[5] , Pi)      #index 4
SE3D.addPoint(sim.bodies[5] , ai_head) #index 5
SE3D.addPoint(sim.bodies[5] , bi_head) #index 6
SE3D.addPoint(sim.bodies[6] , Qj)      #index 2
SE3D.addPoint(sim.bodies[6] , cj_head) #index 3

#hardcode the indecies
PiID = 4; QjID = 2; ai_headID = 5; bi_headID = 6; cj_headID = 3
#add kinematic constraints
rj5 = SE3D.rj(sim,sim.bodies[5],sim.bodies[6],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rj5)

#----------joint 6 - rj6 -----------------
#i L5 , j is ground
#define points
Pi = [0 0 .5]'
Qj = [0 1 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 0 1]'
cj_head = [1 0 0]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
SE3D.addPoint(sim.bodies[6] , Pi)      #index 4
SE3D.addPoint(sim.bodies[6] , ai_head) #index 5
SE3D.addPoint(sim.bodies[6] , bi_head) #index 6
SE3D.addPoint(sim.bodies[1] , Qj)      #index 5
SE3D.addPoint(sim.bodies[1] , cj_head) #index 6

#hardcode the indecies
PiID = 4; QjID = 5; ai_headID = 5; bi_headID = 6; cj_headID = 6
#add kinematic constraints - could not get to work with rj!!!
#rj6 = SE3D.sj(sim,sim.bodies[6],sim.bodies[1],PiID,QjID,ai_headID,bi_headID,cj_headID)
sj6 =  SE3D.sj(sim,sim.bodies[6],sim.bodies[1],PiID,QjID)
dp1 = SE3D.dp1(sim,sim.bodies[6],sim.bodies[1],bi_headID,cj_headID)
SE3D.addConstraint!(sim,sj6)
SE3D.addConstraint!(sim,dp1)

#-----------------------initialize simulation-----------------------------------
SE3D.initForAnalysis(sim)

#determine remainder of system velocities
SE3D.setInitialVelocities(sim)

#----------------------try to get out of this singularity!!---------------------
#SE3D.positionAnalysis(sim) #9.29 S69


#---------------------perform Dynamics Analysis---------------------------------
tstart = 0
tstop = 10
δt = .01

tic()
hist = SE3D.DynamicsAnalysis(sim,tstart,tstop,δt)
toc()

#------------------------------plot---------------------------------------------
penID = 2 #body 2
path = "./unitySim/Assets/data/bricardMechanism/q_rot.csv"
#SE3D.plotReactionTorque(penID,hist)
#SE3D.exportKinematicsToCSV(hist ,path , SE3D.Rx(-pi/2))
