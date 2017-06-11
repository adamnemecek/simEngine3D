# ME 751 SimEngine3D benchmark - 4 bar linkage
#Alex Dawson-Elli

#simulate the 4 bar linkage, and plot the displacement of B0

#includes
include("../src/SimEngine3D.jl")
using SimEngine3D ; SE3D = SimEngine3D;  #alias

#define the simulation
sim = SE3D.Sim(4)

#add the ground body in the simulation as body 1
SE3D.addGround!(sim) #g is in -z

#----------------------define link1 and add to system---------------------------
#initial kinematics
L =.5  #half the length of the pendulum
theta = pi/2 #initial configuration
rot = SE3D.Ry(pi/2)*SE3D.Rz(theta) #account for pi/2 rotation
p = SE3D.A2P(rot)
r = [0; L*sin(theta);-L*cos(theta)]

#mass properties - pendulum 1 (body 2)
ρ = 7800; l = 4; w = .05
volume = l*w^2
mass = volume*ρ
J = eye(3)
J[1,1] = mass/12*(w^2 + w^2)
J[2,2] = mass/12*(l^2 + w^2)
J[3,3] = J[2,2]

#make a body for the pendulum ,body #2
pen1 = SE3D.Body(sim,2,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, pen1, r, p)

#add points to bodies to use for kinematic constraints (see 9.26 - slide 32)
#i is ground, j is pendulum
Pi = [0 0 0]'  #this is the default point
Qj = [-L 0 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 0 -1]'  #ai and bi vectors define the yz plane of body i (ground)
bi_head = [0 1 0]'
cj_head = [0 0 1]'  #vector defining the axis of rotation in the LRF of bodyj
#SE3D.addPoint(sim.bodies[1] , Pi)
SE3D.addPoint(sim.bodies[1] , ai_head) #index 3
SE3D.addPoint(sim.bodies[1] , bi_head) #index 4
SE3D.addPoint(sim.bodies[2] , Qj)      #index 2
SE3D.addPoint(sim.bodies[2] , cj_head) #index 3

#hardcode the indecies
PiID = 1; QjID = 2; ai_headID = 3; bi_headID = 4; cj_headID = 3
#add kinematic constraints
rev = SE3D.rj(sim,sim.bodies[1],sim.bodies[2],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rev)

#--------------define the 2nd pendulum object and add to system-----------------
#initial kinematics
L = 1  #half the length of the pendulum
theta = 0 #initial configuration
rot = SE3D.Ry(pi/2)*SE3D.Rz(theta) #account for pi/2 rotation
p = SE3D.A2P(rot)
r = [0; 4 + L*sin(theta);-L*cos(theta)] #4 accounts for the length of pen1

#mass properties - pendulum 1 (body 2)
ρ = 7800; l = 2; w = .05
volume = l*w^2
mass = volume*ρ
J = eye(3)
J[1,1] = mass/12*(w^2 + w^2)
J[2,2] = mass/12*(l^2 + w^2)
J[3,3] = J[2,2]

#make a body for the pendulum 2 ,body #3
pen2 = SE3D.Body(sim,3,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, pen2, r, p)

#add points to bodies to use for kinematic constraints (see 9.26 - slide 32)
#i is ground, j is pendulum
Pi = [2 0 0]'  #this is the default point
Qj = [-L 0 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [1 0 0]'  #ai and bi vectors define the yz plane of body i (ground)
bi_head = [0 1 0]'
cj_head = [0 0 1]'  #vector defining the axis of rotation in the LRF of bodyj

#SE3D.addPoint(sim.bodies[1] , Pi)
SE3D.addPoint(sim.bodies[2] , Pi) #index 4
SE3D.addPoint(sim.bodies[2] , ai_head) #index 5
SE3D.addPoint(sim.bodies[2] , bi_head) #index 6
SE3D.addPoint(sim.bodies[3] , Qj)      #index 2
SE3D.addPoint(sim.bodies[3] , cj_head) #index 3


#hardcode the indecies
PiID = 4; QjID = 2; ai_headID = 5; bi_headID = 6; cj_headID = 3
#add kinematic constraints
rev = SE3D.rj(sim,sim.bodies[2],sim.bodies[3],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rev)


#-----------------------initialize simulation-----------------------------------
SE3D.initForAnalysis(sim)

#perform Inverse Dynamics Analysis
tstart = 0
tstop = 10
δt = .01
tic()
hist = SE3D.DynamicsAnalysis(sim,tstart,tstop,δt)
toc()
#plot
penID = 2 #body 2
#SE3D.plotReactionTorque(penID,hist)
#SE3D.exportKinematicsToCSV("./")
