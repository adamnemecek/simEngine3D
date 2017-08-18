# ME 751 HW 6 problem 3
#Alex Dawson-Elli

#find the equalibrium position of a spring mass rotational pendulum by waiting
#for the system to damp out.
#includes
include("../src/SimEngine3D.jl")
using SimEngine3D ; SE3D = SimEngine3D;  #alias

#define the simulation
sim = SE3D.Sim(2)

#add the ground body in the simulation as body 1
SE3D.addGround!(sim)

#define the pendulum object and add to system
t = 0  #initial time
L = 2  #half the length of the pendulum
theta = pi/4*cos(2*t) #initial configuration
rot = SE3D.Ry(pi/2)*SE3D.Rz(theta) #account for pi/2 rotation
p = SE3D.A2P(rot)
r = [0,L*sin(theta),-L*cos(theta)]

#make a body for the pendulum ,body #2
pendulum = SE3D.Body(sim,2)

#add to simulation in initial configuration
SE3D.addBody!(sim, pendulum, r, p)

#add points to bodies to use for kinematic constraints (see 9.26 - slide 32)
#i is ground, j is pendulum
Pi = [0 0 0]'  #this is the default point
Qj = [-L 0 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 0 -1]'  #ai and bi vectors define the yz plane of body i (ground)
bi_head = [0 1 0]'
cj_head = [0 0 1]'  #vector defining the axis of rotation in the LRF of bodyj
bodyj_x = [1 0 0]'
#SE3D.addPoint(sim.bodies[1] , Pi)
SE3D.addPoint(sim.bodies[1] , ai_head) #index 3
SE3D.addPoint(sim.bodies[1] , bi_head) #index 4
SE3D.addPoint(sim.bodies[2] , Qj)      #index 2
SE3D.addPoint(sim.bodies[2] , cj_head) #index 3
SE3D.addPoint(sim.bodies[2] , bodyj_x) #index 4

#hardcode the indecies
PiID = 1; QjID = 2; ai_headID = 3; bi_headID = 4; cj_headID = 3 ; bodyj_xID = 4
#add kinematic constraints
rev = SE3D.rj(sim,sim.bodies[1],sim.bodies[2],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rev)

#----------make a rotational Spring-Damper-Actuator and add to system-----------
#points
ai_head = [1 0 0]'
aj_head = [0 0 -1]'
bi_head = [0 1 0]'
bj_head = [1 0 0]'

SE3D.addPoint(sim.bodies[1] , ai_head) #index 5
SE3D.addPoint(sim.bodies[1] , bi_head) #index 6
SE3D.addPoint(sim.bodies[2] , aj_head) #index 5
SE3D.addPoint(sim.bodies[2] , bj_head) #index 6


#hardcode the indecies
ai_headID = 5; aj_headID = 5; bi_headID = 6; bj_headID = 6;
k = 2 ; c = 10 ; Θ₀ = -pi/4 ; h = 0

#add RSDA element
rsda1 = SE3D.RSDA(sim,sim.bodies[1],sim.bodies[2],ai_headID,aj_headID,bi_headID,bj_headID,k,Θ₀,c)
SE3D.addSDA!(sim,rsda1)

#initialize simulation
SE3D.initForAnalysis(sim)

#perform kinematic analysis
tstart = 0
tstop = 10
δt = .01

tic()
hist = SE3D.DynamicsAnalysis(sim,tstart,tstop,δt, "NR")
toc()

#plot
penID = 2 #body 2
#SE3D.plot2DKinematics(penID,hist)

#update data file for unity plot
path = "./unitySim/Assets/data/penRSDAequilibrium/q_rot.csv"
#SE3D.exportKinematicsToCSV(hist ,path , SE3D.Rx(-pi/2))
