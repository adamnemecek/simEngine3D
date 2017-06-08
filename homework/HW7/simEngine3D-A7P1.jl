# ME 751 HW 7 problem 1
#Alex Dawson-Elli

#problem  1 - perform an Inverse Dyanamics analysis and plot the reaction torque
#as a function of time

include("../../src/SimEngine3D.jl")
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
r = [0; L*sin(theta);-L*cos(theta)]

#mass properties
ρ = 7800; l = 4; w = .05
volume = l*w^2
mass = volume*ρ
J = eye(3)
J[1,1] = mass/12*(w^2 + w^2)
J[2,2] = mass/12*(l^2 + w^2)
J[3,3] = J[2,2]

#make a body for the pendulum ,body #2
pendulum = SE3D.Body(sim,2,mass,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, pendulum, r, p)

#add points to bodies to use for kinematic constraints (see 9.26 - slide 32)
#i is ground, j is pendulum
Pi = [0 0 0]'  #this is the default point
Qj = [-L 0 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i (ground)
bi_head = [0 0 -1]'
cj_head = [0 0 1]'  #vector defining the axis of rotation in the LRF of bodyj
bodyj_x = [1 0 0]'  #used for definition of the dot product
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

#add driving constraint using dp1 to specify pendulum angle
f(t) = cos((pi*cos(2*t))/4 - pi/2)
fdot(t) = ((pi*sin(2*t)*sin((pi*cos(2*t))/4 - pi/2))/2)
fddot(t) =(pi*cos(2*t)*sin((pi*cos(2*t))/4 - pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 - pi/2))/4)
# f(t) = cos((pi/4)*cos(2*t))
# fdot(t) = ((pi/2)*sin(2*t)*sin((pi/4)*cos(2*t)))
# fddot(t) = pi*cos(2*t)sin((pi/4)cos(2t)) - (((pi^2)/4)*(sin(2t))^2*cos((pi/4)cos(2t)))
drive = SE3D.dp1(sim,sim.bodies[1],sim.bodies[2],ai_headID,bodyj_xID,1,1,f,fdot,fddot )
SE3D.addConstraint!(sim,drive)

#initialize simulation
SE3D.initForAnalysis(sim)

#perform kinematic analysis
tstart = 0
tstop = 10
δt = .01


hist = SE3D.InverseDynamicsAnalysis(sim,tstart,tstop,δt)

#plot
penID = 2 #body 2
#SE3D.plotReactionTorque(penID,hist)  #develop this!
#SE3D.exportKinematicsToCSV(hist, "test.csv")
