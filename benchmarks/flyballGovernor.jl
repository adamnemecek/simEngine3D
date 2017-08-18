# ME 751 SimEngine3D benchmark - flyballGovernor
#Alex Dawson-Elli

#simulate a flyballGovernor mechanism over time

#includes
include("../src/SimEngine3D.jl")
using SimEngine3D ; SE3D = SimEngine3D;  #alias

#define the simulation
sim = SE3D.Sim(5)

#----------------------add bodies to system-------------------------------------
#add the ground body in the simulation as body 1
SE3D.addGround!(sim) #g is in -z

#general mass properties for system
ρ = 3000 #kg/m³
h = 1; r = .005 ; v = pi*r^2*h
#m = v*ρ ; Jmaj = m/12*(3*r^2 + h^2); Jmin = .5*m*r^2  #mass, major axis inertia, minor axis inertia
m = 10 ; Jmaj = m/12*(3*r^2 + h^2); Jmin = .5*m*r^2
J = zeros(3,3)

#------link 1----------
#initial kinematics
L =.5  #half the length of the rod
r = [0 ; 0 ; .5]
p = [1 ; 0; 0; 0]

#mass properties
J[1,1] = Jmaj
J[2,2] = Jmaj
J[3,3] = Jmin

#add body to system
L1 = SE3D.Body(sim,2,m,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L1, r, p)

#------link 2----------
#initial kinematics
α = pi/6
r = [0 ; -.05 ; 1] + .5*[0 ; -cos(α) ; -sin(α)]
p = SE3D.A2P(SE3D.Rx(α))


#mass properties - slender rod
J[1,1] = Jmaj
J[2,2] = Jmin
J[3,3] = Jmaj

#add body to system
L2 = SE3D.Body(sim,3,m,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L2, r, p)

#------link 3----------
#initial kinematics
r = [0 ; .05 ; 1] + .5*[0 ; cos(α) ; -sin(α)]
p = SE3D.A2P(SE3D.Rx(-α))


#mass properties - slender rod
J[1,1] = Jmaj
J[2,2] = Jmin
J[3,3] = Jmaj

#add body to system
L3 = SE3D.Body(sim,4,m,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, L3, r, p)

#------base----------
#initial kinematics
r = [0 ; 0 ; .5]
p = [1 ; 0; 0; 0]

#mass properties - hollow right circular cylinder
Ro = .05 ; Ri =.005 ; h = .1
#m = pi*ρ*h*(Ro^2 - Ri^2)
m  = 2
Jmaj = m/2*(3*Ro^2 + 3*Ri^2 + h^2) ;   Jmin = m/2*(Ro^2 + Ri^2)
J = zeros(3,3)
J[1,1] = Jmaj
J[2,2] = Jmaj
J[3,3] = Jmin

#add body to system
base = SE3D.Body(sim,5,m,J)

#add to simulation in initial configuration
SE3D.addBody!(sim, base, r, p)

#---------------------add constraints to the system-----------------------------


#----------joint 1 - rj1 -----------------
#i is ground, j is L1
#define points
Pi = [0 0 0]'  #this is the default point
Qj = [0 0 -.5]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [1 0 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 1 0]'
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
Pi = [0 -.05 .5]'
Qj = [0 .5 0]' #Pi and Qj are the points which are coincident in the sphirical joint
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
#i L1 , j is L3
#define points
Pi = [0 .05 .5]'
Qj = [0 -.5 0]' #Pi and Qj are the points which are coincident in the sphirical joint
ai_head = [0 1 0]'  #ai and bi vectors define the yz plane of body i
bi_head = [0 0 1]'
cj_head = [1 0 0]'  #vector defining the axis of rotation in the LRF of bodyj

#add points to bodies
SE3D.addPoint(sim.bodies[2] , Pi)      #index 7
SE3D.addPoint(sim.bodies[2] , ai_head) #index 8
SE3D.addPoint(sim.bodies[2] , bi_head) #index 9
SE3D.addPoint(sim.bodies[4] , Qj)      #index 2
SE3D.addPoint(sim.bodies[4] , cj_head) #index 3

#hardcode the indecies
PiID = 7; QjID = 2; ai_headID = 8; bi_headID = 9; cj_headID = 3
#add kinematic constraints
rj3 = SE3D.rj(sim,sim.bodies[2],sim.bodies[4],PiID,QjID,ai_headID,bi_headID,cj_headID)
SE3D.addConstraint!(sim,rj3)

#----------joint 4 - tj4 -----------------
#i is L1, J is b
#define points
Pi = [0 0 0]'  #Pi and Qj are on the translational axis
Qj = [0 0 0]'
ai_head = [1 0 0]'
bi_head = [0 1 0]'
aj_head = [0 1 0]'
cj_head = [0 0 1]'

#add points to bodies
SE3D.addPoint(sim.bodies[2] , Pi)      #index 10
SE3D.addPoint(sim.bodies[2] , ai_head) #index 11
SE3D.addPoint(sim.bodies[2] , bi_head) #index 12
SE3D.addPoint(sim.bodies[5] , Qj)      #index 2
SE3D.addPoint(sim.bodies[5] , aj_head) #index 3
SE3D.addPoint(sim.bodies[5] , cj_head) #index 4

#hardcode the indecies
PiID = 10; QjID = 2; ai_headID = 11; bi_headID = 12; aj_headID = 3 ; cj_headID = 4
#add kinematic constraints
tj4 = SE3D.tj(sim,sim.bodies[2],sim.bodies[5],PiID,QjID,ai_headID,bi_headID,aj_headID,cj_headID)
SE3D.addConstraint!(sim,tj4)

#---------------------------set up TSDA's---------------------------------------

# #-------------TSDA 1--------------
#points
Pi_head = [0 -.05 0]'
Qj_head = [0 0 0]'

SE3D.addPoint(sim.bodies[5] , Pi_head) #index 5
SE3D.addPoint(sim.bodies[3] , Qj_head) #index 4


#hardcode the indecies
Pi_headID = 5; Qj_headID = 4;
k = 1000 ; c = 20 ; l₀ = .5

#add RSDA element
rsda1 = SE3D.TSDA(sim,sim.bodies[5],sim.bodies[3],Pi_headID,Qj_headID,k,l₀,c)
SE3D.addSDA!(sim,rsda1)

# #-------------TSDA 2--------------
#points
Pi_head = [0 .05 0]'
Qj_head = [0 0 0]'

SE3D.addPoint(sim.bodies[5] , Pi_head) #index 6
SE3D.addPoint(sim.bodies[4] , Qj_head) #index 4


#hardcode the indecies
Pi_headID = 6; Qj_headID = 4;
k = 1000 ; c = 20 ; l₀ = .5

#add RSDA element
rsda2 = SE3D.TSDA(sim,sim.bodies[5],sim.bodies[4],Pi_headID,Qj_headID,k,l₀,c)
SE3D.addSDA!(sim,rsda2)


# #---------------------set system initial velocities-----------------------------
SE3D.set_pdot!(sim.bodies[2], SE3D.ωbar2pdot(sim.bodies[2],[0 0 2]'))  #L2 spinning at 2 rad/s
ω = 4 #rad/seconds
#SE3D.set_rdot!(sim.bodies[4], [-ω*(.05 + cos(pi/6)) 0 0]' )  #L2 spinning at 2 rad/s
#SE3D.set_rdot!(sim.bodies[3], [ω*(.05 + cos(pi/6)) 0 0]' )


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
path = "./unitySim/Assets/data/flyballGovernor/q_rot.csv"
#SE3D.plotReactionTorque(penID,hist)
#SE3D.exportKinematicsToCSV(hist ,path, SE3D.Rx(-pi/2))
