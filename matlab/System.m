classdef System
  %the System class represents the 3D system being modeled
  %it contains all of the bodies and and constrants that make up the
  %system.
  
  %from the attributes of the bodies and constraints, we assemble system level
  %matricies that constain of the kinematic and kinetic properties of the system.
  
  %we then use these matricies to perform analysis on the system. There are
  %three types of analysis that can be conducted on the system.
  %     1)Kinematic Analysis- describes the motion of the system without
  %         regard to the forces creating the motion, and can only happen when
  %         DOF = 0
  %     2)Inverse Dynamic Analysis - evaluates the forces which produced the motions 
  %         observed during kinematics analysis
  %     3)Dynamic Analysis - can only take place where there are one or
  %         more system level degrees of freedom. in dynamic analysis, forces
  %         applied to bodies dictate the time evolution of the system 
  
   %for a birds eye view of system organization, see the properties, and
   %the methods header file
    
  
    
    properties
    end
    
    %----------------------Methods Header----------------------------------
    % build system - 
    %    constructor
    %    addBody
    %    addConstrant
    %    assembleInitialConfig
    
    % KinematicAnalysis
    %   positionAnalysis
    %   velocityAnalysis
    %   accelerationAnalysis
    
    % InverseDynamicsAnalysis
    %  
    
    
   
    
    %% ----------------------build system --------------------------------
    
    methods(Access = public) %
    end
    
end

