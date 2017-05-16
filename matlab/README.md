#simEngine3D

system architecture - The simulation engine is hierarchically structured
with the system being composed of bodies and constraints, each body is 
itself composed of a set of markers. Bodies are defined by their parameters
such as position (r) orientation (p) with euler parameters, as well as mass
and rotational inertia properties. constraints have a type, and act between
two bodies. Makers are defined by position only, as they are points, 
and are defined relative to the local reference frame of the body the belong
to. 
