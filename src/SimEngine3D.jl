
#SimEngine3D provides a set of simulation tools to construct and analyze
#he Kinematics and dynammics of 3D spacial systems

module SimEngine3D

#system objects
include("./sim.jl")
include("./body.jl")
#------------------------geometric Constraints----------------------------------
#low level constraints
include("./GCons/dp1.jl")
include("./GCons/dp2.jl")
include("./GCons/d.jl")
include("./GCons/cd.jl")

#intermediate Constraints
include("./GCons/b1.jl")
include("./GCons/b2.jl")

#high level constraints
include("./GCons/sj.jl")
include("./GCons/cj.jl")
include("./GCons/tj.jl")
include("./GCons/rj.jl")
include("./GCons/uj.jl")

#euler parameter normalization constraints
include("./GCons/ep.jl")

#absolute constraints
include("./GCons/ground.jl")
#-------------------------kinematics, dynamics, utilities-----------------------
include("./kinematics.jl")
include("./dynamics.jl")
include("./utils.jl")

#-------------------------plotting and exporting--------------------------------
include("./snapshot.jl")
include("./plot.jl")

#---------------------------public interface------------------------------------

end
