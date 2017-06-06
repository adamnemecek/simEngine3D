#imports
  using Plots ; gr()


function plot2DKinematics(bodyID,hist)
  """
  plots the position vs time for the bodies requested via their ID Number
  """
    #get time , position, velocity, acceleration for specified body
  t = hist.t
  pos = hist.q[3*(bodyID-1)+1:3*bodyID ,:]
  vel = hist.qdot[3*(bodyID-1)+1:3*bodyID ,:]
  acc = hist.qddot[3*(bodyID-1)+1:3*bodyID ,:]

  #reformat to x,y,z for plotting routine
  x = [pos[1:1,:] ; vel[1:1,:] ; acc[1:1,:]]'
  y = [pos[2:2,:] ; vel[2:2,:] ; acc[2:2,:]]'
  z = [pos[3:3,:] ; vel[3:3,:] ; acc[3:3,:]]'

#setup plot variables
titles = ["body$(bodyID) position"  "velocity"  "acceleration"]
labels = ["x" "x" "x" "y" "y" "y" "z" "z" "z"]
ylabels = ["m" "m/s" "m/s²"]
xlabels = ["t" "t" "t"]

#execute plot
plot(t,[x y z], layout = (3,1), title = titles, label = labels)
plot!(ylabel = ylabels)
end

"""function used for simulating results in unity"""
function exportKinematicsToCSV(hist , path)
  using DataFrames
  kinematics = convert(DataFrame, hist.q)
  writetable(path,kinematics)
end


# """plots the 3D path of a body specified by bodyID"""
# function plot3DKinematics(bodyID,history)
#     # initialize the attractor
#   n = 1500
#   dt = 0.02
#   σ, ρ, β = 10., 28., 8/3
#   x, y, z = 1., 1., 1.
#
#   # initialize a 3D plot with 1 empty series
#   plt = path3d(1, xlim=(-25,25), ylim=(-25,25), zlim=(0,50),
#                   xlab = "x", ylab = "y", zlab = "z",
#                   title = "Lorenz Attractor", marker = 1)
#
#   # build an animated gif, saving every 10th frame
#   @gif for i=1:n
#       dx = σ*(y - x)     ; x += dt * dx
#       dy = x*(ρ - z) - y ; y += dt * dy
#       dz = x*y - β*z     ; z += dt * dz
#       push!(plt, x, y, z)
#   end every 10
# end
