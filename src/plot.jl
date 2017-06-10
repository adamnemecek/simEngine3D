#imports
  using Plots ; gr()
  using DataFrames


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


function plotReactionTorque(bodyID,hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and torque from history
  t = hist.t
  ni = hist.nbarʳ[3*(bodyID-1)+1:3*bodyID , :]

  #setup plot variables
  titles = "torque on body $(bodyID)"
  labels = ["nx" "ny" "nz"]
  ylabels = "n*m"
  xlabels = "t"

  #execute plot
  plot(t, ni', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)

end

function plotVelocityViolations(bodyID,hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and torque from history
  t = hist.t
  vError = hist.νerror[1:1 , :]

  #setup plot variables
  titles = "vError of system"
  labels = "vError"
  ylabels = "(m/s)^2"
  xlabels = "t"

  #execute plot
  plot(t, vError', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)
end



"""function used for simulating results in unity"""
function exportKinematicsToCSV(hist , path, A = eye(3))
  #perform any principle rotations required to get q from our default coordinates
  # (z-up, gravity is negative z) into coordinates for unity


  q = copy(hist.q)
  r = q[1:3*hist.nb, :]
  p = q[3*hist.nb+1:end, :]
  if A != eye(3)  #we need to perform a rotation to get things to look right in unit

    #rotate position vectors
    for bID in 1:hist.nb
      for instant in 1:size(q)[2]
        r[3*(bID-1)+1:3*bID , instant:instant] = A*r[3*(bID-1)+1:3*bID , instant:instant]
      end
    end

    #rotate orientation vecotrs
    for bID in 1:hist.nb
      for instant in 1:size(q)[2]
        p[4*(bID-1)+1:4*bID , instant:instant] = A2P(A*P2A(p[4*(bID-1)+1:4*bID , instant:instant]))
      end
    end
  end
  #put q back together again
  q = [r;p]
  #clean up really small entries, idk how well small numbers are parsed
  for row in 1:size(q)[1]
    for col in 1:size(q)[2]
      if abs(q[row,col]) < .00001
        q[row,col] = 0;
      end
    end
  end
  kinematics = convert(DataFrame, q)
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
