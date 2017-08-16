#imports
  #using Plots ; gr()
  using Plots ; pyplot()
  using DataFrames


function plot2DKinematics(bodyID,hist, sbarP = [0 0 0]')
  """
  plots the position vs time for the bodies requested via their ID Number
  """
  #get time , position, velocity, acceleration for specified body
  t = hist.tgrid

  r = hist.q[3*(bodyID-1)+1:3*bodyID ,:]
  rdot = hist.qdot[3*(bodyID-1)+1:3*bodyID ,:]
  rddot = hist.qddot[3*(bodyID-1)+1:3*bodyID ,:]

  p = hist.q[3*hist.nb+4(bodyID-1)+1:3*hist.nb+4*bodyID ,:]
  pdot = hist.qdot[3*hist.nb+4(bodyID-1)+1:3*hist.nb+4*bodyID ,:]
  pddot = hist.qddot[3*hist.nb+4(bodyID-1)+1:3*hist.nb+4*bodyID ,:]

  #initialize kinematics matricies
  pos = zeros(3,size(hist.tgrid)[1])
  vel = zeros(3,size(hist.tgrid)[1])  #point on body
  acc = zeros(3,size(hist.tgrid)[1])  #point on body

  #calculate translational kinematics of the specified points
  for instant in 1:size(hist.tgrid)[1]
    #calculate angular kinematic terms
    A = P2A(p[:,instant])
    ωbar = pdot2ωbar(p[:,instant], pdot[:,instant])
    ωdotbar = pddot2ωdotbar(p[:,instant], pddot[:,instant])

    #form translational kinematics
    pos[:,instant] = r[:,instant] + A*sbarP
    vel[:,instant] = rdot[:,instant] + A*tilde(ωbar)*sbarP
    acc[:,instant] = rddot[:,instant] + (A*tilde(ωdotbar) + A*tilde(ωbar)*tilde(ωbar))*sbarP
  end

  #reformat to x,y,z for plotting routine
  x = [pos[1:1,:] ; vel[1:1,:] ; acc[1:1,:]]'
  y = [pos[2:2,:] ; vel[2:2,:] ; acc[2:2,:]]'
  z = [pos[3:3,:] ; vel[3:3,:] ; acc[3:3,:]]'

  #setup plot variables
  titles = ["body$(bodyID) point $(sbarP) position"  "velocity"  "acceleration"]
  labels = ["x" "x" "x" "y" "y" "y" "z" "z" "z"]
  ylabels = ["m" "m/s" "m/s²"]
  xlabels = ["t" "t" "t"]

  #execute plot
  plot(t,[x y z], layout = (3,1), title = titles, label = labels)
  plot!(ylabel = ylabels)
end


function plotNetReactionTorque(bodyID, hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and torque from history
  t = hist.tgrid
  nbar = nbarʳ(hist)
  nibar = nbar[3*(bodyID-1)+1:3*bodyID , :]

  #setup plot variables
  titles = "torque on body $(bodyID)"
  labels = ["nx" "ny" "nz"]
  ylabels = "n*m"
  xlabels = "t"

  #execute plot
  plot(t, nibar', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)

end

function plotReactionTorque(bodyID, λID, hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and torque from history
  t = hist.tgrid
  nibar  = rTorque(hist, bodyID , λID)

  #setup plot variables
  titles = "torque on body $(bodyID) from constraint $(λID)"
  labels = ["nx" "ny" "nz"]
  ylabels = "n*m"
  xlabels = "t"

  #execute plot
  plot(t, nibar', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)

end


function plotappliedTorque(bodyID, hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and torque from history
  t = hist.tgrid
  nibar  = hist.nbarᵃ[3*(bodyID-1)+1:3*bodyID , :]

  #setup plot variables
  titles = "torque on body $(bodyID) from RSDA's and TSDA's"
  labels = ["nx" "ny" "nz"]
  ylabels = "n*m"
  xlabels = "t"

  #execute plot
  plot(t, nibar', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)

end

function plotappliedForce(bodyID, hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and force from history
  t = hist.tgrid
  Fi  = hist.Fᵃ[3*(bodyID-1)+1:3*bodyID , :]

  #setup plot variables
  titles = "forces applied to body by TSDA's"
  labels = ["fx" "fy" "fz"]
  ylabels = "n"
  xlabels = "t"

  #execute plot
  plot(t, Fi', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)

end

function plotVelocityViolations(bodyID,hist)
  """
  plots the reaction torque vs time for the bodies requested via their ID Number
  """
  #get time and torque from history
  t = hist.tgrid
  vError = hist.νerror[1:1 , :]

  #setup plot variables
  titles = "vError of system"
  labels = "vError"
  ylabels = "(m/s)^2"
  xlabels = "t"

  #execute plot
  plot(t, vError', title = titles, label = labels, ylabel = ylabels, xlabel = xlabels)
end

"""function used to export simulations results to unity as a csv"""
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
      q[row,col] = round(q[row,col] , 5)
      if abs(q[row,col]) < .00001
        q[row,col] = 0;
      end
    end
  end
  kinematics = convert(DataFrame, q)
  writetable(path,kinematics)
end
