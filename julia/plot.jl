#imports
using Gadfly

plot2DKinematics(bodyID::Array{Int},history::Array{SnapShot}; variable ="position" )
  """
  plots the position vs time for the bodies requested via their ID Number
  """
  #get time
  t = Array{Float64}
  for snapshot in history
    hcat(t,snapshot.q)
  end

  #which kinematic quantity should we look at
  var = Array{Float64}
  if variable == "pos"  || "position"
    for snapshot in history
      var = hcat(var,history.q)
    end
  elseif variable == "vel" || "velocity"
    for snapshot in history
      var = hcat(var,history.q)
    end
  elseif variable == "acc" || "acceleration"
    for snapshot in history
      var = hcat(var,history.q)
    end
  end

#grab only the body we are interested in
 pltVar = var[3*(bodyID-1)+1:3*bodyID , :]

#plot the x,y,z position, acceleration or vel of specifed body
plot(layer(x=t, y = pltVar[1,:], Geom.line),
     layer(x=t, y = pltVar[2,:], Geom.line),
     layer(x=t, y = pltVar[3,:], Geom.line),
     Guide.xLabel("time")
     Guide.yLabel("position"));
end
