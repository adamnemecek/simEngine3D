type SnapShot
  q::Array{Float64}       #[7nb x 1]array of system generalized coordinates = [r;p]
  qdot::Array{Float64}    #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  qddot::Array{Float64}   #[7nb x 1]array of system generalized coordinates = [rdot;pdot]
  t::Float64
end
