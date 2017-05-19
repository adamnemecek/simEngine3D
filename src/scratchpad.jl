include()

type Body
  t::Float64
  a::Array
  c::fobo
end

# how to prpoerly make a minimal argument constructor
type BubbleGum
  f::Int     #flavor
  s::Array   #size
  c::Float64 #cost

  function BubbleGum(f::Int)
    new(f,[1,2],.227)
  end
end
