

# type Body
#   t::Float64
#   a::Array
#   c::fobo
# end

# how to prpoerly make a minimal argument constructor
type BubbleGum
  f::Int     #flavor
  s::Array   #size
  c::Float64 #cost

  function BubbleGum(f::Int)
    self = new(f,[1,2],.227)
    setup(self)
    return self

  end
end

function setup(bb::BubbleGum)
  bb.c = bb.c^2
end
