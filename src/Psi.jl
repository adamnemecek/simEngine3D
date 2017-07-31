"""
this file contains functions specific to the calculation of the full ψ matrix,
used in solving either the dynamics problem with NR or MN approach
"""


"""
build the full ψ Matrix for use in QN and MN methods, as specified on 10.19 slide 6
"""
function buildFullψ(sim::Sim, h::Float64, β₀::Float64 )
  #assemble components if ψUL
  z11 = zeros(3*sim.nb,3*sim.nb)
  z12 = zeros(3*sim.nb,4*sim.nb)

  ɸλ_qq = buildɸλ_qq(sim)
  Inertia  = [sim.M  z12;
              z12'   sim.Jᵖ]
  Jpp = [z11 z12;
         z12' buildJpPddotp(sim)]

  #build ψUL (Upper Left), which is the only part of ψFull different from QN
  ψUL = h^2*β₀^2 * [ɸλ_qq + Jpp]  + Inertia  #10.19 slide 6

  #------------------build the fullΨ matrix--------------------------------
  #extract ψnn from the previously constructed ψUL matrix
  ψ11 = ψUL[1:3*sim.nb,1:3*sim.nb]                     #upper left
  ψ12 = ψUL[1:3*sim.nb,3*sim.nb+1:7*sim.nb]            #upper right
  ψ21 = ψUL[3*sim.nb+1:7*sim.nb,1:3*sim.nb]            #lower left
  ψ22 = ψUL[3*sim.nb+1:7*sim.nb,3*sim.nb+1:7*sim.nb]   #lower right

  #initialize constant matricies useful in construction of ψ
  z31 = zeros(sim.nb,3*sim.nb)
  z33 = zeros(sim.nb,sim.nb)
  z44 = zeros(sim.nc_k,sim.nc_k)
  z43 = zeros(sim.nc_k,sim.nb)


  ψFull =   [ ψ11      ψ21'       z31'   sim.ɸk_r';
              ψ21      ψ22        sim.P' sim.ɸk_p';
              z31      sim.P      z33    z43'     ;
              sim.ɸk_r sim.ɸk_p   z43    z44      ]

  return ψFull
end






"""builds a portion of Ψ upper left hand corner that come from ɸqq"""
function buildɸλ_qq(sim::Sim,)
  #start by flattening the system
end

"""builds [Jᵖpddot]p which is a term in ψ22"""
function buildJpPddotp(sim::Sim)
end

"""returns a list of all constraints broken down to their 4 basic constraints"""
function flattenCons(sim::Sim, )
end
