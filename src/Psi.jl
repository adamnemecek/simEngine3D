#this file contains functions specific to the calculation of the full ψ matrix,
#sed in solving either the dynamics problem with NR or MN approach


"""
build the full ψ Matrix for use in NR and MN methods, as specified on 10.19 slide 6
"""
function buildFullψ(sim::Sim, h::Float64, β₀::Float64 )
  #assemble components if ψUL
  z11 = zeros(3*sim.nb,3*sim.nb)
  z12 = zeros(3*sim.nb,4*sim.nb)

  ɸλ_qq = buildɸλ_qq(sim)
  Inertia  = [sim.M  z12;
              z12'   sim.Jᵖ]
   Jpp = [z11 z12;
          z12' buildJpPddot_p(sim)]

  #build ψUL (Upper Left), which is the only part of ψFull different from QN
  ψUL = h^2*β₀^2 * (ɸλ_qq + Jpp)  + Inertia  #10.19 slide 6

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



"""builds a portion of Ψ upper left hand corner that come from   ɸλ_qq"""
function buildɸλ_qq(sim::Sim)
  #start by flattening the system - we neglect ground terms for now
  flatCons = flattenCons(sim)

  #init ɸλ_qq
  ɸλ_qq = zeros(7*sim.nb,7*sim.nb)

  #loop to fill in ɸλ_qq with information from each constraint
  for (ID,con) in enumerate(flatCons)
    λ = sim.λk[ID + 6,1]  #This a very brittle solution to geting around ground!
    bi = con.bodyi
    bj = con.bodyj

    #calculate all ϕλ_qq values for the current constraint, this is a symmetric matrix

    #       ri                          rj                           pi                           pj
    ϕλ_riri = λ*ϕ_riri(con) ;   ϕλ_rirj = λ*ϕ_rirj(con)  ;   ϕλ_ripi = λ*ϕ_ripi(con)  ;   ϕλ_ripj = λ*ϕ_ripj(con)   #ri
    ϕλ_rjri = ϕλ_rirj'      ;   ϕλ_rjrj = λ*ϕ_rjrj(con)  ;   ϕλ_rjpi = λ*ϕ_rjpi(con)  ;   ϕλ_rjpj = λ*ϕ_rjpj(con)   #rj
    ϕλ_piri = ϕλ_ripi'      ;   ϕλ_pirj = ϕλ_rjpi'       ;   ϕλ_pipi = λ*ϕ_pipi(con)  ;   ϕλ_pipj = λ*ϕ_pipj(con)   #pi
    ϕλ_pjri = ϕλ_ripj'      ;   ϕλ_pjrj = ϕλ_rjpj'       ;   ϕλ_pjpi = ϕλ_pipj'       ;   ϕλ_pjpj = λ*ϕ_pjpj(con)   #pj

    #the above matrix is
    #  ψ_rr | ψ_rp
    #-------|--------
    #  ψ_pr | ψ_pp

    #place value from the current constraint into the system ɸλ_qq matrix

    #ψ_rr
    ɸλ_qq[rr(bi),rr(bi)] = ϕλ_riri ; ɸλ_qq[rr(bi),rr(bj)] = ϕλ_rirj
    ɸλ_qq[rr(bj),rr(bi)] = ϕλ_rjri ; ɸλ_qq[rr(bj),rr(bj)] = ϕλ_rjrj

    #ψ_rp
    ɸλ_qq[rr(bi),pr(bi)] = ϕλ_ripi ; ɸλ_qq[rr(bi),pr(bj)] = ϕλ_ripj
    ɸλ_qq[rr(bj),pr(bi)] = ϕλ_rjpi ; ɸλ_qq[rr(bj),pr(bj)] = ϕλ_rjpj

    #ψ_pr
    ɸλ_qq[pr(bi),rr(bi)] = ϕλ_piri ; ɸλ_qq[pr(bi),rr(bj)] = ϕλ_pirj
    ɸλ_qq[pr(bj),rr(bi)] = ϕλ_pjri ; ɸλ_qq[pr(bj),rr(bj)] = ϕλ_pjrj

    #ψ_pr
    ɸλ_qq[pr(bi),pr(bi)] = ϕλ_pipi ; ɸλ_qq[pr(bi),pr(bj)] = ϕλ_pipj
    ɸλ_qq[pr(bj),pr(bi)] = ϕλ_pjpi ; ɸλ_qq[pr(bj),pr(bj)] = ϕλ_pjpj

  end
  return ɸλ_qq
end

"""builds [Jᵖpddot]p which is a term in ψ22"""  #10.19 slide 8
function buildJpPddot_p(sim::Sim)
  JpPddot_p = zeros(4*sim.nb,4*sim.nb)
  for body in sim.bodies
    a = body.J*G(p(body))*pddot(body)
    T = [0 -a' ;
         a tilde(a)]
    JpPddot_p[4(body.ID-1)+1:4*body.ID,4(body.ID-1)+1:4*body.ID] = -4*G(p(body))'*body.J*G(pddot(body)) + 4*T
  end
  return JpPddot_p
end

"""returns a list of all constraints broken down to their 4 basic constraints"""
function flattenCons(sim::Sim)
  #***ground is currently neglected***
  flatCons =  Array{Any}(0)
  for con in sim.cons
    flatCons = [flatCons ; recCon(con)]
  end
  return flatCons
end

"""recursive con search"""
function recCon(con)
  conList = Array{Any}(0)
  if typeof(con) == ground  #do nothing for grounds
    return Array{Any}(0)
  end
  if isdefined(con, :subGCs) #current GCon is high or intermediate level
    for GC in con.subGCs
      conList = [conList ; recCon(GC)]
    end
    return conList
  end
  return con #we have hit a base Gcon
end
