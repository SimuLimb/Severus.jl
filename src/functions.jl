"""
    severusdir()

# Description 

This function simply returns the string for the Severus path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function severusdir()
    joinpath(@__DIR__, "..")
end

"""
Function to smooth the boundary 

"""

function smoothboundary(Vs,Fs,Eb; Z_thickness_distal=nothing)
    if isnothing(Z_thickness_distal)
        indCurve_boundary = edges2curve(Eb)
        constrained_points =nothing 
        Vs_smooth = smoothmesh_hc(Fs[indCurve_boundary],Vs, 50, 0.1, 0.5; constrained_points=constrained_points)
    else
        Vs = smoothmesh_laplacian(Eb,Vs, 25, 0.5)   
        ZF = [mean(Vs[f])[3] for f in Fs]
        logicSmooth = ZF.<(minimum(ZF)+Z_thickness_distal)
        constrained_points = elements2indices(boundaryedges(Fs[logicSmooth]))
        Vs_smooth = smoothmesh_hc(Fs[logicSmooth],Vs, 50, 0.1, 0.5; constrained_points=constrained_points)
    end       
    return Vs_smooth
end 

"""
function for adding taper to the skin

"""

function taperfunction(Vs,V_tibia,indCurve_tibia,taperFraction,Z_taper_offset,taperThreshold)
    Zs = [v[3] for v in Vs]
    min_Z = minimum(Zs)
    z_level_taper = min_Z + Z_taper_offset
    # Vt_end = V_tibia[indCurve_tibia]
    D= zeros(length(Vs))
    Dp = zeros(length(Vs))
    indMin1 = Int.(zeros(length(Vs)))
    R = fill(Point{3, Float64}(0.0, 0.0, 0.0), length(Vs))
    for (i,v) in enumerate(Vs)
        if v[3] < z_level_taper
            D[i] = ((v[3] - z_level_taper) / (min_Z-z_level_taper))^2        
            d = Vector{Float64}(undef,length(indCurve_tibia))
            for (i_d,j) in enumerate(indCurve_tibia)
                d[i_d] = sqrt((v[1]-V_tibia[j][1])^2 + (v[2]-V_tibia[j][2])^2)
            end
            dp, index = findmin(d) 
            if dp <= taperThreshold
                Dp[i] = taperThreshold
                indMin1[i] = Int(index)
            else 
                Dp[i] = dp
                indMin1[i] = Int(index)
            end
        
        R[i] = V_tibia[indCurve_tibia[Int(index)]] .- v
        #   push!(R1,r)
        end
    end
    Dp ./= maximum(Dp)
    
    R = R.*D.*Dp.*taperFraction
    #process skin taper
    Vs1 = Vs
    Vs_new = [Point{3,Float64}(v[1]+r[1],v[2]+r[2],v[3]) for (v,r) in zip(Vs1,R)]
    return Vs_new 
end 

"""
Function for bezierclose
"""
function bezierclose(Fs,Vs_new,v1,v2,indCurve ; numPointBezier = numPointBezier, p_distal = nothing, numEdgeSmoothSteps = 25, λ = 0.5, pointSpacing=nothing)
    
    if isnothing(pointSpacing)
        pointSpacing = mean(edgelengths(Fs,Vs_new))
    end

    #Bezier
    NV = vertexnormal(Fs,Vs_new; weighting=:area) # May be inefficient to do this to all vertices (/faces)
    NV_indCurve = NV[indCurve] # Vertex normals for the curve
    NE_indCurve = Vector{Point{3,Float64}}(undef,length(indCurve))
    NB_base = Vector{Point{3,Float64}}(undef,length(indCurve))
    m = length(indCurve)
    for i in eachindex(indCurve)
        ii = indCurve[i]
        jj = indCurve[mod1(i+1,m)]
        NE_indCurve[i] = normalizevector(Vs_new[jj] - Vs_new[ii])

        NB_base[i] = cross(NE_indCurve[i],NV_indCurve[i])
    end

    for _ in 1:1:numEdgeSmoothSteps
        for i in eachindex(indCurve)
            i_prev = mod1(i-1,m)
            i_next = mod1(i+1,m)
            NB_base[i] = ((1.0-λ) .* NB_base[i]) .+  λ.*(NB_base[i_prev] .+ NB_base[i_next])./2.0    
            NB_base[i] = normalizevector(NB_base[i])
        end
    end

    VV_B = []
    Vb = Vector{Point{3,Float64}}()
    #pdistal nothing
    if isnothing(p_distal)     
        vz = [v[3] for v in Vs_new]
        len = abs(maximum(vz) - minimum(vz))
        v_low = mean(Vs_new[indCurve])
        p_distal = Point{3,Float64}(v_low[1],v_low[2],v_low[3] -  len*0.03)
        # p_distal = v_low .- len*0.1
    end
     
    for i in eachindex(indCurve)
        nBase = NB_base[i]
        nEnd = normalizevector(Point{3,Float64}(Vs_new[indCurve[i]][1]-p_distal[1],Vs_new[indCurve[i]][2]-p_distal[2],0.0))
        pb = [Vs_new[indCurve[i]], Vs_new[indCurve[i]] + v1/3*nBase, p_distal + v2/3*nEnd,p_distal]
        vb = nbezier(pb,numPointBezier)
        push!(VV_B,vb)
        append!(Vb,vb)
    end   
    Fb = grid2surf(Vb,length(indCurve);periodicity=(false,true),face_type=:quad)   
    Fbq = quad2tri(Fb,Vb)   
    F_n,V_n,_ = joingeom(Fs,Vs_new,Fbq,Vb)
    F_n,V_n,_ = mergevertices(F_n,V_n)
    # Determine number of points to request to closely match desired point spacing
    np = spacing2numvertices(F_n,V_n,pointSpacing) 
    # Remesh with desired number of points 
    Fn,Vn = ggremesh(F_n,V_n; nb_pts=np)
    
    return Fn,Vn
end 
#snap to z z_level

function set_z_level_boundary!(Fn,Vn,z_level;indCurve=nothing)
    #this works only for a single boundary 
    if isnothing(indCurve)
        Eb = boundaryedges(Fn)
        indCurve = edges2curve(Eb)
        indCurve = indCurve[1:end-1]
    end 
    for i in indCurve
        Vn[i] = Point{3,Float64}(Vn[i][1],Vn[i][2],z_level)
    end
    # return Vn 
end 

function getboundaryset(F,V)
    Eb = boundaryedges(F)
    Gb = meshgroup(Eb)
    numGroups = maximum(Gb)
    z_level = zeros(numGroups)
     boundary_set = Vector{Vector{Int}}()
    z = [v[3] for v in V]
    for g in 1:1:numGroups
        ind = edges2curve(Eb[Gb.==g])
        ind = ind[1:end-1]  
        push!(boundary_set,ind)     
        z_level[g] = mean(z[ind])

    end
    # _,indMax = findmax(z_level)
    # if indMax==1
    #     ind1 = edges2curve(Eb[Gb.==1])
    #     ind1 = ind1[1:end-1]
    #     ind2= edges2curve(Eb[Gb.==2])
    #     ind2 = ind2[1:end-1]
    # else
    #     ind1 = edges2curve(Eb[Gb.==2])
    #     ind1 = ind1[1:end-1]
    #     ind2= edges2curve(Eb[Gb.==1])
    #     ind2 = ind2[1:end-1]
    # end
    return boundary_set,z_level
end
