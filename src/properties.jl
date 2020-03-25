using Discretizers
using Distances
using LinearAlgebra


function distances!(
    r::AbstractArray{T,3} where T,
    t::Trajectory{TA,NonPeriodic} where {TA<:AbstractAtomNames},
    ur1::AbstractUnitRange,
    ur2::AbstractUnitRange
    )
    metric = Euclidean()
    for i in 1:length(t)
        pairwise!(view(r,:,:,i), metric, view(t,ur1,i), view(t,ur2,i),dims=2)
    end
    return r
end

function distances!(
    r::AbstractArray{T,3} where T,
    t::Trajectory{TA,TC} where {TA<:AbstractAtomNames, TC<:AbstractOrthorombicCell},
    ur1::AbstractUnitRange,
    ur2::AbstractUnitRange
    )
    metric = (PeriodicEuclidean ∘ celldiag)(t)
    for i in 1:length(t)
        pairwise!(view(r,:,:,i), metric, view(t,ur1,i), view(t,ur2,i),dims=2)
    end
    return r
end

function distances!(
    r::AbstractArray{T,3} where T,
    t::Trajectory{TA,VariableCell} where {TA<:AbstractAtomNames},
    ur1::AbstractUnitRange,
    ur2::AbstractUnitRange
    )
    get_metric(c::NonPeriodic) = Euclidean()
    get_metric(c::AbstractOrthorombicCell) = (PeriodicEuclidean ∘ celldiag)(c)
    for i in 1:length(t)
        pairwise!(view(r,:,:,i), get_metric(c.cell[i]), view(t,ur1,i), view(t,ur2,i),dims=2)
    end
    return r
end


function distances(t::AbstractTrajectory, i1::Integer, i2::Integer)
    out = Vector{Float64}(undef, length(t))
    distances!(out, t, i1, i2)
    return out
end

function distances!(
    r::AbstractVector,
    t::Trajectory{TA,NonPeriodic} where {TA<:AbstractAtomNames},
    i1::Integer,
    i2::Integer
    )
    colwise!(r, Euclidean(), view(t,i1,:), view(t,i2,:))
    return r
end

function distances!(
    r::AbstractVector,
    t::Trajectory{TA,TC} where {TA<:AbstractAtomNames, TC<:AbstractOrthorombicCell},
    i1::Integer,
    i2::Integer
    )
    metric = (PeriodicEuclidean ∘ celldiag)(t)
    colwise!(r, metric, view(t,i1,:), view(t,i2,:))
    return r
end


cellvolume(t::AbstractTrajectory) = volume(t.cell)


"""
    angle(t::AbstractTrajectory, i, j, k) -> Vector{Float64}

Computes angle for atoms i, j, k in degrees
"""
function Base.angle(t::AbstractTrajectory, i, j, k)
    @warn "angle does not understand periodicity"
    r1 = view(t,i,:) .- view(t,j,:)
    r2 = view(t,k,:) .- view(t,j,:)
    acos.(1 .- colwise(CosineDist(),r1,r2)) .* 180 ./ π
end

"""
    dihedral(t::AbstractTrajectory, i,j,k,m) -> Vector{Float64}
"""
function dihedral(t::AbstractTrajectory, i,j,k,m)
    out = zeros(length(t))
    @warn "dihedral does not understand periodicity"
    for n in 1:length(t)
        b1 = view(t,j,n) .- view(t,i,n)
        b2 = view(t,k,n) .- view(t,j,n)
        b3 = view(t,m,n) .- view(t,k,n)
        out[n]=atand((b1×b2)×(b2×b3)⋅b2 / norm(b2), (b1×b2)⋅(b2×b3))
    end
    return out
end



"""
    angletoframe(t::AbstractTrajectory, i, j, k, m; frame=:)

Computes angle betweem two vectors. First from atoms i to j second from k to m.
Only given frame is calculated for first given vector.
If frame=: (default) then calculate for whole trajectory.
"""
function angletoframe(t::AbstractTrajectory, i, j, k, m; frame=:)
    r1 = view(t,j,frame) .- view(t,i,frame)
    r2 = view(t,m,:) .- view(t,k,:)
    acosd.(1 .- colwise(CosineDist(),r1,r2))
end




"""
    sphericalview(t::AbstractTrajectory, i, j) -> Dict

Transfroms vector from atom i to j to spherical coordinates
"""
function sphericalview(t::AbstractTrajectory, i, j)
    function rdis(r)
         rr = sqrt.(sum(r.^2,dims=1)) #TODO fix for periodic
         reshape(rr, (length(rr),1))
    end
    r = view(t,j,:) .- view(t,i,:)
    rr = rdis(r)
    ϕ = atan.(view(r,2,:), view(r,1,:))
    θ = acos.( view(r,3,:) ./ rr  )
    return Dict("r"=>rr, "θ"=>θ.*180 ./ π, "ϕ"=>ϕ.*180 ./ π)
end










"""
    compute_rdf(t, ur1, ur2; mindis=undef, maxdis=9.0, nbins=100) -> Dict

Calculates radial distribution function
"""
function compute_rdf(
    t::Trajectory{TA,TC} where {TA<:AbstractAtomNames, TC<:AbstractOrthorombicCell},
    ur1::AbstractUnitRange,
    ur2::AbstractUnitRange;
    mindis=undef,
    maxdis=9.0,
    nbins=100
    )

    dis = distances(t, ur1, ur2)
    di = DiscretizeUniformWidth(nbins)

    # Collect data that are in defined range [mindis, maxdis]
    data = Dict{Int, Vector{Float64}}()
    for (i,j) in enumerate(ur1)
        tmp = @view dis[i,:,:]
        if mindis != undef
            q1 = tmp .< maxdis
            q2 = tmp .> mindis
            push!(data, j=> tmp[ q1 .& q2 ]  )
        else
            push!(data, j=> tmp[ tmp .< maxdis ]  )
        end
    end

    #Discretize data to form histograms
    ρ = length(ur2)*length(t) / cellvolume(t)  # Average density for ur2  in trajectry
    edges = Dict()
    counts = Dict()
    radius = Dict()
    for (key,val) in data
        if mindis != undef
            e = Array(LinRange(mindis, maxdis, nbins+1))
        else
            e = binedges(di, val)
        end
        r = 0.5 .* (e[1:end-1] .+ e[2:end])
        dr = diff(e)
        disc = LinearDiscretizer(e)
        count = get_discretization_counts(disc, val) ./ (4*π*ρ .* r.^2 .* dr)
        push!(edges, key=> e)
        push!(counts, key=> count)
        push!(radius, key=>r)
    end
    return Dict("r"=>radius, "rdf"=>counts, "edges"=>edges)
end


function compute_rdf(fname::AbstractString, ur1::AbstractUnitRange,
                     ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100)
    try
        t = read_trajectory(fname)
        return compute_rdf(t, ur1, ur2, mindis=mindis, maxdis=maxdis, nbins=nbins)
    catch
        @warn "file $fname failed"
    end
        return Dict()
end


"""
    compute_rdf(ur1, ur2, fnames...; mindis=0.0, maxdis=9.0, nbins=100) -> Dict()

Calculates radial distribution function from given input and adds up results
"""
function compute_rdf(ur1::AbstractUnitRange, ur2::AbstractUnitRange,
                        fnames...; mindis=0.0, maxdis=9.0, nbins=100)
    dtmp = pmap(
        x -> compute_rdf(x, ur1, ur2; mindis=mindis, maxdis=maxdis, nbins=nbins),
        fnames
    )

    # Sum up results
    data = []
    for x in dtmp
        if length(keys(x)) == 0
            continue
        end
        push!(data,x)
    end
    rdf = Dict()
    if length(data) >0
        for (k,v) in data[1]["rdf"]
            push!(rdf, k => deepcopy(v))
        end
    end
    if length(data) > 1
        for x in data[2:end]
            for (k,v) in x["rdf"]
                rdf[k] .+= v
            end
        end
    end
    for k in keys(rdf)
        rdf[k] ./= length(data)
    end
    rk = collect(keys(data[1]["r"]))[1]
    return Dict("r"=>data[1]["r"][rk] , "rdf" => rdf)
end
