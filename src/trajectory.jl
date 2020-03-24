

using Discretizers
using Distances
using LinearAlgebra

abstract type AbstractAtomNames end

abstract type AbstractTrajectory end


struct NoAtomName <: AbstractAtomNames end

struct AtomNames <: AbstractAtomNames
    names::Vector{String}
end

Base.getindex(a::AtomNames,i) = a.names[i]


"""
    Trajectory{TA<:AbstractAtomNames,TC<:AbstractUnitCell} <: AbstractTrajectory

Holds molecular dynamics trajectory

# Fields
- xyz::Array{Float64,3}  : Coordinates for the atoms
- names::TA              : Names of the atoms
- cell::TC               : Unitcell information
"""
struct Trajectory{TA<:AbstractAtomNames,TC<:AbstractUnitCell} <: AbstractTrajectory
    xyz::Array{Float64,3}
    names::TA
    cell::TC
    function Trajectory(xyz::AbstractArray{<:Real,3};
         names::AbstractAtomNames=NoAtomName(), cell::AbstractUnitCell=NonPeriodic())
        if size(xyz,1) != 3
            throw(DimensionMismatch("Trajectory - xyz has wrong dimensions"))
        end
        new{typeof(names), typeof(cell)}(xyz, names, cell)
    end
end

function Trajectory(xyz::AbstractArray{<:Real,3}, names::AbstractVector{String};
    cell=NonPeriodic()
    )
    return Trajectory(xyz; names=AtomNames(names), cell=cell)
end

function Trajectory(
    xyz::AbstractArray{<:Real,3},
    names::AbstractVector{String},
    cell::AbstractVector{T}
    ) where T <: AbstractUnitCell
    if all(x-> cell[1]==x,  cell)
        return Trajectory(xyz; names=AtomNames(names), cell=cell[1])
    end
    return Trajectory(xyz; names=AtomNames(names), cell=VariableCell(cell))
end


natoms(t::AbstractTrajectory) = size(t.xyz,2)

Base.length(t::AbstractTrajectory) = size(t.xyz,3)
function Base.show(io::IO, t::AbstractTrajectory)
    print(io, "Trajectory of ", length(t), " steps and ", natoms(t), " atoms" )
end

Base.getindex(t::AbstractTrajectory, frame) = t.xyz[:,:,frame]
Base.getindex(t::AbstractTrajectory, atom, frame) = t.xyz[:,atom,frame]
Base.lastindex(t::AbstractTrajectory) = length(t)
Base.firstindex(t::AbstractTrajectory) = 1

Base.view(t::AbstractTrajectory, frame) = view(t.xyz,:,:,frame)
Base.view(t::AbstractTrajectory, atom, frame) = view(t.xyz,:,atom,frame)

Base.setindex!(t::AbstractTrajectory, X, frame) = t.xyz[:,:,frame] = X
Base.setindex!(t::AbstractTrajectory, X, atom, frame) = t.xyz[:,atom,frame] = X

Base.iterate(t::AbstractTrajectory, state=1) = state > length(t) ? nothing : (view(t,state), state+1)

function Base.eltype(::Type{Trajectory{TA,TC}}) where {TA<:AbstractAtomNames, TC<:AbtractPeriodicCell}
    return Array{Float64,2}
end

function celldiag(t::Trajectory{TA,TC}, i) where {TA<:AbstractAtomNames, TC<:AbtractPeriodicCell}
    return celldiag(t)
end

function celldiag(t::Trajectory{TA,VariableCell}, i) where TA<:AbstractAtomNames
    return celldiag(t.cell[i])
end

function celldiag(t::Trajectory{TA,TC}) where {TA<:AbstractAtomNames, TC<:AbtractPeriodicCell}
    return celldiag(t.cell)
end


function distances(t::AbstractTrajectory, ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    out = Array{Float64}(undef, length(ur1), length(ur2), length(t))
    distances!(out,t,ur1,ur2)
    return out
end

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
    compute_rdf(t::AbstactPeriodicCellTrajectory, ur1::AbstractUnitRange,
                     ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100) -> Dict

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
