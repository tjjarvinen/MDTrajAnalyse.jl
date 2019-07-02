module trajectory

using Distances, LinearAlgebra, Discretizers

export AbstractTrajectory,
       AbstractTrajectoryWithNames,
       Trajectory,
       TrajectoryWithNames,
       PeriodicCellTrajectory,
       natoms,
       distances,
       compute_rdf


abstract type AbstractTrajectory end
abstract type AbstractTrajectoryWithNames <: AbstractTrajectory end
abstract type AbstactPeriodicCellTrajectory <: AbstractTrajectory end


mutable struct Trajectory <: AbstractTrajectory
    xyz::Array{Float64,3}
    function Trajectory(xyz::AbstractArray{<:Real,3})
        if size(xyz,1) != 3
            throw(DimensionMismatch("Trajectory - xyz has wrong dimensions"))
        end
        new(xyz)
    end
end


mutable struct TrajectoryWithNames <: AbstractTrajectoryWithNames
    xyz::Array{Float64,3}
    anames::Vector{String}
    function TrajectoryWithNames(xyz::AbstractArray{<:Real,3}, anames::Vector{String})
        if size(xyz,1) != 3
            throw(DimensionMismatch("TrajectoryWithNames - xyz has wrong dimensions"))
        end
        if size(xyz,2) != length(anames)
            throw(DimensionMismatch("TrajectoryWithNames - xyz and anames have different ammount of atoms"))
        end
        new(xyz, anames)
    end
end

mutable struct PeriodicCellTrajectory <: AbstactPeriodicCellTrajectory
    xyz::Array{Float64,3}
    cell::Array{Float64,3}
    function PeriodicCellTrajectory(xyz::AbstractArray{<:Real,3}, cell::AbstractArray{<:Real,3})
        if size(xyz,1) != 3
            throw(DimensionMismatch("PeriodicCellTrajectory - xyz has wrong dimensions"))
        end
        if size(cell, 1) != 3 || size(cell, 2) != 3
            throw(DimensionMismatch("PeriodicCellTrajectory - cell has wrong dimensions"))
        end
        if size(cell,3) != size(xyz, 3)
            throw(DimensionMismatch("PeriodicCellTrajectory - cell and xyz have different dimensions"))
        end
        new(xyz,cell)
    end
end


natoms(t::AbstractTrajectory) = size(t.xyz,2)

Base.length(t::AbstractTrajectory) = size(t.xyz,3)
Base.show(io::IO, t::AbstractTrajectory) = print(io, "Trajectory of ",
            length(t), " steps and ", natoms(t), " atoms" )

Base.getindex(t::AbstractTrajectory, frame) = t.xyz[:,:,frame]
Base.getindex(t::AbstractTrajectory, atom, frame) = t.xyz[:,atom,frame]
Base.lastindex(t::AbstractTrajectory) = length(t)

Base.view(t::AbstractTrajectory, frame) = view(t.xyz,:,:,frame)
Base.view(t::AbstractTrajectory, atom, frame) = view(t.xyz,:,atom,frame)

Base.setindex!(t::AbstractTrajectory, X, frame) = t.xyz[:,:,frame] = X
Base.setindex!(t::AbstractTrajectory, X, atom, frame) = t.xyz[:,atom,frame] = X

function distances(t::AbstractTrajectory, ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    L = Euclidean()
    out = Array{Float64}(undef, length(ur1), length(ur2), length(t))
    for i in 1:length(t)
        out[:,:,i]=out,pairwise(L, view(t,i)[:,ur1], view(t,i)[:,ur2],dims=2)
    end
    return out
end

function distances(t::AbstactPeriodicCellTrajectory, ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    out = Array{Float64}(undef, length(ur1), length(ur2), length(t))
    for i in 1:length(t)
        out[:,:,i]=pairwise(PeriodicEuclidean(diag(t.cell[:,:,i])), view(t,i)[:,ur1], view(t,i)[:,ur2],dims=2)
    end
    return out
end

function compute_rdf(t::AbstactPeriodicCellTrajectory, ur1::AbstractUnitRange, ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100)
    dis = distances(t, ur1, ur2)
    vd = diag(t.cell[:,:,1])
    ρ = length(ur2) / (vd[1]*vd[2]*vd[3])
    di = DiscretizeUniformWidth(nbins)
    data = Dict()
    for i in ur1
        tmp = @view dis[i,:,:]
        if mindis != undef
            q1 = tmp .< maxdis
            q2 = tmp .> mindis
            push!(data, i=> tmp[ q1 .& q2 ]  )
        else
            push!(data, i=> tmp[ tmp .< maxdis ]  )
        end
    end

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
    return Dict("r"=>radius, "rdf"=>counts)
end

end  # module trajectory
