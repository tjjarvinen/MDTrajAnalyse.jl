module trajectory

using Distances, LinearAlgebra, Discretizers
using ..cell

export AbstractTrajectory,
       AbstractTrajectoryWithNames,
       Trajectory,
       TrajectoryWithNames,
       PeriodicCellTrajectory,
       natoms,
       distances,
       distances!,
       cellvolume,
       angletoframe,
       sphericalview,
       compute_rdf


abstract type AbstractTrajectory end
abstract type AbstractTrajectoryWithNames <: AbstractTrajectory end
abstract type AbstactPeriodicCellTrajectory{T} <: AbstractTrajectory  where T<:AbstractUnitCell  end


struct Trajectory <: AbstractTrajectory
    xyz::Array{Float64,3}
    function Trajectory(xyz::AbstractArray{<:Real,3})
        if size(xyz,1) != 3
            throw(DimensionMismatch("Trajectory - xyz has wrong dimensions"))
        end
        new(xyz)
    end
end


struct TrajectoryWithNames <: AbstractTrajectoryWithNames
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

struct PeriodicCellTrajectory{T} <: AbstactPeriodicCellTrajectory{T}
    xyz::Array{Float64,3}
    cell::Vector{T}
    function PeriodicCellTrajectory(xyz::AbstractArray{<:Real,3}, c::AbstractVector{T}) where T <: AbstractUnitCell
        if size(xyz,1) != 3
            throw(DimensionMismatch("PeriodicCellTrajectory - xyz has wrong dimensions"))
        end
        if length(c) != size(xyz, 3)
            throw(DimensionMismatch("PeriodicCellTrajectory - cell and xyz have different dimensions"))
        end
        new{T}(xyz, c)
    end
end

struct PeriodicConstCellTrajectory{T} <: AbstactPeriodicCellTrajectory{T}
    xyz::Array{Float64,3}
    cell::T
    function PeriodicCellTrajectory(xyz::AbstractArray{<:Real,3}, c::T)  where T <: AbstractUnitCell
        if size(xyz,1) != 3
            throw(DimensionMismatch("PeriodicCellTrajectory - xyz has wrong dimensions"))
        end
        new{T}(xyz, c)
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

cell.celldiag(t::AbstactPeriodicCellTrajectory, i) = celldiag(t.cell[i])
cell.celldiag(t::PeriodicConstCellTrajectory, i) = celldiag(t.cell)


function distances(t::AbstractTrajectory, ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    out = Array{Float64}(undef, length(ur1), length(ur2), length(t))
    distances!(out,t,ur1,ur2)
end

function distances!(r::AbstractArray{T,3} where T, t::AbstractTrajectory,
                     ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    L = Euclidean()
    for i in 1:length(t)
        pairwise!(view(r,:,:,i),L, view(t,ur1,i), view(t,ur2,i),dims=2)
    end
    return r
end

function distances!(r::AbstractArray{T,3} where T, t::AbstactPeriodicCellTrajectory,
                     ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    for i in 1:length(t)
        pairwise!(view(r,:,:,i),PeriodicEuclidean(celldiag(t,i)), view(t,ur1,i), view(t,ur2,i),dims=2)
    end
    return r
end

function distances!(r::AbstractArray{T,3} where T, t::AbstactPeriodicCellTrajectory{TriclinicCell},
                     ur1::AbstractUnitRange, ur2::AbstractUnitRange)
    @error("Not implemented")
end


function distances(t::AbstractTrajectory, i1::Integer, i2::Integer)
    colwise(Euclidean(), view(t,i1,:), view(t,i2,:))
end

function distances!(r::AbstractVector, t::AbstractTrajectory, i1::Integer, i2::Integer)
    colwise!(r, Euclidean(), view(t,i1,:), view(t,i2,:))
end

function distances(t::PeriodicConstCellTrajectory, i1::Integer, i2::Integer)
    M = PeriodicEuclidean(celldiag(t,1))
    colwise(M, view(t,i1,:), view(t,i2,:))
end

function distances(t::AbstactPeriodicCellTrajectory, i1::Integer, i2::Integer)
    [ evaluate(PeriodicEuclidean(celldiag(t,j)), view(t,i1,j), view(t,i2,j)) for j in 1:length(t) ]
end

function distances(t::AbstactPeriodicCellTrajectory{TriclinicCell}, i1::Integer, i2::Integer)
    @error("Not implemented")
end

function distances!(r::AbstractVector, t::PeriodicConstCellTrajectory, i1::Integer, i2::Integer)
    M = PeriodicEuclidean(celldiag(t,1))
    colwise!(r, M, view(t,i1,:), view(t,i2,:))
end

function distances!(r::AbstractVector, t::AbstactPeriodicCellTrajectory, i1::Integer, i2::Integer)
    r = [ evaluate(PeriodicEuclidean(celldiag(t,j)), view(t,i1,j), view(t,i2,j)) for j in 1:length(t) ]
end

function distances!(r::AbstractVector, t::AbstactPeriodicCellTrajectory{TriclinicCell}, i1::Integer, i2::Integer)
    @error("Not implemented")
end

cellvolume(t::PeriodicConstCellTrajectory) = volume(t.cell)
cellvolume(t::PeriodicConstCellTrajectory,i) = cellvolume(t)

cellvolume(t::PeriodicCellTrajectory) = sum(volume.(t.cell)) / length(t)
cellvolume(t::PeriodicCellTrajectory, i) = volume(t.cell[i])


"""
    angle(t::AbstractTrajectory, i, j, k)

Computes angle for atoms i, j, k in degrees
"""
function Base.angle(t::AbstractTrajectory, i, j, k)
    r1 = view(t,i,:) .- view(t,j,:)
    r2 = view(t,k,:) .- view(t,j,:)
    acos.(1 .- colwise(CosineDist(),r1,r2)) .* 180 ./ π
end

"""
    angletoframe(t::AbstractTrajectory, i, j, k; frame=:)

Computes angle two vectors. First from atoms i to j second from k to m.
Only given frame is calculated for first given vector.
If frame=: (default) then calculate for whole trajectory.
"""
function angletoframe(t::AbstractTrajectory, i, j, k, m; frame=:)
    r1 = view(t,j,frame) .- view(t,i,frame)
    r2 = view(t,m,:) .- view(t,k,:)
    acos.(1 .- colwise(CosineDist(),r1,r2)) .* 180 ./ π
end




"""
    sphericalview(t::AbstractTrajectory, i, j)

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
                     ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100)

Calculates radial distribution function
"""
function compute_rdf(t::AbstactPeriodicCellTrajectory, ur1::AbstractUnitRange,
                     ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100)
    dis = distances(t, ur1, ur2)
    vd = cellvolume(t)  # this is mean volume which is ok
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

    ρ = length(ur2)*length(t) / cellvolume(t)
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
