
import Base.+
import Base.*

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

function Trajectory(
    xyz::AbstractArray{<:Real,3},
    names::AbstractVector{String},
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

Base.length(names::AtomNames) = length(names.names)
Base.getindex(names::AtomNames) = names.names[i]
Base.lastindex(names::AtomNames) = length(names)
Base.firstindex(names::AtomNames) = 1

Base.getindex(t::AbstractTrajectory, frame) = t.xyz[:,:,frame]
Base.getindex(t::AbstractTrajectory, atom, frame) = t.xyz[:,atom,frame]
Base.lastindex(t::AbstractTrajectory) = length(t)
Base.firstindex(t::AbstractTrajectory) = 1

Base.view(t::AbstractTrajectory, frame) = view(t.xyz,:,:,frame)
Base.view(t::AbstractTrajectory, atom, frame) = view(t.xyz,:,atom,frame)

Base.setindex!(t::AbstractTrajectory, X, frame) = t.xyz[:,:,frame] = X
Base.setindex!(t::AbstractTrajectory, X, atom, frame) = t.xyz[:,atom,frame] = X

Base.iterate(t::AbstractTrajectory, state=1) = state > length(t) ? nothing : (view(t,state), state+1)

Base.iterate(names::AtomNames, state=1) = state > length(names) ? nothing : (names[state],state+1)

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



function subtrajectory(t::AbstractTrajectory, i)
    return Trajectory(t[i]; names=t.names, cell=t.cell)
end


function subtrajectory(t::AbstractTrajectory, i::Integer)
    return Trajectory(reshape(t[i], (3,size(t[i], 2),1)); names=t.names, cell=t.cell)
end

function subtrajectory(
    t::Trajectory{T,VariableCell} where {T<:AbstractAtomNames},
    i::Integer
    )
    return Trajectory(reshape(t[i], (3,size(t[i], 2),1)); names=t.names, cell=t.cell[i])
end

function subtrajectory(
    t::Trajectory{NoAtomName,T} where {T<:AbstractUnitCell},
    i
    )
    return Trajectory(t[i]; names=t.names, cell=t.cell)
end

function subtrajectory(
    t::Trajectory{NoAtomName,VariableCell},
    i
    )
    return Trajectory(t[i]; names=t.names, cell=t.cell[i])
end

function subtrajectory(
    t::Trajectory{TA,VariableCell} where {TA<:AbstractAtomNames},
    i
    )
    return Trajectory(t[i]; names=t.names, cell=t.cell[i])
end

(+)(a1::AtomNames, a2::AtomNames) = AtomNames(vcat(a1.names, a2.names))

function (+)(t1::AbstractTrajectory, t2::AbstractTrajectory)
    @assert t1.names == t2.names
    @assert t1.cell == t2.cell
    return Trajectory(cat(t1.xyz, t2.xyz; dims=3); names=t1.names, cell=t1.cell)
end

function (+)(
    t1::Trajectory{TA,VariableCell},
    t2::Trajectory{TA,VariableCell}
    ) where {TA<:AbstractAtomNames}

    @assert t1.names == t2.names
    return Trajectory(cat(t1.xyz, t2.xyz; dims=2);
     names=t1.names+t2.names, cell=vcat(t1.cell, t2.cell)
     )
end

function (*)(t1::Trajectory{AtomNames,TC}, t2::Trajectory{AtomNames,TC}) where {TC<:AbstractUnitCell}
    @assert length(t1) == length(t2)
    @assert t1.cell == t2.cell
    return Trajectory(cat(t1.xyz, t2.xyz; dims=2); names=t1.names+t2.names, cell=t1.cell)
end

function (*)(t1::Trajectory{NoAtomName,TC}, t2::Trajectory{NoAtomName,TC}) where {TC<:AbstractUnitCell}
    @assert length(t1) == length(t2)
    @assert t1.cell == t2.cell
    return Trajectory(cat(t1.xyz, t2.xyz; dims=2); names=NoAtomName(), cell=t1.cell)
end

getatoms(t::AbstractTrajectory, i) = Trajectory(t.xyz[:,i,:], t.names[i], t.cell)

function getatoms(t::Trajectory{NoAtomName,T} where {T<:AbstractUnitCell}, i)
    return Trajectory(t.xyz[:,i,:]; names=NoAtomName(), cell=t.cell)
end
