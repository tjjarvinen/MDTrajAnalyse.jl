

using LinearAlgebra
using StaticArrays



abstract type AbstractUnitCell end
abstract type AbtractPeriodicCell <: AbstractUnitCell end
abstract type AbstractOrthorombicCell <: AbtractPeriodicCell end

struct NonPeriodic <: AbstractUnitCell end

struct VariableCell <: AbtractPeriodicCell
    cell::Vector{AbtractPeriodicCell}
end

struct CubicCell <: AbstractOrthorombicCell
    a::Float64
end

struct OrthorombicCell <: AbstractOrthorombicCell
    abc::SVector{3,Float64}
    OrthorombicCell(a::Real,b::Real,c::Real) = new([a,b,c])
    OrthorombicCell(a::Real) = new([a,a,a])
    function OrthorombicCell(a::AbstractVector)
        @assert length(a) == 3
        new(a)
    end
end

struct TriclinicCell <: AbtractPeriodicCell
    abc::SMatrix{3,3,Float64}
    function TriclinicCell(abc::AbstractMatrix)
        @assert size(abc) == (3,3)
        new(abc)
    end
    function TriclinicCell(a::AbstractVector)
        @assert length(a) == 3
        new(Diagonal(a))
    end
    function TriclinicCell(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
        tmp = zeros(3,3)
        tmp[1,1] = sind(γ)*a
        tmp[2,1] = cosd(γ)*a
        tmp[2,2] = b
        tmp[2,3] = cosd(α)*c
        tmp[1,3] = c*(cosd(β)-cosd(γ)*cosd(α))/sind(γ)
        tmp[3,3] = sqrt(c^2-tmp[1,3]^2-tmp[2,3]^2)
        new(tmp)
    end
end

parsecell(a::Number) = CubicCell(a)
function parsecell(a::Real, b::Real, c::Real)
    a == b == c && return CubicCell(a)
    return OrthorombicCell(a,b,c)
end

function parsecell(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    if α == β == γ == 90
        a == b == c && return CubicCell(a)
        return OrthorombicCell(a,b,c)
    else
        return TriclinicCell(a,b,c,α,β,γ)
    end
end


Base.show(io::IO, c::CubicCell) = print(io, "CubicCell a=",c.a)
Base.show(io::IO, c::OrthorombicCell) = print(io, "OrthorombicCell")
Base.show(io::IO, c::TriclinicCell) = print(io, "TriclinicCell")
Base.show(io::IO, c::VariableCell) = print(io, "Variable cell of length ", length(c))

Base.length(c::VariableCell) = length(c.cell)
Base.getindex(c::VariableCell, i) = c.cell[i]
Base.lastindex(c::VariableCell) = length(c)
Base.firstindex(c::VariableCell) = 1

Base.iterate(c::VariableCell, state=1) = state > length(c) ? nothing : (c[state], state+1)

celldiag(c::AbtractPeriodicCell) = diag(cellmatrix(c))

cellmatrix(c::CubicCell) = Diagonal([c.a, c.a, c.a])
cellmatrix(c::OrthorombicCell) = Diagonal(c.abc)
cellmatrix(c::TriclinicCell) = c.abc
cellmatrix(c::VariableCell) = [ cellmatrix(x) for x in c]
cellvectorlengths(c::AbtractPeriodicCell) = [norm(x) for x in eachcol(cellmatrix(c))]
cellvectorlengths(c::AbstractOrthorombicCell) = (Diagonal ∘ cellmatrix)(c)
cellvectorlengths(c::VariableCell) = [ cellvectorlengths(x) for x in c ]


volume(c::CubicCell) = c.a^3
volume(c::OrthorombicCell) = prod(c.abc)
volume(c::AbtractPeriodicCell) = (abs ∘ det ∘ cellmatrix)(c)
volume(c::VariableCell) = [ volume(x) for x in c ]
