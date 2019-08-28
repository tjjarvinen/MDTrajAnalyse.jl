module cell

using LinearAlgebra


export AbstractUnitCell,
       AbstractOrthorombicCell,
       CubicCell,
       OrthorombicCell,
       TriclinicCell,
       celldiag,
       cellmatrix,
       volume


abstract type AbstractUnitCell end
abstract type AbstractOrthorombicCell <: AbstractUnitCell end

struct CubicCell <: AbstractOrthorombicCell
    abc::Float64
end

struct OrthorombicCell <: AbstractOrthorombicCell
    abc::Vector{Float64}
    OrthorombicCell(a::Number,b::Number,c::Number) = new([a,b,c])
    OrthorombicCell(a::Number) = new([a,a,a])
    function OrthorombicCell(a::AbstractVector)
        @assert length(a) == 3
        new(a)
    end
end

struct TriclinicCell <: AbstractUnitCell
    abc::Matrix{Float64}
    function TriclinicCell(abc::AbstractMatrix)
        @assert size(abc) == 3
        new(abc)
    end
    function OrthorombicCell(a::AbstractVector)
        @assert length(a) == 3
        new(Diagonal(a))
    end
end

Base.show(io::IO, c::CubicCell) = print(io, "CubicCell a=",c.abc)
Base.show(io::IO, c::OrthorombicCell) = print(io, "OrthorombicCell")
Base.show(io::IO, c::TriclinicCell) = print(io, "TriclinicCell")

celldiag(c::CubicCell) = [c.abc, c.abc, c.abc]
celldiag(c::OrthorombicCell) = c.abc
celldiag(c::TriclinicCell) = diag(c.abc)

cellmatrix(c::CubicCell) = Diagonal(celldiag(c))
cellmatrix(c::OrthorombicCell) = Diagonal(c.abc)
cellmatrix(c::TriclinicCell) = c.abc


volume(c::CubicCell) = c.abc^3
volume(c::OrthorombicCell) = c.abc[1]*c.abc[2]*c.abc[3]
volume(c::TriclinicCell) = abs(det(c.abc))

end #module
