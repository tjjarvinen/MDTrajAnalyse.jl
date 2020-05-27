
using Distributed


"""
    read_trajectory(fname::AbstractString) -> Trajectory

Reads trajectory from given filename.
Type of file determined based on sufix.

## List of suported file formats
- XYZ
- PDB
"""
function read_trajectory(fname::AbstractString)
    sufix=split(fname, ".")[end]
    if sufix in ("pdb", "PDB")
        return read_pdb(fname)
    elseif sufix in ("xyz", "XYZ")
        return read_xyz(fname)
    elseif sufix in ("dcd", "DCD")
        return read_dcd(fname)
    else
        error("Trajectory file format not recognized")
    end
end


"""
    read_trajectory(fnames::AbstractString...) -> Vector{Trajectory}

Reads given files using threading
"""
function read_trajectory(fnames::AbstractString...)
    return fetch.( [Threads.@spawn(read_trajectory(x)) for x in fnames] )
end

"""
    read_dcd(fname::AbstractString) -> Trajectory

Reads DCD trajectory file format file.

Format definition is based on CP2K.
[Reference](https://github.com/cp2k/cp2k/blob/master/src/motion/dumpdcd.F)
"""
function read_dcd(fname::AbstractString)
    function read_cell(io::IO)
        a = read(io, Float64)
        γ = read(io, Float64)
        b = read(io, Float64)
        β = read(io, Float64)
        α = read(io, Float64)
        c = read(io, Float64)
        return parsecell(a,b,c,α,β,γ)
    end
    x = Float32[]
    y = Float32[]
    z = Float32[]
    natom_dcd = zero(Int32)
    cell = Vector{AbstractUnitCell}()
    frame = 0
    open(fname,"r") do io
        read(io,4) # Skip empty space

        id_dcd = String(read(io,4))
        istep = read(io, Int)
        iskip = read(io, Int32)

        read(io,24) # Skip empty space

        dt = read(io, Float32)
        have_unit_cell = read(io, Int32)

        read(io, 44) # Skip empty space

        nremark = read(io, Int32)
        remark_dcd = String(read(io, 80))
        remark_dcd2 = String(read(io, 80))

        read(io, 8) # Skip empty space

        natom_dcd = read(io, Int32)

        @debug "1st phase of DCD-file" natom_dcd, dt, id_dcd, istep, iskip, have_unit_cell, nremark

        while true
            try
                read(io, 8)
                push!(cell, read_cell(io))
                read(io, 8)
                append!(x , [read(io, Float32) for _ in 1:natom_dcd])
                read(io, 8)
                append!(y, [read(io, Float32) for _ in 1:natom_dcd])
                read(io, 8)
                append!(z, [read(io, Float32) for _ in 1:natom_dcd])
                frame += 1
            catch e
                if isa(e, EOFError)
                    break
                else
                    throw(e)
                end
            end
        end
    end
    @debug "Frame = $frame, natoms = $natom_dcd, length = $(length(x))"
    @assert length(x) == length(y) == length(z) "x-, y- or z coordinate lenghts dont match"
    @debug "array" Array(hcat(x,y,z)')
    return Trajectory( reshape(Array(hcat(x,y,z)'), 3, natom_dcd, frame), cell)
end


read_xyz(fname::AbstractString) = open(read_xyz, fname, "r")


function read_xyz(io::IO)
    lines = Vector{String}()
    atoms = Vector{String}()
    xyz = Vector{Float64}()

    na = parse(Int,readline(io))  # Number of atoms
    line = readline(io)   # Comment line
    for i in 1:na
        line = readline(io)
        cont = split(line)
        push!(atoms,cont[1])
        append!(xyz, parse.(Float64, cont[2:4]))
    end
    i = 1
    for line in eachline(io)
        cont = split(line)
        if i == na+2
            i = 1
            continue
        elseif i > 2
            length(cont) == 4 && append!(xyz, parse.(Float64, cont[2:4]))
        end
        i += 1
    end
    return Trajectory( reshape(xyz, 3, na, Int(length(xyz)/(3*na))), atoms)
end


function write_xyz(io::IO, traj::Trajectory{AtomNames, TC}) where {TC<:AbstractUnitCell}
    for x in traj
        println(io, "  ", length(traj.names),"\n")
        for (n, r) in zip(traj.names, eachcol(x))
            println(io, n, "   ", r[1], "  ", r[2], "  ", r[3])
        end
    end
    return nothing
end

function write_xyz(fname::AbstractString, traj::Trajectory)
    open(fname, "w") do file
        write_xyz(file, traj)
    end
    return nothing
end


function read_pdb(fname::AbstractString)
    xyz = Vector{Float64}()
    atoms = Vector{String}()
    crystal = Vector{AbstractUnitCell}()
    open(fname, "r") do file
        lineiterator = eachline(file)

        # Read data that need to read only once, like atom names
        for line in lineiterator
            if occursin(r"^CRYST1", line)
                terms = split(line)
                push!(crystal, parsecell(parse.(Float64, terms[2:7])...))
            elseif occursin(r"^ATOM", line) || occursin(r"^HETATM", line)
                push!(atoms, line[77:78])
                push!(xyz, parse(Float64, line[31:38]))
                push!(xyz, parse(Float64, line[39:46]))
                push!(xyz, parse(Float64, line[47:54]))
            elseif occursin(r"^END", line)
                break
            end
        end

        # Read bulk data
        for line in lineiterator
            if occursin(r"^ATOM", line) || occursin(r"^HETATM", line)
                push!(xyz, parse(Float64, line[31:38]))
                push!(xyz, parse(Float64, line[39:46]))
                push!(xyz, parse(Float64, line[47:54]))
            elseif occursin(r"^CRYST1", line)
                terms = split(line)
                push!(crystal, parsecell(parse.(Float64, terms[2:7])...))
            end
        end
    end

    @assert length(atoms)*3*length(crystal) == length(xyz) "Number of atoms and parsed coordinates do not match."

    return Trajectory(reshape(xyz ,3, length(atoms), length(crystal)), atoms, crystal)
end
