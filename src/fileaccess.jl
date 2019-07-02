module fileaccess

using ..trajectory

export read_xyz,
       read_pdb

function read_xyz(fname)
    lines = Vector{String}()
    open(fname, "r") do file
        lines = readlines(file)
    end

    # How many atoms
    natoms = parse(Int, lines[1])

    # How many clusters - use of floor allows extra empty lines at end
    nclusters = Int(floor(length(lines)/(natoms+2)))
    @debug "Type of nclusters $(typeof(nclusters))"

    xyz = zeros(Float64, 3, natoms)
    atoms = Vector{String}(undef, natoms)
    data = Vector{Float64}()

    for na in 1:natoms
        cont = split(lines[na+2])
        atoms[na] = cont[1]
        xyz[1,na] = parse(Float64, cont[2])
        xyz[2,na] = parse(Float64, cont[3])
        xyz[3,na] = parse(Float64, cont[4])
    end
    append!(data,xyz)

    for nc in 2:nclusters
        for na in 1:natoms
            cont = split(lines[(nc-1)*(natoms+2)+na+2])
            xyz[1,na] = parse(Float64, cont[2])
            xyz[2,na] = parse(Float64, cont[3])
            xyz[3,na] = parse(Float64, cont[4])
        end
        append!(data,xyz)
    end

    return TrajectoryWithNames( reshape(data, 3, natoms, Int(length(data)/(3*natoms))), atoms)
end

function read_pdb(fname)
    #NOTE only orthorombic, tetragonal and cubic crystals
    lines = Vector{String}()
    open(fname, "r") do file
        lines = readlines(file)
    end
    xyz = Vector{Float64}()
    atoms = Vector{String}()
    crystal = Vector{Float64}()
    i = 1

    while i<length(lines)
        if occursin("CRYST1", lines[i])
            terms = split(lines[i])
            append!(crystal,parse.(Float64, terms[2:4]))
        elseif occursin("ATOM", lines[i])
            terms = split(lines[i])
            push!(atoms, terms[end])
            append!(xyz, parse.(Float64, terms[4:6]) )
        elseif occursin("END", lines[i])
            i+=1
            break
        end
        i+=1
    end

    for line in lines[i:end]
        if occursin("ATOM", line)
            terms = split(line)
            append!(xyz, parse.(Float64, terms[4:6]) )
        elseif occursin("CRYST1", line)
            terms = split(line)
            append!(crystal,parse.(Float64, terms[2:4]))
        end
    end
    natoms = length(atoms)
    nframes = Int(length(xyz)/(3*natoms))
    cell = zeros(3,3,nframes)

    cell[1:9:end] = crystal[1:3:end]
    cell[5:9:end] = crystal[2:3:end]
    cell[9:9:end] = crystal[3:3:end]
    return PeriodicCellTrajectory( reshape(xyz ,3, natoms, Int(length(xyz)/(3*natoms))  ), cell )
end


end  # module fileaccess
