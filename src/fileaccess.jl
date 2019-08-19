module fileaccess

using ..trajectory
using Distributed


export read_xyz,
       read_pdb,
       parallel_rdf,
       parallel_rdf_from_files

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
        elseif occursin("ATOM", lines[i]) || occursin("HETATM", lines[i])
            push!(atoms, lines[i][77:78])
            push!(xyz, parse(Float64, lines[i][31:38]))
            push!(xyz, parse(Float64, lines[i][39:46]))
            push!(xyz, parse(Float64, lines[i][47:54]))
        elseif occursin("END", lines[i])
            i+=1
            break
        end
        i+=1
    end

    for line in lines[i:end]
        if occursin("ATOM", line) || occursin("HETATM", lines[i])
            push!(xyz, parse(Float64, line[31:38]))
            push!(xyz, parse(Float64, line[39:46]))
            push!(xyz, parse(Float64, line[47:54]))
        elseif occursin("CRYST1", line)
            terms = split(line)
            append!(crystal,parse.(Float64, terms[2:4]))
        end
    end
    natoms = length(atoms)

    #This is a hack!!!
    nframes = Int( floor(length(xyz)/(3*natoms)))
    cell = zeros(3,3,nframes)
    if length(crystal)*3 == nframes
        cell[1:9:end] = crystal[1:3:end]
        cell[5:9:end] = crystal[2:3:end]
        cell[9:9:end] = crystal[3:3:end]
    elseif length(crystal)*3 < nframes
        #TODO fix this cheating !!!
        cell[1:9:3*length(crystal)] = crystal[1:3:end]
        cell[5:9:3*length(crystal)] = crystal[2:3:end]
        cell[9:9:3*length(crystal)] = crystal[3:3:end]
    else
        cell[1:9:end] = crystal[1:3:3*nframes]
        cell[5:9:end] = crystal[2:3:3*nframes]
        cell[9:9:end] = crystal[3:3:3*nframes]
    end
    #Hack again!!!
    out = xyz[1:3*natoms*nframes]
    return PeriodicCellTrajectory( reshape(out ,3, natoms, nframes), cell )
end




function rdf_from_file(fname, ur1::AbstractUnitRange,
                     ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100)
    try
        t = read_pdb(fname)
        return compute_rdf(t, ur1, ur2, mindis=mindis, maxdis=maxdis, nbins=nbins)
    catch
        @warn "file $fname failed"
    end
        return Dict()
end

function parallel_rdf_from_files(fnames::AbstractVector{String}, ur1::AbstractUnitRange,
                     ur2::AbstractUnitRange; mindis=undef, maxdis=9.0, nbins=100)
    dtmp = pmap( x -> rdf_from_file(x, ur1, ur2, mindis=mindis, maxdis=maxdis, nbins=nbins)   , fnames)
    data = []
    for x in dtmp
        if length(keys(x)) == 2
            push!(data,x)
        end
    end
    rdf = Dict()
    for (k,v) in data[1]["rdf"]
        push!(rdf, k => deepcopy(v))
    end
    for x in data[2:end]
        for (k,v) in x["rdf"]
            rdf[k] .+= v
        end
    end
    for k in keys(rdf)
        rdf[k] ./= length(data)
    end
    rk = collect(keys(data[1]["r"]))[1]
    return Dict("r"=>data[1]["r"][rk] , "rdf" => rdf)
end

end  # module fileaccess
