
using Distributed


"""
    read_trajectory(fname::AbstractString) -> Trajectory

Reads trajectory from given filename.
Type of file

## List of suported file formats
"""
function read_trajectory(fname::AbstractString)
    sufix=split(fname, ".")[end]
    if sufix in ["pdb", "PDB"]
        return read_pdb(fname)
    elseif sufix in ["xyz", "XYZ"]
        return read_xyz(fname)
    else
        error("Trajectory file format not recognized")
    end
end

function read_dcd(fname::AbstractString)
    #TODO
end


function read_xyz(fname::AbstractString)
    lines = Vector{String}()
    na = 0
    atoms = Vector{String}()
    xyz = Vector{Float64}()
    open(fname, "r") do file
        na = parse(Int,readline(file))  # Number of atoms
        line = readline(file)
        for i in 1:na
            line = readline(file)
            cont = split(line)
            push!(atoms,cont[1])
            append!(xyz, parse.(Float64, cont[2:4]))
        end
        for line in eachline(file)
            cont = split(line)
            length(cont) == 4 && append!(xyz, parse.(Float64, cont[2:4]))
        end
    end
    return Trajectory( reshape(xyz, 3, na, Int(length(xyz)/(3*na))), atoms)
end


function write_xyz(fname::AbstractString, traj::Trajectory)

    #TODO
end


function read_pdb(fname::AbstractString)
    xyz = Vector{Float64}()
    atoms = Vector{String}()
    crystal = Vector{AbstractUnitCell}()
    open(fname, "r") do file
        lineiterator = eachline(file)

        for line in lineiterator
            if occursin("CRYST1", line)
                terms = split(line)
                push!(crystal, parsecell(parse.(Float64, terms[2:7])...))
            elseif occursin("ATOM", line) || occursin("HETATM", line)
                push!(atoms, line[77:78])
                push!(xyz, parse(Float64, line[31:38]))
                push!(xyz, parse(Float64, line[39:46]))
                push!(xyz, parse(Float64, line[47:54]))
            elseif occursin("END", line)
                break
            end
        end

        for line in lineiterator
            if occursin("ATOM", line) || occursin("HETATM", line)
                push!(xyz, parse(Float64, line[31:38]))
                push!(xyz, parse(Float64, line[39:46]))
                push!(xyz, parse(Float64, line[47:54]))
            elseif occursin("CRYST1", line)
                terms = split(line)
                push!(crystal, parsecell(parse.(Float64, terms[2:7])...))
            end
        end
    end

    @assert length(atoms)*3*length(crystal) == length(xyz) "Number of atoms and parsed coordinates do not match."

    return Trajectory(reshape(xyz ,3, length(atoms), length(crystal)), atoms, crystal)
end




function _rdf_from_file(fname::AbstractString, ur1::AbstractUnitRange,
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
    rdf_from_files(ur1::AbstractUnitRange, ur2::AbstractUnitRange, fnames...;
                     mindis=0.0, maxdis=9.0, nbins=100)
"""
function rdf_from_files(ur1::AbstractUnitRange, ur2::AbstractUnitRange,
                        fnames...; mindis=0.0, maxdis=9.0, nbins=100)
    dtmp = pmap( x -> _rdf_from_file(x, ur1, ur2, mindis=mindis, maxdis=maxdis, nbins=nbins) , fnames)

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
