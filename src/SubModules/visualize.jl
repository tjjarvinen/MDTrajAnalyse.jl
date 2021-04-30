module Visualize

export visualize,
       external_visualize

using AbstractPlotting
using ..MDTrajAnalyse

const colours = Dict(
    "C"=>:black,
    " C"=>:black,
    "O"=>:red,
    " O"=>:red,
    "H"=>:white,
    " H"=>:white,
    "Ar"=>:orange,
    "N"=>:blue,
    " N"=>:blue
    )

const sizes = Dict(
    "C"=>0.9,
    " C"=>0.9,
    "O"=>0.9,
    " O"=>0.9,
    "H"=>0.6,
    " H"=>0.6,
    "Ar"=>1.0,
    "N"=>0.9,
    " N"=>0.9
    )

function visualize(traj::AbstractTrajectory, frame)
    acol = map( x -> colours[x], traj.names)
    asizes = map( x -> sizes[x], traj.names)
    sc = meshscatter(traj[frame]; markersize=asizes,color=acol)
    return sc
end


function visualize(traj::AbstractTrajectory; step=100)
    acol = map( x -> colours[x], traj.names)
    asizes = map( x -> sizes[x], traj.names)
    fig = Figure()
    ax = Axis3(fig[1,1])
    i = Slider(fig[2,1]; range=1:step:length(traj), startvalue=1)
    x = lift( j -> traj[j][1,:], i.value)
    y = lift( j -> traj[j][2,:], i.value)
    z = lift( j -> traj[j][3,:], i.value)
    meshscatter!(ax,x,y,z; markersize=asizes, color=acol)
    return fig
end


function external_visualize(traj::AbstractTrajectory;
    cmd="vlc", stdout=devnull, stderr=devnull, command="vmd"
    )
    fname= tempname()*".xyz"
    write_xyz(fname, traj)
    run(pipeline(`$(command) $(fname)`, stdout=stdout, stderr=stderr))
    rm(fname)
end



end # module
