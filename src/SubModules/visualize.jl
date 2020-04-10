module Visualize

export visualize

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

sizes = Dict(
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

function visualize(traj::AbstractTrajectory, frame; rotationspeed=0.05)
    acol = map( x -> colours[x], traj.names)
    asizes = map( x -> sizes[x], traj.names)
    sc = meshscatter(traj[frame]'; markersize=asizes,color=acol)
    cameracontrols(sc).rotationspeed[] = rotationspeed
    return sc
end

function visualize(traj::AbstractTrajectory; step=100, rotationspeed=0.05)
    acol = map( x -> colours[x], traj.names)
    asizes = map( x -> sizes[x], traj.names)
    si, i  = textslider(1:step:length(traj), "frame"; start=1)
    x = lift( j -> traj[j][1,:], i)
    y = lift( j -> traj[j][2,:], i)
    z = lift( j -> traj[j][3,:], i)

    sc = Scene()

    meshscatter!(sc,x,y,z; markersize=asizes, color=acol)
    cameracontrols(sc).rotationspeed[] = rotationspeed
    final=hbox(si,sc)
end

end # module
