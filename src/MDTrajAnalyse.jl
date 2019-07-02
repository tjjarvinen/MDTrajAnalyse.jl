module MDTrajAnalyse

export AbstractTrajectory,
       AbstractTrajectoryWithNames,
       Trajectory,
       TrajectoryWithNames,
       PeriodicCellTrajectory,
       natoms,
       distances,
       compute_rdf,

       read_xyz,
       read_pdb


include("trajectory.jl")
include("fileaccess.jl")

using .trajectory
using .fileaccess

greet() = print("Hello World!")

end # module
