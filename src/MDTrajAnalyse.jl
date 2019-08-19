module MDTrajAnalyse

export AbstractTrajectory,
       AbstractTrajectoryWithNames,
       Trajectory,
       TrajectoryWithNames,
       PeriodicCellTrajectory,
       natoms,
       distances,
       distances!,
       cellvolume,
       angletoframe,
       sphericalview,
       compute_rdf,

       read_xyz,
       read_pdb,
       parallel_rdf,
       parallel_rdf_from_files


include("trajectory.jl")
include("fileaccess.jl")

using .trajectory
using .fileaccess

greet() = print("Hello World!")

end # module
