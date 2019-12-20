module MDTrajAnalyse

export AbstractTrajectory,
       AbstractTrajectoryWithNames,
       Trajectory,
       TrajectoryWithNames,
       PeriodicCellTrajectory,
       PeriodicConstCellTrajectory,
       natoms,
       dihedral,
       distances,
       distances!,
       cellvolume,
       angletoframe,
       sphericalview,
       compute_rdf,

       rdf_from_files,
       read_pdb,
       read_trajectory,
       read_xyz,

       AbstractUnitCell,
       AbstractOrthorombicCell,
       CubicCell,
       NonPeriodic,
       OrthorombicCell,
       TriclinicCell,
       celldiag,
       cellmatrix,
       parsecell,
       volume

include("cell.jl")
include("trajectory.jl")
include("fileaccess.jl")

using .cell
using .trajectory
using .fileaccess


end # module
