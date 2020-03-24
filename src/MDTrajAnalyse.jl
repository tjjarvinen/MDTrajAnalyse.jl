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

       AbstractAtomNames,
       AtomNames,
       NoAtomName,


       rdf_from_files,
       read_pdb,
       read_trajectory,
       read_xyz,

       AbstractUnitCell,
       AbtractPeriodicCell,
       AbstractOrthorombicCell,
       CubicCell,
       NonPeriodic,
       OrthorombicCell,
       TriclinicCell,
       VariableCell,
       celldiag,
       cellmatrix,
       parsecell,
       volume

include("cell.jl")
include("trajectory.jl")
include("fileaccess.jl")


end # module
