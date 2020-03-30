module MDTrajAnalyse

export AbstractTrajectory,
       Trajectory,
       natoms,
       getatoms,
       dihedral,
       distances,
       distances!,
       cellvolume,
       angletoframe,
       sphericalview,
       subtrajectory,
       compute_rdf,
       get_close_atoms,

       AbstractAtomNames,
       AtomNames,
       NoAtomName,

       read_pdb,
       read_trajectory,
       read_xyz,
       write_xyz,

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
include("properties.jl")
include("SubModules/visualize.jl")

end # module
