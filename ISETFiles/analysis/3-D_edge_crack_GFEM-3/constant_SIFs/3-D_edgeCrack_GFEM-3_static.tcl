# ************************** 3-D_edgeCrack_GFEM-3_static.tcl **************************

set crackFile "edge_cracked_bar_crack_surface"
set resultsFolder "res_static"
file mkdir $resultsFolder

set domainFile "edge_cracked_bar_domain"

set isVisco false
set generateParaviewFiles true

# If the following is true, the code will solve all 3rd-scale local problems in parallel.
# Use with caution! This has not been extensively tested.
set solveLocalProbsInParallel false

# ************************** 3-D_edgeCrack_GFEM-3_static.tcl **************************

# Turn on creation of extract file
set env(DO_QA) ON

set num_threads 16
# Set number of threads for "Pardiso Solver"
set env(MKL_NUM_THREADS) $num_threads
# Set number of threads for "Apple Accelerate Framework BLAS and LAPACK"
set env(VECLIB_MAXIMUM_THREADS) $num_threads
# Set number of threads for assemble
set env(OMP_NUM_THREADS) $num_threads

# Create a linear analysis object
createAnalysis linear
analysis setParallelAssemble

# Use pFEM-GFEM approximation space:
# It gives better SIFs than GFEM but implementation is
# not as robust as GFEM for crack propagation at this time.
parSet useCompElPfemGfem true

# All non-polynomial enrichments are shifted by their nodal values:
# This improves conditioning.
parSet ShiftedSpBasisEnrichments true

# Read the model
readFile "$domainFile.grf" phfile

# These coordinates will define a bounding box
# for the second-scale, coarser local problem.
set x_min -0.01
set y_min -1.01
set z_min -1.51
set x_max  2.01
set y_max  1.01
set z_max  1.51

# Set quantities to be postprocessed
graphMesh appendName Solution
graphMesh appendName VonMises

# If wants to see input geomesh and BCs
set printInputGMesh true

if { $printInputGMesh == true } {
  
  if { $generateParaviewFiles == true } {
    
    graphMesh create "$resultsFolder/3-D_edgeCrack_geomesh.vtu" VTKStyle
    graphMesh postProcess resolution 0
    # Uncomment if want to quit
    # exit

  }

}

# Refine mesh globally
set nref 0
for { set iref 0 } { $iref < $nref } { incr iref } {
  
  refine all

}
meshReport

# Set resulting p-order:
# Going over p = 2 may cause conditioning problems.
set p_glob 2
set p_loc 2

# CPU time used for complete simulation
set totalTime 0

# Solve the Initial Global Problem (IGP)
puts "\n\n@@@@ Initial Global Problem (IGP) ***\n\n"
enrichApprox iso approxOrder $p_glob
assemble
solve

# Creates the GraphMesh
if { $generateParaviewFiles == true } {

  graphMesh create "$resultsFolder/PostProcessMesh_3-D_edgeCrack_IG.vtu" VTKStyle
  graphMesh postProcess resolution 0 sampleInside

}

# ************************** Setup for crack **************************
# The following parameters are related to the edge-crack problem:
# Remember to change this when running the pavement problem.

# Limit for length of crack front edges (defined by edges of facets along crack front):
# TIP: Around 1/20 of length of initial crack front
set initialCrackFrontLength 3.0
set initalCrackLength_a_o 1.0
set frontEdgeLengthLimit [ expr $initialCrackFrontLength / 61.0 ]

# To use crack surface remeshing (set to true or false):
# Remeshing can be done for the input crack surface, IF needed.
set useRemeshCrackSurface false
set remeshInputCrackSurface false

# Limit for length of edges of crack surface facets (enforced when remeshing):
# This regards the size of the crack surface facets NOT at the front.
# It is normally bigger than $frontEdgeLengthLimit
set mycrackSurfEdgeLength [ expr 2.0 * $frontEdgeLengthLimit ]

# Whether to preserve current crack front when remeshing:
# Needed only in mixed-mode problems with sharp turns that need to be preserved.
# If set to true, code below will preserve only the crack front at the first step.
set storeCrackFront false

# ************************** Setup for crack in the Global Problem **************************

crackMgr crackMgrID 1 create crackFile "$crackFile.crf"

# Post-process first crack file
if { $generateParaviewFiles == true } {
  
  crackMgr crackMgrID 1 postSurface "$resultsFolder/Crack_initial_surf.vtk"
  
}

# ************************** Parameters for remeshing the crack **************************

# To use crack surface remeshing (set to true or false)
crackMgr crackMgrID 1 setOptions useRemeshCrackSurface $useRemeshCrackSurface

# Limit for crack surface edge size when remeshing
crackMgr crackMgrID 1 setOptions crackSurfEdgeLength $mycrackSurfEdgeLength

# Whether to preserve current crack front when remeshing
crackMgr crackMgrID 1 setOptions storeCrackFront $storeCrackFront

# Front will be refined if edges are longer than this limit
crackMgr crackMgrID 1 setOptions frontEdgeLengthLimit $frontEdgeLengthLimit

# Remesh initial crf file:
# The use of "classified nodes" here is a workaround for remeshing of crack surface (CGAL) function.
# WARNING:
# I can use this because of classified nodes defined in .crf.
# Otherwise 2 facets will be missing.

if { $remeshInputCrackSurface == true } {

  crackMgr crackMgrID 1 remeshCrackSurface 
  if { $generateParaviewFiles == true } {
    
    crackMgr crackMgrID 1 postSurface "$resultsFolder/Crack_initial_surf_after_remesh.vtk"
  
  }

}

# Do NOT extract at crack front points on the domain boundary:
# This is usually required for surface breaking cracks since extraction 
# at boundary points is not always possible.
crackMgr crackMgrID 1 setOptions numFrontEndVertToSkip 1
crackMgr crackMgrID 1 setOptions useMLSForSkippedFrontEndVert true

# Apply Moving Least Square (MLS) approximation to extracted SIFs:
# This must be used if we use option $numFrontEndVertToSkip.
crackMgr crackMgrID 1 setOptions useMLSForExtraction true

# The following is performed only if we need the code to extract 1 SIF only!
# Create crackGrowthLaw object:
crackMgr crackMgrID 1 crGrowthPhysicsLaw create
crackMgr crackMgrID 1 crGrowthPhysicsLaw crGrowthType FatigueCrackGrowth
set delta_a_max_law 1.0
crackMgr crackMgrID 1 crGrowthPhysicsLaw crFrontIncrementType fixedDeltaA delta_a_max $delta_a_max_law

# Scaling law
crackMgr crackMgrID 1 crGrowthPhysicsLaw crFrontScalingLaw ParisLaw C 1.0 m 1.0 R 1.0

# The following option is used for SIMPLIFIED simulations only:
# If false, SIFs will not be constant along crack front.
crackMgr crackMgrID 1 crGrowthPhysicsLaw crFrontScalingLaw setConstFrontDispl true

# ************************** Options for VISCOELASTIC extraction **************************
# Set all visco options, method, load control, load type, etc.

if { $isVisco == true } {

  set loadFnType "freq"
  # t_peak = 0.5 * (load period)
  # This next t_peak corresponds, for instance, to a period equals to 1800s, or 30min
  set t_peak 900.0

  # Tells the crack manager the name of the viscoelastic material defined in the .grf file
  crackMgr crackMgrID 1 setOptions useViscoExtractionCP viscoMatName "isoVisco4CP"

  # Method used to compute inverse Laplace transform:
  # If viscoFourierTransf is false: use Zakian's method (NOT recommended).
  crackMgr crackMgrID 1 setOptions viscoFourierTransf true

  # Used at ViscoCorrespondPrinc::get_lambda_s() method
  # viscoLoadType { constant | ramp | creep | hat | freq | haversine }
  crackMgr crackMgrID 1 setOptions viscoLoadType $loadFnType

  # Set the t_peak used for frequency or haversine loads
  crackMgr crackMgrID 1 setOptions viscoPeakTime $t_peak

  # Whether it is force or displacement controlled:
  # This is needed to compute ERR(t).
  crackMgr crackMgrID 1 setOptions viscoForceControlled true

  # CAD: This is the default. Setting just for clarity.
  # NOTE:
  # LaplaceTransf is the recommended option. Use others at your own risk...
  # Used at ViscoCorrespondPrinc::get_TimeExtract() method
  # viscoIntegrType { LaplaceTransf | GaussIntegr | incremental }
  crackMgr crackMgrID 1 setOptions viscoIntegrType LaplaceTransf

}

# ************************** Finished setting up cracks at global scale **************************

set totalTime_cp 0.0
set totalTime_pa 0.0

set totalTime_sol 0.0

puts stdout "\n\n*********************************************"
puts stdout "\nStarting 2nd-Scale Local Problem 10"
puts stdout "\n*********************************************\n\n"

createLocalProblem probID 10 \
                   xyzMin $x_min $y_min $z_min xyzMax $x_max $y_max $z_max \
                   numLayers 1 \
                   userBC "spr_local"

crackMgr localProb 10 crackMgrID 11 create crackSurfaceFromCrackMgr 1

# This is not yet allowed:
# parSet localProb 10 ShiftedSpBasisEnrichments true

# Refine mesh globally
set nref 0
for { set iref 0 } { $iref < $nref } { incr iref } {
  
  puts stdout "Beginning Step = $iref"
  refine localProb 10 all

}

# ************************** Parameters for PROCESSING the crack 11 **************************

# Use branch-heaviside pair in elements fully cut by crack surface and in the geometrical enrichment zone:
# This improves accuracy of solution without adding too many dofs.
parSet localProb 10 useBranchHeavisideFnPairWithGeomEnrich true

# Setting to use planar cuts away from crack front
crackMgr localProb 10 crackMgrID 11 setOptions usePlanarCutsSignedDist true

# Setting branch function type
crackMgr localProb 10 crackMgrID 11 setOptions useBranchFn true
crackMgr localProb 10 crackMgrID 11 setOptions branchFunctionType curvedFrontOD6

# Enrich cylinder around crack front:
# Enrich all nodes within radius R of the front
# TIP: Radius = 4 * maxEdgeLen
# crackMgr localProb 10 crackMgrID 11 branchFnGeometricEnrich inCylinder radius 0.0013
#
# OR
#
# Enrich nodes in layers around the crack front:
crackMgr localProb 10 crackMgrID 11 branchFnGeometricEnrich numElemLayers 4

# Refine 3-D GFEM elements cut by crack surface
set nrefSurface 0
for { set iref 0 } { $iref < $nrefSurface } { incr iref } {
  
  crackMgr localProb 10 crackMgrID 11 refineMesh crackSurface
  
}

# Refine 3-D GFEM mesh around crack front:
# Try to refine this second-scale local problem as much as you can.
# The idea of this new GFEM-3 is that this local problem will provide a good displacement field to the global problem.
# NOTE: I am not refining this problem at all. This is only to test the methodology.
set maxEdgeElemLenReq [ expr 0.30 * $initalCrackLength_a_o ]
# crackMgr localProb 10 crackMgrID 11 refineMesh crackFronts maxEdgeLen $maxEdgeElemLenReq

# Post-process crack surface at this step:
# It is the initial crack.
if { $generateParaviewFiles == true } {
  
  crackMgr localProb 10 crackMgrID 11 postSurface "$resultsFolder/Crack_surf.vtk" VTKStyle

}

set initTime [ clock clicks -milliseconds ]

# The following is set to true.
# Then, node snapping can be used along with GFEM-gl.
# Use this with caution!
parSet localProb 10 solveInKamikazeMode true

# Snap to the crack surface nodes that are close to it:
# This improves conditioning of global matrix.
# This was NOT extensively tested yet for GFEM-gl:
crackMgr localProb 10 crackMgrID 11 setOptions snapNodesToSurfAndFront snappingTolSurf 0.025 snappingTolFront 0.025

# Process crack:
# Cut elements, select enrichments, etc.
crackMgr localProb 10 crackMgrID 11 process

# Get min edge length
set minEdgeLen [ parGet localProb 10 minEdgeLen ]
puts stdout "\n ** minEdgeLen (2nd-scale LP) = $minEdgeLen"

set finalTime [ clock clicks -milliseconds ]
set totalTime_cp [ expr $totalTime_cp + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for crack process command: $totalTime_cp s \n"

# Set polynomial order of the approximation
enrichApprox localProb 10 iso approxOrder $p_loc

# Enforce the search of local descendants while calculating global-local BCs
# since the enriched global solution is computed based on the other local problem
parSet localProb 10 enforceDescSearchGlobLocBC true

# Do NOT extract at crack front points on the domain boundary:
# This is usually required for surface breaking cracks since extraction 
# at boundary points is not always possible.
crackMgr localProb 10 crackMgrID 11 setOptions numFrontEndVertToSkip 1
crackMgr localProb 10 crackMgrID 11 setOptions useMLSForSkippedFrontEndVert true

# Apply Moving Least Square (MLS) approximation to extracted SIFs:
# This must be used if we use option $numFrontEndVertToSkip.
crackMgr localProb 10 crackMgrID 11 setOptions useMLSForExtraction true

# The following is performed only if we need the code to extract 1 SIF only!
# Create crackGrowthLaw object
crackMgr localProb 10 crackMgrID 11 crGrowthPhysicsLaw create
crackMgr localProb 10 crackMgrID 11 crGrowthPhysicsLaw crGrowthType FatigueCrackGrowth
crackMgr localProb 10 crackMgrID 11 crGrowthPhysicsLaw crFrontIncrementType fixedDeltaA delta_a_max $delta_a_max_law

# Scaling law
crackMgr localProb 10 crackMgrID 11 crGrowthPhysicsLaw crFrontScalingLaw ParisLaw C 1.0 m 1.0 R 1.0

# The following option is used for SIMPLIFIED simulations only:
# If false, SIFs will not be constant along crack front.
crackMgr localProb 10 crackMgrID 11 crGrowthPhysicsLaw crFrontScalingLaw setConstFrontDispl true

# ************************** Solve Local Problem (LP) 10 **************************

# Assemble
set initTime [ clock clicks -milliseconds ]
assemble localProb 10
set finalTime [ clock clicks -milliseconds ]
set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for assemble command: $totalTime_pa s \n"

# Solve
set initTime [ clock clicks -milliseconds ]
solve localProb 10
set finalTime [ clock clicks -milliseconds ]
set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for solve command: $totalTime_sol s \n"

# Generate post-process VTK of local mesh
if { $generateParaviewFiles == true } {

  graphMesh localProb 10 create "$resultsFolder/PostProcessMesh_3-D_edgeCrack_LOC.vtu" VTKStyle
  graphMesh localProb 10 postProcess resolution 0 sampleInside

}

# ************************** Solve Enriched Global Problem (EGP) **************************

puts stdout "\n\n*********************************************"
puts stdout "\nStarting Enriched Global Problem (EGP)"
puts stdout "\n*********************************************\n\n"

# Set global-local enrichments:
setCompNodSpBasis inBBox \
                  xyzMin $x_min $y_min $z_min xyzMax $x_max $y_max $z_max \
                  spBasis 10

enrichApprox iso approxOrder $p_glob

# Use lower integration order to speed-up enriched global assembly
parSet lowerPorderForIntegEG true

# Assemble
set initTime [ clock clicks -milliseconds ]
assemble
set finalTime [ clock clicks -milliseconds ]
set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for assemble command: $totalTime_pa s \n"

# Solve
set initTime [ clock clicks -milliseconds ]
solve
set finalTime [ clock clicks -milliseconds ]
set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for solve command: $totalTime_sol s \n"

# ************************** G-L ITERATION for the FIRST STEP **************************
# The G-L iteration at the begining is necessary to improve the overall accuracy
# TBD: Test how many iterations are really needed in the pavement problems.
# For this edge-crack problem, the results considering 1 or 2 iterations were very similar to each other.
set max_iter_gl 2
for { set iter_gl 1 } { $iter_gl <= $max_iter_gl } { incr iter_gl } {

  puts stdout "\n\n*********************************************"
  puts stdout "\nStarting Enriched Global Problem (EGP) G-L Iter. $iter_gl"
  puts stdout "\n*********************************************\n\n"

  # ************************** Solve Local Problem (LP) 10 **************************

  # Assemble
  set initTime [ clock clicks -milliseconds ]
  assemble localProb 10 
  set finalTime [ clock clicks -milliseconds ]
  set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for assemble command: $totalTime_pa s \n"

  # Solve
  set initTime [ clock clicks -milliseconds ]
  solve localProb 10 
  set finalTime [ clock clicks -milliseconds ]
  set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for solve command: $totalTime_sol s \n"

  # Generate post-process VTK of local mesh
  if { $generateParaviewFiles == true } {

    # This will overwrite what it was written before:
    graphMesh localProb 10 create "$resultsFolder/PostProcessMesh_3-D_edgeCrack_LOC.vtu" VTKStyle
    graphMesh localProb 10 postProcess resolution 0 sampleInside

  }

  # Assemble
  set initTime [ clock clicks -milliseconds ]
  assemble
  set finalTime [ clock clicks -milliseconds ]
  set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for assemble command: $totalTime_pa s \n"

  # Solve
  set initTime [ clock clicks -milliseconds ]
  solve
  set finalTime [ clock clicks -milliseconds ]
  set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for solve command: $totalTime_sol s \n"

}

# Generate post-process VTK of enriched global mesh
if { $generateParaviewFiles == true } {

  # This will overwrite what it was written before:
  graphMesh create "$resultsFolder/PostProcessMesh_3-D_edgeCrack_EG.vtu" VTKStyle
  graphMesh postProcess resolution 0 sampleInside

}

# New GFEM-3:
# Defining third-scale local problem:
    
# ************************** Create 3rd-scale local problems 100 **************************
set localProbId 100

puts stdout "\n\n*********************************************"
puts stdout "\nStarting 3rd-scale Local Problem $localProbId"
puts stdout "\n*********************************************\n\n"

# The following coordinates define a bounding box for
# creating the 3rd-scale LP 100.
# This is problem dependent.
# TBD: Change it when solving the pavement problem.
set x_min_auxlp  0.24
set y_min_auxlp -0.49
set x_max_auxlp  1.76
set y_max_auxlp  0.49

set z_min_auxlp -0.71
set z_max_auxlp  0.71

createLocalProblem probID $localProbId \
                   xyzMin $x_min_auxlp $y_min_auxlp $z_min_auxlp xyzMax $x_max_auxlp $y_max_auxlp $z_max_auxlp \
                   numLayers 1 \
                   userBC "spr_local"

set crackMgrId 111
crackMgr localProb $localProbId crackMgrID $crackMgrId create crackSurfaceFromCrackMgr 1

# Since we are not using the 3rd-scale local problem solution to enrich
# the global problem, this is allowed.
parSet localProb $localProbId ShiftedSpBasisEnrichments true

# Refine mesh globally
set nref 0
for { set iref 0 } { $iref < $nref } { incr iref } {
  
  puts stdout "Beginning Step = $iref"
  refine localProb $localProbId all

}

# ************************** Parameters for PROCESSING the crack 111, 222, ... **************************

# Use branch-heaviside pair in elements fully cut by crack surface and in the geometrical enrichment zone:
# This improves accuracy of solution without adding too many dofs.
parSet localProb $localProbId useBranchHeavisideFnPairWithGeomEnrich true

# Setting to use planar cuts away from crack front
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions usePlanarCutsSignedDist true

# Setting branch function type
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions useBranchFn true
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions branchFunctionType curvedFrontOD6

# Enrich cylinder around crack front:
# Enrich all nodes within radius R of the front
# TIP: Radius = 4 * maxEdgeLen
# crackMgr localProb $localProbId crackMgrID $crackMgrId branchFnGeometricEnrich inCylinder radius 0.0013
#
# OR
#
# Enrich nodes in layers around the crack front:
crackMgr localProb $localProbId crackMgrID $crackMgrId branchFnGeometricEnrich numElemLayers 4

# Refine 3-D GFEM elements cut by crack surface
set nrefSurface 0
for { set iref 0 } { $iref < $nrefSurface } { incr iref } {
  
  crackMgr localProb $localProbId crackMgrID $crackMgrId refineMesh crackSurface
  
}

# Refine 3-D GFEM mesh around crack front:
# These 3rd-scale local problems must be refined enough in order to get good SIFs.
# TIP: Choose close to 3% of characteristic crack size.
# For example, radius of circular crack or length of an edge crack.
set maxEdgeElemLenReqAuxLP [ expr 0.03 * $initalCrackLength_a_o ]
crackMgr localProb $localProbId crackMgrID $crackMgrId refineMesh crackFronts maxEdgeLen $maxEdgeElemLenReqAuxLP

# Post-process crack surface at this step:
# It is the initial crack.
if { $generateParaviewFiles == true } {
  
  crackMgr localProb $localProbId crackMgrID $crackMgrId postSurface "$resultsFolder/Crack_surf.vtk" VTKStyle

}

set initTime [ clock clicks -milliseconds ]

# The following is set to true.
# Then, node snapping can be used along with GFEM-gl.
# Use this with caution!
parSet localProb $localProbId solveInKamikazeMode true

# Snap to the crack surface nodes that are close to it:
# This improves conditioning of global matrix.
# This was NOT extensively tested yet for GFEM-gl:
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions snapNodesToSurfAndFront snappingTolSurf 0.025 snappingTolFront 0.025

# Process crack:
# Cut elements, select enrichments, etc.
crackMgr localProb $localProbId crackMgrID $crackMgrId process

# Get min edge length
set minEdgeLen [ parGet localProb $localProbId minEdgeLen ]
puts stdout "\n ** minEdgeLen (3rd-scale LP $localProbId) = $minEdgeLen"

set finalTime [ clock clicks -milliseconds ]
set totalTime_cp [ expr $totalTime_cp + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for crack process command: $totalTime_cp s \n"

# Set polynomial order of the approximation
enrichApprox localProb $localProbId iso approxOrder $p_loc

# Enforce the search of local descendants while calculating global-local BCs
# since the enriched global solution is computed based on the other local problem
parSet localProb $localProbId enforceDescSearchGlobLocBC true

# Do NOT extract at crack front points on the domain boundary:
# This is usually required for surface breaking cracks since extraction 
# at boundary points is not always possible.
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions numFrontEndVertToSkip 1
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions useMLSForSkippedFrontEndVert true

# Apply Moving Least Square (MLS) approximation to extracted SIFs:
# This must be used if we use option $numFrontEndVertToSkip.
crackMgr localProb $localProbId crackMgrID $crackMgrId setOptions useMLSForExtraction true

# The following is performed only if we need the code to extract 1 SIF only!
# Create crackGrowthLaw object
crackMgr localProb $localProbId crackMgrID $crackMgrId crGrowthPhysicsLaw create
crackMgr localProb $localProbId crackMgrID $crackMgrId crGrowthPhysicsLaw crGrowthType FatigueCrackGrowth
crackMgr localProb $localProbId crackMgrID $crackMgrId crGrowthPhysicsLaw crFrontIncrementType fixedDeltaA delta_a_max $delta_a_max_law

# Scaling law
crackMgr localProb $localProbId crackMgrID $crackMgrId crGrowthPhysicsLaw crFrontScalingLaw ParisLaw C 1.0 m 1.0 R 1.0

# The following option is used for SIMPLIFIED simulations only:
# If false, SIFs will not be constant along crack front.
crackMgr localProb $localProbId crackMgrID $crackMgrId crGrowthPhysicsLaw crFrontScalingLaw setConstFrontDispl true

# ************************** Solve 3rd-scale local problems 100, 200, ... **************************

# Assemble
set initTime [ clock clicks -milliseconds ]
assemble localProb $localProbId
set finalTime [ clock clicks -milliseconds ]
set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for assemble command: $totalTime_pa s \n"

# Solve
set initTime [ clock clicks -milliseconds ]
solve localProb $localProbId
set finalTime [ clock clicks -milliseconds ]
set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for solve command: $totalTime_sol s \n"

# Generate post-process VTK of local mesh
if { $generateParaviewFiles == true } {

  graphMesh localProb $localProbId create "$resultsFolder/PostProcessMesh_3-D_edgeCrack_3rdScale_LOC$localProbId.vtu" VTKStyle
  graphMesh localProb $localProbId postProcess resolution 0 sampleInside

}

# ************************** Extract SIFs and compute viscoelastic G **************************

set initTime [ clock clicks -milliseconds ]

# Values for DCM extraction
set dist [ expr 4.0 * $minEdgeLen ]
set deltaD [ expr 0.1 * $minEdgeLen ]
set num 11

# Extract SIFs with DCM

puts stdout "\n Using DCM to extract SIF in 3rd-scale Local Problems:"
puts stdout " extrDist_min = $dist "
puts stdout " deltaDist = $deltaD "

if { $isVisco == true } {

  set s_time 0
  set e_time [ expr 10.0 * $t_peak ]

  # Recommended number of time steps
  set n_steps [ expr int( ( $e_time - $s_time ) / $t_peak * 10 ) ]
    
  # crackMgr crackMgrID 1 extract DCM \
  #                       extrDist_min $dist \
  #                       deltaDist $deltaD \
  #                       numExtrPoints $num \
  #                       start_time $s_time \
  #                       end_time $e_time \
  #                       nSteps $n_steps
  # 
  # crackMgr localProb 10 crackMgrID 11 extract DCM \
  #                                     extrDist_min $dist \
  #                                     deltaDist $deltaD \
  #                                     numExtrPoints $num \
  #                                     start_time $s_time \
  #                                     end_time $e_time \
  #                                     nSteps $n_steps

  # New GFEM-3:
  crackMgr localProb 100 crackMgrID 111 extract DCM \
                                        extrDist_min $dist \
                                        deltaDist $deltaD \
                                        numExtrPoints $num \
                                        start_time $s_time \
                                        end_time $e_time \
                                        nSteps $n_steps
} else {

  # crackMgr crackMgrID 1 extract DCM \
  #                       extrDist_min $dist \
  #                       deltaDist $deltaD \
  #                       numExtrPoints $num
  # 
  # crackMgr localProb 10 crackMgrID 11 extract DCM \
  #                                     extrDist_min $dist \
  #                                     deltaDist $deltaD \
  #                                     numExtrPoints $num

  # New GFEM-3:
  crackMgr localProb 100 crackMgrID 111 extract DCM \
                                        extrDist_min $dist \
                                        deltaDist $deltaD \
                                        numExtrPoints $num
}

set finalTime [ clock clicks -milliseconds ]
set totalTime_afp [ expr ( $finalTime - $initTime ) / 1.0e+03 ]
puts stdout "\nTime spent for extract command: $totalTime_afp s \n"
  
# ************************** Finished extracting SIFs with DCM **************************

# Compute total time estimate for simulation step
set totalTime [ expr $totalTime_cp + $totalTime_pa + $totalTime_sol + $totalTime_afp ]
puts stdout "\nTime spent for current step: $totalTime s\n"

exit