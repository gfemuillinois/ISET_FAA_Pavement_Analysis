# ************************** pavement_GFEMgl_propagation_constant_vel.tcl **************************

set crackFile "Fully_Crack"
set resultsFolder "res_propagation"
file mkdir $resultsFolder

set domainFile "D200_ImpDisp"

set isVisco true
set generateParaviewFiles true

# ************************** pavement_GFEMgl_propagation_constant_vel.tcl **************************

# Turn on creation of extract file
set env(DO_QA) ON

set num_threads 24
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

set x_min -0.05649
set y_min -0.02
set z_min  0.0
set x_max  0.05649
set y_max  0.13
set z_max  6.1

# Set quantities to be postprocessed
graphMesh appendName Solution
# graphMesh appendName VonMises

# If wants to see input geomesh and BCs
set printInputGMesh true

if { $printInputGMesh == true } {

  if { $generateParaviewFiles == true } {

    graphMesh create "$resultsFolder/pavement_geomesh.vtu" VTKStyle
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

# Maximum number of crack propagation steps
set num_steps 250

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

  graphMesh create "$resultsFolder/PostProcessMesh_pavement_IG.vtu" VTKStyle
  graphMesh postProcess resolution 0 sampleInside

}

# ************************** Setup for crack **************************
# The following parameters are related to the pavement's crack

# Limit for length of crack front edges (defined by edges of facets along crack front):
# Around 1/20 of length of initial crack front
set initialCrackFrontLength 6.1
set initalCrackLength_a_o 0.013
set frontEdgeLengthLimit [ expr $initialCrackFrontLength / 61.0 ]

# Frequency of surface remesh (if doing remesh):
# After first step and after remesh_surface_freq steps thereafter.
set remesh_surface_freq 50

# Maximum allowed number of remeshings to be performed during entire simulation
set max_num_remeshings $num_steps

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

# Parameters for scaling law along crack front (Paris type law):
# myC parameter must involve unit conversion here, IF needed.
set mym 1.6296
set myC 3.0e-01
set myR 0.0

# Fracture thoughness of material:
# Only relevant for stable crack propagation problems.
# This value of 1.0 will lead to crack arrest before running all simulations steps.
# set myGc 1.0

# Maximum DeltaA for each step
set delta_a_max_law [ expr 0.05 * $initalCrackLength_a_o ]

# Index of starting step:
# This is useful when re-starting a simulation using a .crf created by a previous run.
set starting_step 0

# ************************** Setup for crack for the Global Problem **************************

if { $starting_step == 0 } {

  crackMgr crackMgrID 1 create crackFile "$crackFile.crf"

} else {
  
  crackMgr crackMgrID 1 create crackFile "$crackFile.crf_$starting_step"

}

# Post-process first crack file
if { $starting_step == 0 } {

  crackMgr crackMgrID 1 postSurface "$resultsFolder/Crack_initial_surf.vtk"

} else {

  crackMgr crackMgrID 1 postSurface "$resultsFolder/Crack_initial_surf_$starting_step.vtk"

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

# Don't let very stretched surface facets appear along the crack front:
# Set this to activate use of PAS (Propagate And Stretch).
# Otherwise only PAE (Propagate And Extrude == add layer of facets along crack front) is used.
#
# If displacement at a front vertex is larger than this tolerance, a new layer
# of facets and vertices are created. Otherwise the front vertex is stretched.
#
# Setting high values next will enforce only PAS, which is fine for pavement simulations,
# and do not allow very strechted new facets to appear.
#
# The next parameters enforce PAS only:

crackMgr crackMgrID 1 setOptions frontStretchTol 110.0
crackMgr crackMgrID 1 setOptions accumulatedCrackFrontAdvanceLimit [ expr $num_steps * $delta_a_max_law * 2 ]
crackMgr crackMgrID 1 setOptions numStepsAccumDeltaA [ expr $num_steps + 1]

# ************************** Parameters for advancing the crack **************************

# Do NOT extract at crack front points on the boundary:
# This is usually required for surface breaking cracks since extraction
# at boundary points is not always possible.
crackMgr crackMgrID 1 setOptions numFrontEndVertToSkip 1
crackMgr crackMgrID 1 setOptions useMLSForSkippedFrontEndVert false

# Apply Moving Least Square (MLS) approximation to extracted SIFs:
# This must be used if we use option $numFrontEndVertToSkip.
crackMgr crackMgrID 1 setOptions useMLSForExtraction false

# Use MLS for the geometrical position of the front vertices:
# It improves the shape of the front and distributes the crack front vertices with a equal spacing.
crackMgr crackMgrID 1 setOptions useMLSForCrackFrontVertices false

# Set to ignore KII since we know there is no Mode II
# crackMgr crackMgrID 1 setOptions ignoreKII true

# Set to ignore KIII since we know there is no twisting (Mode III) in this example
# crackMgr crackMgrID 1 setOptions ignoreKIII true

# Create crackGrowthLaw object
crackMgr crackMgrID 1 crGrowthPhysicsLaw create
crackMgr crackMgrID 1 crGrowthPhysicsLaw crGrowthType FatigueCrackGrowth
crackMgr crackMgrID 1 crGrowthPhysicsLaw crFrontIncrementType fixedDeltaA delta_a_max $delta_a_max_law

# Scaling law
crackMgr crackMgrID 1 crGrowthPhysicsLaw crFrontScalingLaw ParisLaw C $myC m $mym R $myR

# The following option is usef for SIMPLIFIED simulations only:
crackMgr crackMgrID 1 crGrowthPhysicsLaw crFrontScalingLaw setConstFrontDispl true

# Kinking angle less than this will be ignored:
# Small angles are often the effect of noise KII.
# Default = 2 degrees.
# Not needed if we set ignoreKII == true.
# Large values will enforce propagation on a plane even if KII is not small.
crackMgr crackMgrID 1 crGrowthPhysicsLaw frontKinkingAngleLimitDEG 120.0

# Twisting angle less than this will be ignored:
# Small angles are often the effect of noise KIII.
# Not needed if we set ignoreKIII == true.
# Large values will enforce propagation on a plane even if KIII is not small.
crackMgr crackMgrID 1 crGrowthPhysicsLaw frontTwistingAngleLimitDEG 120.0

# Print this so we can check it
crackMgr crackMgrID 1 crGrowthPhysicsLaw print

# ************************** Options for VISCOELASTIC extraction **************************
# Set all visco options, method, load control, load type, etc.

if { $isVisco == true } {

  set loadFnType "haversine"
  # t_peak = 0.5 *(load period)
  # This next t_peak corresponds, for instance, to a period equals to 1800 s, or 30 min
  set t_peak 900.0

  # Tells the crack manager the name of the viscoelastic material defined in the .grf file
  crackMgr crackMgrID 1 setOptions useViscoExtractionCP viscoMatName "isoVisco4CP"

  # Method used to compute inverse Laplace transform:
  # If viscoFourierTransf is false: use Zakian's method (NOT recommended).
  crackMgr crackMgrID 1 setOptions viscoFourierTransf true

  # Used at ViscoCorrespondPrinc::get_lambda_s() method
  # viscoLoadType { constant | ramp | creep | hat | freq | haversine }
  crackMgr crackMgrID 1 setOptions viscoLoadType $loadFnType

  # Set the t_peak used for frequency loads
  crackMgr crackMgrID 1 setOptions viscoPeakTime $t_peak

  # Whether it is force or displacement controlled:
  # This is needed to compute ERR(t).
  crackMgr crackMgrID 1 setOptions viscoForceControlled false

  # CAD: This is the default. Setting just for clarity.
  # NOTE:
  # LaplaceTransf is the recommended option. Use others at your own risk...
  # Used at ViscoCorrespondPrinc::get_TimeExtract() method
  # viscoIntegrType { LaplaceTransf | GaussIntegr | incremental }
  crackMgr crackMgrID 1 setOptions viscoIntegrType LaplaceTransf

}

# ************************** Finish setting up cracks at global scale **************************
 
# Counter of number of remeshes performed
set remesh_counter 0

# Used to name the crack surface file of the last propagation step or arrested crack
set final_crack_prop_step $starting_step

set totalTime_cp  0.0
set totalTime_pa  0.0
set totalTime_sol 0.0

# Loop over the crack front advancement steps
for { set i_step $starting_step } { $i_step < $num_steps } { incr i_step } {

  puts stdout "\n\n*********************************************"
  puts stdout "\nStarting crack propagation step $i_step"
  puts stdout "\n*********************************************\n\n"

  if { [ expr ( $i_step - $starting_step ) % 2 ] == 0 } {

    puts stdout "\n\n*********************************************"
    puts stdout "\nStarting Local Problem (LP) 10 step $i_step"
    puts stdout "\n*********************************************\n\n"

    if { $i_step == $starting_step } {

        # Create a local problem from the global mesh
        createLocalProblem probID 10 \
                           xyzMin $x_min $y_min $z_min xyzMax $x_max $y_max $z_max \
                           numLayers 1 \
                           userBC "spr_local"

        crackMgr localProb 10 crackMgrID 11 create crackSurfaceFromCrackMgr 1

    }
    # This is not yet allowed:
    # parSet localProb 10 ShiftedSpBasisEnrichments true

    # Refine mesh globally
    set nref 1
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

    # Refine 3-D GFEM elements cut by crack surface
    set nrefSurface 1
    for { set iref 0 } { $iref < $nrefSurface } { incr iref } {

	crackMgr localProb 10 crackMgrID 11 refineMesh crackSurface
	
    }
    
    # REFINING AROUND BBOX
    # Important to reduce number of dofs.
    set crackLength [ expr $initalCrackLength_a_o + $delta_a_max_law * $i_step ]
    
    set xboxsize 0.0132
    set yboxsize 0.0132
    set zboxsize 0.1

    set z_min_refbox [ expr 3.05 - 0.5 * $zboxsize ]
    set z_max_refbox [ expr 3.05 + 0.5 * $zboxsize ]

    set nrefbbox 10
    for { set iref 0 } { $iref < $nrefbbox } { incr iref } {
      
      # BBox x dimensions:
      set x_min_refbox [ expr -0.5 * $xboxsize ]
      set x_max_refbox [ expr  0.5 * $xboxsize ]
      
      # BBox y dimensions:
      set y_min_refbox [ expr $crackLength - 0.5 * $yboxsize ]
      set y_max_refbox [ expr $crackLength + 0.5 * $yboxsize ]

      # Refine!
      refine localProb 10 hasNodeInBBox xyzMin $x_min_refbox $y_min_refbox $z_min_refbox \
                                        xyzMax $x_max_refbox $y_max_refbox $z_max_refbox

      # Resizing the bounding box to just refine elements that are close to the crack front
      set xboxsize [ expr $xboxsize / 1.3 ]
      set yboxsize [ expr $yboxsize / 1.3 ]
    }

    # Enrich cylinder around crack front:
    # Enrich all nodes within radius R of the front
    # TIP: Radius = 4 * maxEdgeLen
    crackMgr localProb 10 crackMgrID 11 branchFnGeometricEnrich inCylinder radius 0.0013
    #
    # OR
    #
    # Enrich nodes in layers around the crack front:
    # crackMgr localProb 10 crackMgrID 11 branchFnGeometricEnrich numElemLayers 4

    # Refine 3-D GFEM mesh around crack front:
    # TIP: Choose close to 3% of characteristic crack size.
    # For example, radius of circular crack or length of an edge crack.
    # set maxEdgeElemLenReq [ expr 0.08 * $initalCrackLength_a_o ]
    # crackMgr localProb 10 crackMgrID 11 refineMesh crackFronts maxEdgeLen $maxEdgeElemLenReq

    # Post-process crack surface at this step:
    # It is the initial crack or the one that was propagated at the previous step.
    crackMgr localProb 10 crackMgrID 11 postSurface "$resultsFolder/Crack_surf_$i_step.vtk" VTKStyle

    set initTime [ clock clicks -milliseconds ]

    # Snap to the crack surface nodes that are close to it:
    # This improves conditioning of global matrix.
    # This was NOT extensively tested yet for GFEM-gl,
    # but is important to be used for the pavement simulations!
    crackMgr localProb 10 crackMgrID 11 setOptions snapNodesToSurfAndFront snappingTolSurf 0.025 snappingTolFront 0.025
    
    # Process crack:
    # Cut elements, select enrichments, etc.
    crackMgr localProb 10 crackMgrID 11 process 

    # Get min edge length
    set minEdgeLen [ parGet localProb 10 minEdgeLen ]
    puts stdout "\n minEdgeLen = $minEdgeLen"

    set finalTime [ clock clicks -milliseconds ]
    set totalTime_cp [ expr $totalTime_cp + ( $finalTime - $initTime ) / 1.0e+03 ]
    puts stdout "\nTime spent for crack process command (step $i_step): $totalTime_cp s \n"

    # Set polynomial order of the approximation
    enrichApprox localProb 10 iso approxOrder $p_loc

    # Enforce the search of local descendants while calculating global-local BCs
    # since the enriched global solution is computed based on the other local problem
    parSet localProb 10 enforceDescSearchGlobLocBC true
    
    # ************************** Solve Local Problem (LP) 10 **************************
    
    # Assemble
    set initTime [ clock clicks -milliseconds ]
    assemble localProb 10 
    set finalTime [ clock clicks -milliseconds ]
    set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
    puts stdout "\nTime spent for assemble command (step $i_step): $totalTime_pa s \n"

    # Solve
    set initTime [ clock clicks -milliseconds ]
    solve localProb 10
    set finalTime [ clock clicks -milliseconds ]
    set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
    puts stdout "\nTime spent for solve command (step $i_step): $totalTime_sol s \n"

    # Generate post-process VTK of local mesh
    if { $generateParaviewFiles == true } {
      
      graphMesh localProb 10 create "$resultsFolder/PostProcessMesh_pavement_LOC_$i_step.vtu" VTKStyle
      graphMesh localProb 10 postProcess resolution 0 sampleInside

    }

    if { $i_step > $starting_step } {
        
        # Prepare crack manager for next step
        crackMgr localProb 100 crackMgrID 111 cleanUp afterPropagation

        # Unrefine local mesh
        unRefine localProb 100 toInitialMesh

    }

  } else {

    puts stdout "\n\n*********************************************"
    puts stdout "\nStarting Local Problem (LP) 100 Step $i_step"
    puts stdout "\n*********************************************\n\n"

    if { $i_step == [ expr $starting_step + 1 ] } {
        
        # Create a local problem from the global mesh
        createLocalProblem probID 100 \
                           xyzMin $x_min $y_min $z_min xyzMax $x_max $y_max $z_max \
                           numLayers 1 \
                           userBC "spr_local"

        crackMgr localProb 100 crackMgrID 111 create crackSurfaceFromCrackMgr 1
    }
    # This is not yet allowed:
    # parSet localProb 100 ShiftedSpBasisEnrichments true

    # Refine mesh globally
    set nref 1
    for { set iref 0 } { $iref < $nref } { incr iref } {
      
      puts stdout "Beginning Step = $iref"
      refine localProb 100 all

    }

    # ************************** Parameters for PROCESSING the crack 111 **************************

    # Use branch-heaviside pair in elements fully cut by crack surface and in the geometrical enrichment zone:
    # This improves accuracy of solution without adding too many dofs.
    parSet localProb 100 useBranchHeavisideFnPairWithGeomEnrich true

    # Setting to use planar cuts away from crack front
    crackMgr localProb 100 crackMgrID 111 setOptions usePlanarCutsSignedDist true

    # Setting branch function type
    crackMgr localProb 100 crackMgrID 111 setOptions useBranchFn true
    crackMgr localProb 100 crackMgrID 111 setOptions branchFunctionType curvedFrontOD6

    # Refine 3-D GFEM elements cut by crack surface
    set nrefSurface 1
    for { set iref 0 } { $iref < $nrefSurface } { incr iref } {

	crackMgr localProb 100 crackMgrID 111 refineMesh crackSurface
	
    }
    
    # REFINING AROUND BBOX
    # Important to reduce number of dofs.
    set crackLength [ expr $initalCrackLength_a_o + $delta_a_max_law * $i_step ]
    
    set xboxsize 0.0132
    set yboxsize 0.0132
    set zboxsize 0.1

    set z_min_refbox [ expr 3.05 - 0.5 * $zboxsize ]
    set z_max_refbox [ expr 3.05 + 0.5 * $zboxsize ]
    
    set nrefbbox 10
    for { set iref 0 } { $iref < $nrefbbox } { incr iref } {
      
      # BBox x dimensions:
      set x_min_refbox [ expr -0.5 * $xboxsize ]
      set x_max_refbox [ expr  0.5 * $xboxsize ]
      
      # BBox y dimensions:
      set y_min_refbox [ expr $crackLength - 0.5 * $yboxsize ]
      set y_max_refbox [ expr $crackLength + 0.5 * $yboxsize ]

      # Refine!
      refine localProb 100 hasNodeInBBox xyzMin $x_min_refbox $y_min_refbox $z_min_refbox \
                                         xyzMax $x_max_refbox $y_max_refbox $z_max_refbox

      # Resizing the bounding box to just refine elements that are close to the crack front
      set xboxsize [ expr $xboxsize / 1.3 ]
      set yboxsize [ expr $yboxsize / 1.3 ]
    }

    # Enrich cylinder around crack front:
    # Enrich all nodes within radius R of the front
    # TIP: Radius = 4 * maxEdgeLen
    crackMgr localProb 100 crackMgrID 111 branchFnGeometricEnrich inCylinder radius 0.0013
    #
    # OR
    #
    # Enrich nodes in layers around the crack front:
    # crackMgr localProb 100 crackMgrID 111 branchFnGeometricEnrich numElemLayers 4

    # Refine 3-D GFEM mesh around crack front:
    # TIP: Choose close to 3% of characteristic crack size.
    # For example, radius of circular crack or length of an edge crack.
    # set maxEdgeElemLenReq [ expr 0.08 * $initalCrackLength_a_o ]
    # crackMgr localProb 100 crackMgrID 111 refineMesh crackFronts maxEdgeLen $maxEdgeElemLenReq

    # Post-process crack surface at this step:
    # It is the initial crack or the one that was propagated at the previous step.
    crackMgr localProb 100 crackMgrID 111 postSurface "$resultsFolder/Crack_surf_$i_step.vtk" VTKStyle

    set initTime [ clock clicks -milliseconds ]

    # Snap to the crack surface nodes that are close to it:
    # This improves conditioning of global matrix.
    # This was NOT extensively tested yet for GFEM-gl,
    # but is important to be used for the pavement simulations!
    crackMgr localProb 100 crackMgrID 111 setOptions snapNodesToSurfAndFront snappingTolSurf 0.025 snappingTolFront 0.025

    # Process crack:
    # Cut elements, select enrichments, etc.
    crackMgr localProb 100 crackMgrID 111 process

    # Get min edge length
    set minEdgeLen [ parGet localProb 100 minEdgeLen ]
    puts stdout "\n minEdgeLen = $minEdgeLen"

    set finalTime [ clock clicks -milliseconds ]
    set totalTime_cp [ expr $totalTime_cp + ( $finalTime - $initTime ) / 1.0e+03 ]
    puts stdout "\nTime spent for crack process command (step $i_step): $totalTime_cp s \n"

    # Set polynomial order of the approximation
    enrichApprox localProb 100 iso approxOrder $p_loc

    # Enforce the search of local descendants while calculating global-local BCs
    # since the enriched global solution is computed based on the other local problem
    parSet localProb 100 enforceDescSearchGlobLocBC true

    # ************************** Solve Local Problem (LP) 100 **************************

    # Assemble
    set initTime [ clock clicks -milliseconds ]
    assemble localProb 100
    set finalTime [ clock clicks -milliseconds ]
    set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
    puts stdout "\nTime spent for assemble command (step $i_step): $totalTime_pa s \n"

    # Solve
    set initTime [ clock clicks -milliseconds ]
    solve localProb 100
    set finalTime [ clock clicks -milliseconds ]
    set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
    puts stdout "\nTime spent for solve command (step $i_step): $totalTime_sol s \n"

    # Generate post-process VTK of local mesh
    if { $generateParaviewFiles == true } {

      graphMesh localProb 100 create "$resultsFolder/PostProcessMesh_pavement_LOC_$i_step.vtu"  VTKStyle
      graphMesh localProb 100 postProcess resolution 0 sampleInside

    }

    # Prepare crack manager for next step
    crackMgr localProb 10 crackMgrID 11 cleanUp afterPropagation

    # Unrefine local mesh
    unRefine localProb 10 toInitialMesh

  }

  # ************************** Solve Enriched Global Problem (EGP) **************************

  puts stdout "\n\n*********************************************"
  puts stdout "\nStarting Enriched Global Problem (EGP) Step $i_step"
  puts stdout "\n*********************************************\n\n"

  # Set global-local enrichments:

  if { [ expr ( $i_step - $starting_step ) % 2 ] == 0 } {

    # delete spBasis 100
    setCompNodSpBasis inBBox \
                      xyzMin $x_min $y_min $z_min xyzMax $x_max $y_max $z_max \
                      spBasis 10

  } else {

    # delete spBasis 10
    setCompNodSpBasis inBBox \
                      xyzMin $x_min $y_min $z_min xyzMax $x_max $y_max $z_max \
                      spBasis 100

  }

  enrichApprox iso approxOrder $p_glob

  # Use lower integration order to speed-up enriched global assembly
  parSet lowerPorderForIntegEG true

  # Assemble
  set initTime [ clock clicks -milliseconds ]
  assemble
  set finalTime [ clock clicks -milliseconds ]
  set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for assemble command (step $i_step): $totalTime_pa s \n"

  # Solve
  set initTime [ clock clicks -milliseconds ]
  solve
  set finalTime [ clock clicks -milliseconds ]
  set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for solve command (step $i_step): $totalTime_sol s \n"

  # ************************** G-L ITERATION for the FIRST STEP **************************

  # The G-L iteration at the begining is necessary to improve the overall accuracy
  if { $i_step == $starting_step } {

    set max_iter_gl 2
    for { set iter_gl 1 } { $iter_gl <= $max_iter_gl } { incr iter_gl } {

        puts stdout "\n\n*********************************************"
        puts stdout "\nStarting Enriched Global Problem (EGP) Step $i_step G-L Iter $iter_gl"
        puts stdout "\n*********************************************\n\n"

        # ************************** Solve Local Problem (LP) 10 **************************

        # Assemble
        set initTime [ clock clicks -milliseconds ]
        assemble localProb 10 
        set finalTime [ clock clicks -milliseconds ]
        set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
        puts stdout "\nTime spent for assemble command (step $i_step): $totalTime_pa s \n"

        # Solve
        set initTime [ clock clicks -milliseconds ]
        solve localProb 10
        set finalTime [ clock clicks -milliseconds ]
        set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
        puts stdout "\nTime spent for solve command (step $i_step): $totalTime_sol s \n"

        # Generate post-process VTK of local mesh
        if { $generateParaviewFiles == true } {

          graphMesh localProb 10 create "$resultsFolder/PostProcessMesh_pavement_LOC_$i_step.vtu" VTKStyle
          graphMesh localProb 10 postProcess resolution 0 sampleInside

        }

        # Assemble
        set initTime [ clock clicks -milliseconds ]
        assemble
        set finalTime [ clock clicks -milliseconds ]
        set totalTime_pa [ expr $totalTime_pa + ( $finalTime - $initTime ) / 1.0e+03 ]
        puts stdout "\nTime spent for assemble command (step $i_step): $totalTime_pa s \n"

        # Solve
        set initTime [ clock clicks -milliseconds ]
        solve
        set finalTime [ clock clicks -milliseconds ]
        set totalTime_sol [ expr $totalTime_sol + ( $finalTime - $initTime ) / 1.0e+03 ]
        puts stdout "\nTime spent for solve command (step $i_step): $totalTime_sol s \n"

    }

  }

  # Generate post-process VTK of enriched global mesh
  if { $generateParaviewFiles == true } {

    graphMesh create "$resultsFolder/PostProcessMesh_pavement_EG_$i_step.vtu" VTKStyle
    graphMesh postProcess resolution 0 sampleInside

  }

  # ************************** ADVANCE FRONT WITH PHYSICS **************************

  if { $i_step != 0 } {
    
    # Do not store front unless this is the initial step
    crackMgr crackMgrID 1 setOptions storeCrackFront false

  }

  crackMgr crackMgrID 1 setOptions useRemeshCrackSurface false

  if { $useRemeshCrackSurface == true } {
    
    if { $i_step == 0 || [ expr $i_step % $remesh_surface_freq ] == 0 } {

      if { $remesh_counter < $max_num_remeshings } {
        
        puts stdout "\n\nCrack surface will be remeshed at this step = $i_step\n"
        crackMgr crackMgrID 1 setOptions useRemeshCrackSurface true
        incr remesh_counter

      }

    }

  }

  set initTime [ clock clicks -milliseconds ]

  # Values for DCM extraction
  set dist [ expr 4.0 * $minEdgeLen ] 
  set deltaD [ expr 0.1 * $minEdgeLen ]
  set num 11

  # Advance front with DCM
    
  puts stdout "\n Using DCM to extract SIF in enriched global problem"
  puts stdout " extrDist_min = $dist"
  puts stdout " deltaDist = $deltaD"

  if { $isVisco == true } {

    set s_time 0
    # This next e_time ( = twice the t_peak ) is required by the haversine load function!
    set e_time [ expr 2.0 * $t_peak ]

    # Recommended number of time steps
    set n_steps [ expr int( ( $e_time - $s_time ) / $t_peak * 10 ) ]
    
    crackMgr crackMgrID 1 advanceFront DCM \
                          extrDist_min $dist \
                          deltaDist $deltaD \
                          numExtrPoints $num \
                          start_time $s_time \
                          end_time $e_time \
                          nSteps $n_steps

  } else {

    crackMgr crackMgrID 1 advanceFront DCM \
                          extrDist_min $dist \
                          deltaDist $deltaD \
                          numExtrPoints $num

  }

  set finalTime [ clock clicks -milliseconds ]
  set totalTime_afp [ expr ( $finalTime - $initTime ) / 1.0e+03 ]
  puts stdout "\nTime spent for advanceFront command (step $i_step): $totalTime_afp s \n"

  # Crack surface was advanced using solution computed at this step:
  # This is the index for the new surface.
  set final_crack_prop_step [ expr $i_step + 1 ]
  
  # ************************** Finished ADVANCING crack **************************

  # Checking if the fracture has arrested
  set hasCrackArrested [ crackMgr crackMgrID 1 parGet isCrackArrested ]
  
  if { $hasCrackArrested == "true" } {
  
    puts stdout "\n********************************************************"
    puts stdout "The fracture has arrested. Finishing code and exiting...\n"
    break

  }

  # Compute total time estimate for simulation step
  set totalTime_step [ expr $totalTime_cp + $totalTime_pa + $totalTime_sol + $totalTime_afp ]
  puts stdout "\nTime spent for current step (step $i_step): $totalTime_step s\n"

  # Compute total time estimate for entire simulation until current step
  # (accounting for major steps only) and print to screen
  set totalTime [ expr $totalTime + $totalTime_step ]
  puts stdout "\nTime spent for entire simulation until current step (step $i_step): $totalTime s\n"

}

# Post-process the crack surface at the last step
crackMgr crackMgrID 1 postSurface "$resultsFolder/Crack_surf_$final_crack_prop_step.vtk" VTKStyle

# Print total time estimate for entire simulation step:
# Accounting for major steps only.
puts stdout "\nTime spent for entire simulation ($final_crack_prop_step steps): $totalTime s\n"

# END
exit
