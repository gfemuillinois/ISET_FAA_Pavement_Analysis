from re import S
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import *
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# ---------------------- Added: Load Configurations
conf_name = 'D200' # Configuration Name
load_case_x = 2 # Loading Scenario (Longitudinal Path)
load_case_y = 5 # Loading Scenario (Transverse Path)

# X & Y locations for the load conf
loadcase_x_dic = {'S75': [0.2832,1.2072, 2.1312, 3.0552, 3.9792, 4.9032, 5.8200],
                'D200': [0.2541, 1.1881, 2.1221, 3.0561, 3.9901, 4.9241, 5.8500]}
loadcase_y_dic = {'S75': [0, 0.719, 1.438, 2.157, 2.877],
                'D200': [0,0.7238, 1.4476, 2.1714, 2.8952]}
                
# -----------------------------------------
ModelName = 'Pavement_{}_Case{}{}'.format(conf_name, load_case_y, load_case_x)
JobName = '{}_Case{}{}'.format(conf_name, load_case_y, load_case_x)
GearsFileName = '{}.txt'.format(conf_name)

# -------------------- Geometry Properties
HMAThickness = 0.127
ConcreteThickness = 0.432
BaseThickness = 0.254
SlabWidth = 6.1
SlabLength = 6.1
Joint = 0.013

# -------------------- Material Properties
E_HMA = 1.0e010; nu_HMA = 0.3
E_Concrete = 2.76e010; nu_Concrete = 0.15
E_Base = 3.0e08; nu_Base = 0.35

# -------------------- Impose displacements in joints
isImposeDisp = False

# -------------------- Load Position
# Assume a two-dimensional coordinate system located at the geometric 
# center of the top face of the pavement, where the x-direction represents 
# the length of the pavement (2*SlabLength + Joint) and the y-direction 
# the depth (SlabWidth).
# Gear_x_center - Distance from the center of the gear in x-direction
#                 to the center of the pavement
# Gear_y_center - Distance from the center of the gear in y-direction
#                 to the center of the pavement

if not isImposeDisp: 
    #Gear_x_center = 0.5*SlabLength + 0.5*Joint
    #Gear_y_center = 0.0
    Gear_x_center = loadcase_x_dic[conf_name][load_case_x-1]
    Gear_y_center = loadcase_y_dic[conf_name][load_case_y-1]
    isElliptic = False
    # -------------------- Reading gears file
    f = open(GearsFileName,"r")
    f.seek(0,2)               
    endLocation = f.tell()  
    f.seek(0)                
    Tire_Box_x = 0.0
    Tire_Box_y = 0.0
    Wheels_in_gears_coordinates = []
    Pressure = 0.0
    line = "n"
    while(line != '' and f.tell() != endLocation):
        line = f.readline()
        if('BEGIN TIRE ' in line):
            beginLine = line
            line = f.readline()
            if('DIMENSIONS' in beginLine):
                while('END' not in line):
                    if(line!='\n'):
                        dimensions = line.split(";")
                        Tire_Box_x = float(dimensions[0])
                        Tire_Box_y = float(dimensions[1])
                    line = f.readline()
            elif('COORDINATES' in beginLine):
                while('END' not in line):
                    if(line!='\n'):
                        coords = line.split(";")
                        x_Coord = float(coords[0])
                        y_Coord = float(coords[-1])
                        Wheels_in_gears_coordinates.append([x_Coord, y_Coord])
                    line = f.readline()
            elif('PRESSURE' in beginLine):
                while('END' not in line):
                    if(line!='\n'):
                        Pressure = float(line)
                    line = f.readline()

 # -------------------- Mesh controls
LocalZoneRefSize = 0.1
TransitionZone = 0.25

maxElemSize = 0.3
thicknessRefSize = 0.06
transitionZoneMaxSize = 0.06
localZoneMaxSize = 0.025

# -------------------- additional variables
PavementThickness = HMAThickness+ConcreteThickness+BaseThickness
PavementLength = 2.0*SlabLength+Joint

# -------------------- Defining geometry
pavement = mdb.Model(modelType=STANDARD_EXPLICIT, name=ModelName)
pavement.ConstrainedSketch(name='__profile__', sheetSize=
    200.0)
pavement.sketches['__profile__'].rectangle(point1=(
    -SlabLength-0.5*Joint, -ConcreteThickness-BaseThickness), point2=(SlabLength+0.5*Joint, HMAThickness))
pavement.Part(dimensionality=THREE_D, name='Slabs', type=
    DEFORMABLE_BODY)
slabs = pavement.parts['Slabs']
slabs.BaseSolidExtrude(depth=SlabWidth, 
    sketch=pavement.sketches['__profile__'])
del pavement.sketches['__profile__']

# -------------------- Defining layers
temp = 0.5*PavementThickness+HMAThickness
pavement.ConstrainedSketch(gridSpacing=64.56, name=
    '__profile__', sheetSize=2582.71, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[4], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[7], 
    sketchOrientation=RIGHT, origin=(0.0, -temp, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-0.5*PavementLength, 
    temp), point2=(0.5*PavementLength, temp))
pavement.sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[6])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#10 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[7])
del pavement.sketches['__profile__']

temp = 0.5*(ConcreteThickness+BaseThickness)
pavement.ConstrainedSketch(gridSpacing=57.17, name=
    '__profile__', sheetSize=2287.01, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[0], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[3], 
    sketchOrientation=RIGHT, origin=(0.0, -temp, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-0.5*PavementLength, 
    temp-ConcreteThickness), point2=(0.5*PavementLength, temp-ConcreteThickness))
pavement.sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[9])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#1 ]', ), ), sketch=pavement.sketches['__profile__']
    , sketchUpEdge=slabs.edges[3])
del pavement.sketches['__profile__']
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#1 ]', ), ), edges=(
    slabs.edges[5], ), sweepPath=
    slabs.edges[12])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#1 ]', ), ), edges=(
    slabs.edges[13], ), sweepPath=
    slabs.edges[3])

# -------------------- Defining transition zone refinement
# Base Layer
temp = ConcreteThickness + 0.5*BaseThickness
sizeT = 0.5*TransitionZone + 0.5*Joint
pavement.ConstrainedSketch(gridSpacing=57.08, name=
    '__profile__', sheetSize=2283.28, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[8], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[10], 
    sketchOrientation=RIGHT, origin=(0.0, -temp, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-sizeT, 
    -0.5*temp), point2=(-sizeT, 0.5*temp))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[14])
pavement.sketches['__profile__'].Line(point1=(sizeT, -0.5*temp)
    , point2=(sizeT, 0.5*temp))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[15])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#100 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[10])
del pavement.sketches['__profile__']

# Concrete Layer
pavement.ConstrainedSketch(gridSpacing=57.12, name=
    '__profile__', sheetSize=2284.87, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[11], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[25], 
    sketchOrientation=RIGHT, origin=(0.0, -0.5*ConcreteThickness, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(sizeT, -0.5*ConcreteThickness)
    , point2=(sizeT, 0.5*ConcreteThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[22])
pavement.sketches['__profile__'].ParallelConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[2], entity2=
    pavement.sketches['__profile__'].geometry[22])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[12], 
    entity2=pavement.sketches['__profile__'].geometry[8])
pavement.sketches['__profile__'].Line(point1=(-sizeT, 
    -0.5*ConcreteThickness), point2=(-sizeT, 0.5*ConcreteThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[23])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[3], entity2=
    pavement.sketches['__profile__'].geometry[23])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[13], 
    entity2=pavement.sketches['__profile__'].geometry[8])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#800 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[25])
del pavement.sketches['__profile__']

# HMA Layer
pavement.ConstrainedSketch(gridSpacing=57.1, name=
    '__profile__', sheetSize=2284.21, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[18], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[37], 
    sketchOrientation=RIGHT, origin=(0.0, 0.5*HMAThickness, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(sizeT, -0.5*HMAThickness)
    , point2=(sizeT, 0.5*HMAThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[28])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[12], 
    entity2=pavement.sketches['__profile__'].geometry[28])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[14], 
    entity2=pavement.sketches['__profile__'].geometry[20])
pavement.sketches['__profile__'].Line(point1=(-sizeT, 
    -0.5*HMAThickness), point2=(-sizeT, 0.5*HMAThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[29])
pavement.sketches['__profile__'].ParallelConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[2], entity2=
    pavement.sketches['__profile__'].geometry[29])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[15], 
    entity2=pavement.sketches['__profile__'].geometry[20])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#40000 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[37])
del pavement.sketches['__profile__']

slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#1 ]', ), ), edges=(
    slabs.edges[14], 
    slabs.edges[17]), sweepPath=
    slabs.edges[22])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#8 ]', ), ), edges=(
    slabs.edges[29], 
    slabs.edges[31]), sweepPath=
    slabs.edges[3])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#40 ]', ), ), edges=(
    slabs.edges[35], 
    slabs.edges[39]), sweepPath=
    slabs.edges[7])

# -------------------- Defining local zone refinement
# HMA Layer
sizeT = 0.5*LocalZoneRefSize + 0.5*Joint
pavement.ConstrainedSketch(gridSpacing=0.01, name=
    '__profile__', sheetSize=0.58, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[40], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[6], 
    sketchOrientation=RIGHT, origin=(0.0, 0.5*HMAThickness, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-sizeT, 
    0.5*HMAThickness), point2=(-sizeT, -0.5*HMAThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[38])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[5], entity2=
    pavement.sketches['__profile__'].geometry[38])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[16], 
    entity2=pavement.sketches['__profile__'].geometry[5])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[17], 
    entity2=pavement.sketches['__profile__'].geometry[9])
pavement.sketches['__profile__'].Line(point1=(sizeT, 
    0.5*HMAThickness), point2=(sizeT, -0.5*HMAThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[39])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[5], entity2=
    pavement.sketches['__profile__'].geometry[39])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[18], 
    entity2=pavement.sketches['__profile__'].geometry[5])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[19], 
    entity2=pavement.sketches['__profile__'].geometry[9])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#0 #100 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[6])
del pavement.sketches['__profile__']

# Concrete layer
pavement.ConstrainedSketch(gridSpacing=0.03, name=
    '__profile__', sheetSize=1.23, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[37], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[29], 
    sketchOrientation=RIGHT, origin=(0.0, -0.5*ConcreteThickness, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-sizeT, 
    0.5*ConcreteThickness), point2=(-sizeT, -0.5*ConcreteThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[46])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[20], 
    entity2=pavement.sketches['__profile__'].geometry[31])
pavement.sketches['__profile__'].Line(point1=(sizeT, 
    0.5*ConcreteThickness), point2=(sizeT, -0.5*ConcreteThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[47])
pavement.sketches['__profile__'].ParallelConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[2], entity2=
    pavement.sketches['__profile__'].geometry[47])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[21], 
    entity2=pavement.sketches['__profile__'].geometry[31])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#0 #20 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[29])
del pavement.sketches['__profile__']

# Base layer
temp = 0.5*BaseThickness + ConcreteThickness
pavement.ConstrainedSketch(gridSpacing=0.03, name=
    '__profile__', sheetSize=1.46, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[29], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[44], 
    sketchOrientation=RIGHT, origin=(0.0, -temp, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-sizeT, 
    0.5*BaseThickness), point2=(-sizeT, -0.5*BaseThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[52])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[5], entity2=
    pavement.sketches['__profile__'].geometry[52])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[22], 
    entity2=pavement.sketches['__profile__'].geometry[27])
pavement.sketches['__profile__'].Line(point1=(sizeT, 0.5*BaseThickness)
    , point2=(sizeT, -0.5*BaseThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[53])
pavement.sketches['__profile__'].ParallelConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[2], entity2=
    pavement.sketches['__profile__'].geometry[53])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[23], 
    entity2=pavement.sketches['__profile__'].geometry[27])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#20000000 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[44])
del pavement.sketches['__profile__']

slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#2 ]', ), ), edges=(
    slabs.edges[14], 
    slabs.edges[17]), sweepPath=
    slabs.edges[25])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#10 ]', ), ), edges=(
    slabs.edges[30], 
    slabs.edges[32]), sweepPath=
    slabs.edges[49])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#200 ]', ), ), edges=(
    slabs.edges[36], 
    slabs.edges[40]), sweepPath=
    slabs.edges[67])

# -------------------- Defining joint
# HMA Layer
pavement.ConstrainedSketch(gridSpacing=0.03, name=
    '__profile__', sheetSize=1.39, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[66], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[34], 
    sketchOrientation=RIGHT, origin=(0.0, 0.5*HMAThickness, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-0.5*Joint, 
    0.5*HMAThickness), point2=(-0.5*Joint, -0.5*HMAThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[62])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[29], 
    entity2=pavement.sketches['__profile__'].geometry[62])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[24], 
    entity2=pavement.sketches['__profile__'].geometry[29])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[25], 
    entity2=pavement.sketches['__profile__'].geometry[33])
pavement.sketches['__profile__'].Line(point1=(0.5*Joint, 
    0.5*HMAThickness), point2=(0.5*Joint, -0.5*HMAThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[63])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[29], 
    entity2=pavement.sketches['__profile__'].geometry[63])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[26], 
    entity2=pavement.sketches['__profile__'].geometry[29])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[27], 
    entity2=pavement.sketches['__profile__'].geometry[33])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#0:2 #4 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[34])
del pavement.sketches['__profile__']

# Concrete layer
pavement.ConstrainedSketch(gridSpacing=0.02, name=
    '__profile__', sheetSize=1.01, transform=
    pavement.parts['Slabs'].MakeSketchTransform(
    sketchPlane=pavement.parts['Slabs'].faces[63], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=pavement.parts['Slabs'].edges[28], 
    sketchOrientation=RIGHT, origin=(0.0, -0.5*ConcreteThickness, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-0.5*Joint, 
    0.5*ConcreteThickness), point2=(-0.5*Joint, -0.5*ConcreteThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[70])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[5], entity2=
    pavement.sketches['__profile__'].geometry[70])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[28], 
    entity2=pavement.sketches['__profile__'].geometry[13])
pavement.sketches['__profile__'].Line(point1=(0.5*Joint, 
    0.5*ConcreteThickness), point2=(0.5*Joint, -0.5*ConcreteThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[71])
pavement.sketches['__profile__'].ParallelConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[2], entity2=
    pavement.sketches['__profile__'].geometry[71])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[29], 
    entity2=pavement.sketches['__profile__'].geometry[13])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#0 #80000000 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[28])
del pavement.sketches['__profile__']

# Base layer
temp = ConcreteThickness + 0.5*BaseThickness
pavement.ConstrainedSketch(gridSpacing=0.03, name=
    '__profile__', sheetSize=1.39, transform=
    slabs.MakeSketchTransform(
    sketchPlane=slabs.faces[55], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=slabs.edges[16], 
    sketchOrientation=RIGHT, origin=(0.0, -temp, SlabWidth)))
slabs.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, sketch=
    pavement.sketches['__profile__'])
pavement.sketches['__profile__'].Line(point1=(-0.5*Joint, 
    0.5*BaseThickness), point2=(-0.5*Joint, -0.5*BaseThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[76])
pavement.sketches['__profile__'].PerpendicularConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[5], entity2=
    pavement.sketches['__profile__'].geometry[76])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[30], 
    entity2=pavement.sketches['__profile__'].geometry[53])
pavement.sketches['__profile__'].Line(point1=(0.5*Joint, 
    0.5*BaseThickness), point2=(0.5*Joint, -0.5*BaseThickness))
pavement.sketches['__profile__'].VerticalConstraint(
    addUndoState=False, entity=
    pavement.sketches['__profile__'].geometry[77])
pavement.sketches['__profile__'].ParallelConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].geometry[2], entity2=
    pavement.sketches['__profile__'].geometry[77])
pavement.sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    pavement.sketches['__profile__'].vertices[31], 
    entity2=pavement.sketches['__profile__'].geometry[53])
slabs.PartitionFaceBySketch(faces=
    slabs.faces.getSequenceFromMask((
    '[#0 #800000 ]', ), ), sketch=
    pavement.sketches['__profile__'], sketchUpEdge=
    slabs.edges[16])
del pavement.sketches['__profile__']

slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#80 ]', ), ), edges=(
    slabs.edges[14], 
    slabs.edges[17]), sweepPath=
    slabs.edges[37])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#10 ]', ), ), edges=(
    slabs.edges[30], 
    slabs.edges[32]), sweepPath=
    slabs.edges[17])
slabs.PartitionCellBySweepEdge(cells=
    slabs.cells.getSequenceFromMask((
    '[#20 ]', ), ), edges=(
    slabs.edges[36], 
    slabs.edges[40]), sweepPath=
    slabs.edges[49])
slabs.RemoveFaces(deleteCells=False, 
    faceList=
    slabs.faces.getSequenceFromMask(
    mask=('[#0:2 #800000 ]', ), ))
slabs.RemoveFaces(deleteCells=False, 
    faceList=
    slabs.faces.getSequenceFromMask(
    mask=('[#0:2 #10000000 ]', ), ))

# -------------------- Defining materials
pavement.Material(name='HMA')
pavement.materials['HMA'].Elastic(table=((E_HMA, nu_HMA), ))
pavement.Material(name='Concrete')
pavement.materials['Concrete'].Elastic(table=((E_Concrete, nu_Concrete), 
    ))
pavement.Material(name='Base')
pavement.materials['Base'].Elastic(table=((E_Base, nu_Base), ))
pavement.HomogeneousSolidSection(material='HMA', name=
    'HMA', thickness=None)
pavement.HomogeneousSolidSection(material='Concrete', 
    name='Concrete', thickness=None)
pavement.HomogeneousSolidSection(material='Base', name=
    'Base', thickness=None)
slabs.Set(cells=
    slabs.cells.getSequenceFromMask((
    '[#81e18 ]', ), ), name='HMA')
slabs.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    slabs.sets['HMA'], sectionName='HMA'
    , thicknessAssignment=FROM_SECTION)
slabs.Set(cells=
    slabs.cells.getSequenceFromMask((
    '[#46184 ]', ), ), name='Concrete')
slabs.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    slabs.sets['Concrete'], sectionName=
    'Concrete', thicknessAssignment=FROM_SECTION)
slabs.Set(cells=
    slabs.cells.getSequenceFromMask((
    '[#38063 ]', ), ), name='Base')
slabs.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    slabs.sets['Base'], sectionName=
    'Base', thicknessAssignment=FROM_SECTION)

# -------------------- Sets and surfaces 
# Defining SURFACES for springs BCs
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#0:2 #800000 ]', ), ), name='Base_BC_0')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#0:2 #10000 ]', ), ), name='Concrete_BC_0')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#0:2 #200000 ]', ), ), name='HMA_BC_0')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#300000 #c000 #403000 ]', ), ), name='Base_BC_1')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#c00000 #30000 #c00 ]', ), ), name='Concrete_BC_1')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#3000000 #c0000 #8000300 ]', ), ), name='HMA_BC_1')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#0:2 #20000 ]', ), ), name='Base_BC_2')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#0:2 #100000 ]', ), ), name='Concrete_BC_2')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#0:2 #2000000 ]', ), ), name='HMA_BC_2')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#c0000030 #0 #800c ]', ), ), name='Base_BC_3')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#c00 #30000030 ]', ), ), name='Concrete_BC_3')
slabs.Surface(side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#c000 #3000300 #80000 ]', ), ), name='HMA_BC_3')
slabs.Surface(name='Bottom_face_BC', 
    side1Faces=
    slabs.faces.getSequenceFromMask(
    ('[#3000000c #0 #40000c0 ]', ), ))

# Defining Sets for mesh refinement
slabs.Set(edges=
    slabs.edges.findAt(
    ((-0.5*PavementLength, HMAThickness, 0.5*SlabWidth), ), ((-0.5*PavementLength, 0.0, 0.5*SlabWidth), ), 
    ((-0.5*PavementLength, -ConcreteThickness, 0.5*SlabWidth), ), ((0.5*PavementLength, -ConcreteThickness, 0.5*SlabWidth), ),
    ((-0.5*PavementLength, -(ConcreteThickness+BaseThickness), 0.5*SlabWidth), ), ((0.5*PavementLength, HMAThickness, 0.5*SlabWidth), ), 
    ((0.5*PavementLength, 0.0, 0.5*SlabWidth), ), ((0.5*PavementLength, -(ConcreteThickness+BaseThickness), 0.5*SlabWidth), ), ), 
    name='Ext_Z_Lines')
temp = 0.5*TransitionZone + 0.5*Joint + 0.5*(SlabLength - 0.5*TransitionZone)
slabs.Set(edges=
    slabs.edges.findAt(
    ((-temp, 0.0, 0.0), ), ((temp, HMAThickness, 0.0), ), ((-temp, -ConcreteThickness, 0.0), ), 
    ((temp, -ConcreteThickness, SlabWidth), ), ((temp, HMAThickness, SlabWidth), ), ((temp, 0.0, SlabWidth), ), 
    ((-temp, -(ConcreteThickness+BaseThickness), SlabWidth), ), ((-temp, -(ConcreteThickness+BaseThickness), 
    0.0), ), ), name='Ext_X_Lines_dir1')
slabs.Set(edges=
    slabs.edges.findAt(
    ((-temp, HMAThickness, SlabWidth), ), ((-temp, HMAThickness, 0.0), ), ((temp, 0.0, 0.0), ), 
    ((temp, -ConcreteThickness, 0.0), ), ((-temp, 0.0, SlabWidth), ), ((temp, -(ConcreteThickness+BaseThickness), 0.0), ), 
    ((-temp, -ConcreteThickness, SlabWidth), ), ((temp, -(ConcreteThickness+BaseThickness), SlabWidth), ), ), 
    name='Ext_X_Lines_dir2')
slabs.Set(edges=
    slabs.edges.findAt(
    ((-0.5*PavementLength, 0.5*HMAThickness, 0.0), ),((-0.5*PavementLength, -0.5*ConcreteThickness, 0.0), ), ((-0.5*PavementLength, -ConcreteThickness-0.5*BaseThickness, 0.0), ), 
    ((-0.5*(TransitionZone + Joint), 0.5*HMAThickness, 0.0), ),((-0.5*(TransitionZone + Joint), -0.5*ConcreteThickness, 0.0), ), ((-0.5*(TransitionZone + Joint), -ConcreteThickness-0.5*BaseThickness, 0.0), ), 
    ((0.5*(TransitionZone + Joint), 0.5*HMAThickness, 0.0), ),((0.5*(TransitionZone + Joint), -0.5*ConcreteThickness, 0.0), ), ((0.5*(TransitionZone + Joint), -ConcreteThickness-0.5*BaseThickness, 0.0), ), 
    ((0.5*PavementLength, 0.5*HMAThickness, 0.0), ),((0.5*PavementLength, -0.5*ConcreteThickness, 0.0), ), ((0.5*PavementLength, -ConcreteThickness-0.5*BaseThickness, 0.0), ), 
    ((-0.5*PavementLength, 0.5*HMAThickness, SlabWidth), ),((-0.5*PavementLength, -0.5*ConcreteThickness, SlabWidth), ), ((-0.5*PavementLength, -ConcreteThickness-0.5*BaseThickness, SlabWidth), ), 
    ((-0.5*(TransitionZone + Joint), 0.5*HMAThickness, SlabWidth), ),((-0.5*(TransitionZone + Joint), -0.5*ConcreteThickness, SlabWidth), ), ((-0.5*(TransitionZone + Joint), -ConcreteThickness-0.5*BaseThickness, SlabWidth), ), 
    ((0.5*(TransitionZone + Joint), 0.5*HMAThickness, SlabWidth), ),((0.5*(TransitionZone + Joint), -0.5*ConcreteThickness, SlabWidth), ), ((0.5*(TransitionZone + Joint), -ConcreteThickness-0.5*BaseThickness, SlabWidth), ), 
    ((0.5*PavementLength, 0.5*HMAThickness, SlabWidth), ),((0.5*PavementLength, -0.5*ConcreteThickness, SlabWidth), ), ((0.5*PavementLength, -ConcreteThickness-0.5*BaseThickness, SlabWidth), ),  ), 
    name='Ext_Y_Lines')
slabs.Set(edges=
    slabs.edges.findAt(
    ((-0.5*Joint, 0.5*HMAThickness, 0.0), ), ((0.5*Joint, 0.5*HMAThickness, 0.0), ), 
    ((-0.5*Joint, 0.5*HMAThickness, SlabWidth), ), ((0.5*Joint, 0.5*HMAThickness, SlabWidth), ),
    ((-0.5*Joint, HMAThickness, 0.5*SlabWidth), ), ((0.5*Joint, HMAThickness, 0.5*SlabWidth), ),
    ((-0.5*Joint, 0.0, 0.5*SlabWidth), ), ((0.5*Joint, 0.0, 0.5*SlabWidth), ),
    ((0.0, 0.0, 0.0), ), ((0.0, HMAThickness, 0.0), ), 
    ((0.0, 0.0, SlabWidth), ), ((0.0, HMAThickness, SlabWidth), ), ), 
    name='Joint_Lines')
temp = 0.5*(LocalZoneRefSize + Joint)
slabs.Set(edges=
    slabs.edges.findAt(
    ((-temp, 0.5*HMAThickness, 0.0), ), ((temp, 0.5*HMAThickness, 0.0), ), 
    ((-temp, 0.5*HMAThickness, SlabWidth), ), ((temp, 0.5*HMAThickness, SlabWidth), ),
    ((-temp, HMAThickness, 0.5*SlabWidth), ), ((temp, HMAThickness, 0.5*SlabWidth), ),
    ((-temp, 0.0, 0.5*SlabWidth), ), ((temp, 0.0, 0.5*SlabWidth), ), ), 
    name='Ext_Local_Ref_Lines')
temp = 0.5*(TransitionZone + Joint)
slabs.Set(edges=
    slabs.edges.findAt(
    ((-temp, HMAThickness, 0.5*SlabWidth), ), ((-temp, 0.0, 0.5*SlabWidth), ), 
    ((-temp, -ConcreteThickness, 0.5*SlabWidth), ), ((temp, -ConcreteThickness, 0.5*SlabWidth), ),
    ((-temp, -(ConcreteThickness+BaseThickness), 0.5*SlabWidth), ), ((temp, HMAThickness, 0.5*SlabWidth), ), 
    ((temp, 0.0, 0.5*SlabWidth), ), ((temp, -(ConcreteThickness+BaseThickness), 0.5*SlabWidth), ), ), 
    name='Int_Z_Lines')
temp = -ConcreteThickness - 0.5*BaseThickness
slabs.Set(edges=
    slabs.edges.findAt(
    ((-0.5*Joint, temp, 0.0), ), ((0.5*Joint, temp, 0.0), ), 
    ((-0.5*Joint, temp, SlabWidth), ), ((0.5*Joint, temp, SlabWidth), ),
    ((-0.5*Joint-0.5*LocalZoneRefSize, temp, 0.0), ), ((0.5*Joint+0.5*LocalZoneRefSize, temp, 0.0), ), 
    ((-0.5*Joint-0.5*LocalZoneRefSize, temp, SlabWidth), ), ((0.5*Joint+0.5*LocalZoneRefSize, temp, SlabWidth), ), 
    ((-0.5*Joint-0.5*LocalZoneRefSize-0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness, 0.0), ), ((-0.5*Joint-0.5*LocalZoneRefSize-0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness-BaseThickness, 0.0), ), 
    ((-0.5*Joint-0.5*LocalZoneRefSize-0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness, SlabWidth), ), ((-0.5*Joint-0.5*LocalZoneRefSize-0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness-BaseThickness, SlabWidth), ),
    ((-0.5*Joint-0.25*LocalZoneRefSize, -ConcreteThickness, 0.0), ), ((-0.5*Joint-0.25*LocalZoneRefSize, -ConcreteThickness-BaseThickness, 0.0), ), 
    ((-0.5*Joint-0.25*LocalZoneRefSize, -ConcreteThickness, SlabWidth), ), ((-0.5*Joint-0.25*LocalZoneRefSize, -ConcreteThickness-BaseThickness, SlabWidth), ),  
    ((0.5*Joint+0.5*LocalZoneRefSize+0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness, 0.0), ), ((0.5*Joint+0.5*LocalZoneRefSize+0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness-BaseThickness, 0.0), ), 
    ((0.5*Joint+0.5*LocalZoneRefSize+0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness, SlabWidth), ), ((0.5*Joint+0.5*LocalZoneRefSize+0.25*(TransitionZone-LocalZoneRefSize), -ConcreteThickness-BaseThickness, SlabWidth), ),
    ((0.5*Joint+0.25*LocalZoneRefSize, -ConcreteThickness, 0.0), ), ((0.5*Joint+0.25*LocalZoneRefSize, -ConcreteThickness-BaseThickness, 0.0), ), 
    ((0.5*Joint+0.25*LocalZoneRefSize, -ConcreteThickness, SlabWidth), ), ((0.5*Joint+0.25*LocalZoneRefSize, -ConcreteThickness-BaseThickness, SlabWidth), ),
    ((-0.5*Joint, -ConcreteThickness-BaseThickness, 0.5*SlabWidth), ), ((0.5*Joint, -ConcreteThickness-BaseThickness, 0.5*SlabWidth), ),
    ((-0.5*LocalZoneRefSize-0.5*Joint, -ConcreteThickness-BaseThickness, 0.5*SlabWidth), ), ((0.5*LocalZoneRefSize+0.5*Joint, -ConcreteThickness-BaseThickness, 0.5*SlabWidth), ),
    ((-0.5*Joint, -ConcreteThickness, 0.5*SlabWidth), ), ((0.5*Joint, -ConcreteThickness, 0.5*SlabWidth), ),
    ((-0.5*LocalZoneRefSize-0.5*Joint, -ConcreteThickness, 0.5*SlabWidth), ), ((0.5*LocalZoneRefSize+0.5*Joint, -ConcreteThickness, 0.5*SlabWidth), ), ), 
    name='Base_RefZone_Lines')
temp=0.5*Joint + 0.25*LocalZoneRefSize
slabs.Set(edges=
    slabs.edges.findAt(
    ((-temp, 0.0, 0.0), ), ((temp, HMAThickness, 0.0), ), ((temp, 0.0, SlabWidth), ), ((
    temp, HMAThickness, SlabWidth), ), ), 
    name='LocalRefZone_XLines_dir1')
slabs.Set(edges=
    slabs.edges.findAt(
    ((temp, 0.0, 0.0), ), ((-temp, HMAThickness, 0.0), ), ((-temp, 0.0, SlabWidth), ), ((
    -temp, HMAThickness, SlabWidth), ), ), 
    name='LocalRefZone_XLines_dir2')
temp=0.5*Joint + 0.25*LocalZoneRefSize + 0.25*(TransitionZone-LocalZoneRefSize)
slabs.Set(edges=
    slabs.edges.findAt(
    ((-temp, 0.0, 0.0), ), ((temp, HMAThickness, 0.0), ), ((temp, 0.0, SlabWidth), 
    ), ((temp, HMAThickness, SlabWidth), ), ), 
    name='TransZone_XLines_dir1')
slabs.Set(edges=
    slabs.edges.findAt(
    ((temp, 0.0, 0.0), ), ((-temp, HMAThickness, 0.0), ), ((-temp, 0.0, SlabWidth), 
    ), ((-temp, HMAThickness, SlabWidth), ), ), 
    name='TransZone_XLines_dir2')
slabs.Set(edges=
    slabs.edges.findAt(
    ((0.5*Joint, -0.5*ConcreteThickness, SlabWidth), ), ((-0.5*Joint, -0.5*ConcreteThickness, SlabWidth), ), ),
    name='Joint_Y_dir1')
slabs.Set(edges=
    slabs.edges.findAt(
    ((0.5*Joint, -0.5*ConcreteThickness, 0.0), ), ((-0.5*Joint, -0.5*ConcreteThickness, 0.0), ), ),  
    name='Joint_Y_dir2')
slabs.Set(edges=
    slabs.edges.findAt(
    ((0.5*Joint+0.5*LocalZoneRefSize, -0.5*ConcreteThickness, SlabWidth), ), ((-0.5*Joint-0.5*LocalZoneRefSize, -0.5*ConcreteThickness, SlabWidth), ), ), 
    name='LocalZone_Y_dir1')
slabs.Set(edges=
    slabs.edges.findAt(
    ((-0.5*Joint-0.5*LocalZoneRefSize, -0.5*ConcreteThickness, 0.0), ), ((0.5*Joint+0.5*LocalZoneRefSize, -0.5*ConcreteThickness, 0.0), ), ),  
    name='LocalZone_Y_dir2')
slabs.Surface(side1Faces=
    slabs.faces.findAt(((-0.5*Joint, 
    -0.5*ConcreteThickness, 0.5*SlabWidth), )), name='Joint_Face_1')
slabs.Surface(side1Faces=
    slabs.faces.findAt(((0.5*Joint, 
    -0.5*ConcreteThickness, 0.5*SlabWidth), )), name='Joint_Face_2')

# Defining SURFACES for pressure BCs
# Drawing pressure surfaces
if not isImposeDisp:
    pavement.ConstrainedSketch(gridSpacing=0.3, name=
        '__profile__', sheetSize=12.2, transform=
        slabs.MakeSketchTransform(
        sketchPlane=slabs.faces[18], 
        sketchPlaneSide=SIDE1, 
        sketchUpEdge=slabs.edges[94], 
        sketchOrientation=RIGHT, origin=(0.0, HMAThickness, 0.5*SlabWidth)))
    slabs.projectReferencesOntoSketch(
        filter=COPLANAR_EDGES, sketch=
        pavement.sketches['__profile__'])
    for i in range(len(Wheels_in_gears_coordinates)):
        x_wheel_center = Gear_x_center + Wheels_in_gears_coordinates[i][0]
        y_wheel_center = Gear_y_center + Wheels_in_gears_coordinates[i][1]
        x_min = x_wheel_center - 0.5*Tire_Box_x
        x_max = x_wheel_center + 0.5*Tire_Box_x
        y_min = y_wheel_center - 0.5*Tire_Box_y
        y_max = y_wheel_center + 0.5*Tire_Box_y
        if(isElliptic):
            pavement.sketches['__profile__'].EllipseByCenterPerimeter(
            axisPoint1=(x_max, y_wheel_center), axisPoint2=(x_wheel_center, y_max), center=(x_wheel_center, 
            y_wheel_center))
        else:
            pavement.sketches['__profile__'].rectangle(point1=(x_min, 
                y_min), point2=(x_max, y_max))
    slabs.PartitionFaceBySketch(faces=
        slabs.faces.getSequenceFromMask((
        '[#c0000 #c03000 #1000000 ]', ), ), sketch=
        pavement.sketches['__profile__'], sketchUpEdge=
        slabs.edges[94])
    del pavement.sketches['__profile__']
    # Defining pressure surface
    wheels = []
    wheels_lines = []
    wheels_lines_Joint = []
    wheels_lines_LocalProb = []
    wheels_lines_TransZone = []
    HalfTransZoneSize = 0.5*TransitionZone+0.5*Joint
    HalfLocalZoneSize = 0.5*LocalZoneRefSize+0.5*Joint
    Domains = [[-0.5*PavementLength, -HalfTransZoneSize],
            [-HalfTransZoneSize, -HalfLocalZoneSize],
            [-HalfLocalZoneSize, -0.5*Joint],
            [-0.5*Joint, 0.5*Joint],
            [ 0.5*Joint, HalfLocalZoneSize],
            [ HalfLocalZoneSize, HalfTransZoneSize],
            [ HalfTransZoneSize, 0.5*PavementLength]]
    for i in range(len(Wheels_in_gears_coordinates)):
        x_wheel_center = Gear_x_center + Wheels_in_gears_coordinates[i][0]
        y_wheel_center = - Gear_y_center - Wheels_in_gears_coordinates[i][1] + SlabWidth*0.5
        x_min = x_wheel_center - 0.5*Tire_Box_x
        x_max = x_wheel_center + 0.5*Tire_Box_x
        y_min = y_wheel_center - 0.5*Tire_Box_y
        y_max = y_wheel_center + 0.5*Tire_Box_y
        x_min_l = x_min
        x_max_l = x_max
        y_min_l = y_min
        y_max_l = y_max
        if(abs(x_min_l)>PavementLength*0.5 and abs(x_max_l)<PavementLength*0.5):
            x_min_l = -PavementLength*0.5
        if(abs(x_max_l)>PavementLength*0.5 and abs(x_min_l)<PavementLength*0.5):
            x_max_l = PavementLength*0.5
        if(y_min_l<0.0 and y_max_l>0.0):
            y_min_l = 0.0
        if(y_max_l>SlabWidth and y_min_l<SlabWidth):
            y_max_l = SlabWidth
        if(y_min<0.0):
            y_wheel_center = 0.5*y_max
        if(y_max>SlabWidth):
            y_wheel_center = 0.5*(SlabWidth-y_min)
        if(abs(x_min_l) <= 0.5*Joint):
            wheels_lines_Joint.append(((x_min_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        elif(abs(x_min_l) <= 0.5*Joint + 0.5*LocalZoneRefSize):
            wheels_lines_LocalProb.append(((x_min_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        elif(abs(x_min_l) <= 0.5*Joint + 0.5*TransitionZone):
            wheels_lines_TransZone.append(((x_min_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        else:
            wheels_lines.append(((x_min_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        if(abs(x_max_l) <= 0.5*Joint):
            wheels_lines_Joint.append(((x_max_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        elif(abs(x_max_l) <= 0.5*Joint + 0.5*LocalZoneRefSize):
            wheels_lines_LocalProb.append(((x_max_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        elif(abs(x_max_l) <= 0.5*Joint + 0.5*TransitionZone):
            wheels_lines_TransZone.append(((x_max_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        else:
            wheels_lines.append(((x_max_l, HMAThickness, (y_max - y_min)*0.5+y_min),))
        minDist = min(abs(x_min), abs(x_max), abs(x_wheel_center))
        for domain in Domains:
            minDomain = domain[0]
            maxDomain = domain[1]
            if(x_min <= maxDomain):
                if(x_max <= maxDomain):
                    x_wheel_center = 0.5*(x_min - x_max) + x_max
                else:
                    x_wheel_center = 0.5*(x_min - maxDomain) + maxDomain
                    x_min = maxDomain
                wheels.append(((x_wheel_center, HMAThickness, y_wheel_center),))
                if(abs(x_wheel_center) <= 0.5*Joint):
                    wheels_lines_Joint.append(((x_wheel_center, HMAThickness, y_min_l),))
                    wheels_lines_Joint.append(((x_wheel_center, HMAThickness, y_max_l),))
                elif(abs(x_wheel_center) <= 0.5*Joint + 0.5*LocalZoneRefSize):
                    wheels_lines_LocalProb.append(((x_wheel_center, HMAThickness, y_min_l),))
                    wheels_lines_LocalProb.append(((x_wheel_center, HMAThickness, y_max_l),))
                elif(abs(x_wheel_center) <= 0.5*Joint + 0.5*TransitionZone):
                    wheels_lines_TransZone.append(((x_wheel_center, HMAThickness, y_min_l),))
                    wheels_lines_TransZone.append(((x_wheel_center, HMAThickness, y_max_l),))
                else:
                    wheels_lines.append(((x_wheel_center, HMAThickness, y_min_l),))
                    wheels_lines.append(((x_wheel_center, HMAThickness, y_max_l),))
    wheels = tuple(wheels)
    slabs.Surface(name='Gears', side1Faces=
        slabs.faces.findAt(*wheels))
    wheels_lines = tuple(wheels_lines)
    if(len(wheels_lines) != 0):
        slabs.Set(edges=
            slabs.edges.findAt(*wheels_lines),  
            name='Gears_Lines')
    wheels_lines_Joint = tuple(wheels_lines_Joint)
    if(len(wheels_lines_Joint) != 0):
        slabs.Set(edges=
            slabs.edges.findAt(*wheels_lines_Joint),  
            name='Gears_Lines_Joint')
    wheels_lines_LocalProb = tuple(wheels_lines_LocalProb)
    if(len(wheels_lines_LocalProb) != 0):
        slabs.Set(edges=
            slabs.edges.findAt(*wheels_lines_LocalProb),  
            name='Gears_Lines_LocalProb')
    wheels_lines_TransZone = tuple(wheels_lines_TransZone)
    if(len(wheels_lines_TransZone) != 0):
        slabs.Set(edges=
            slabs.edges.findAt(*wheels_lines_TransZone),  
            name='Gears_Lines_TransZone')

# -------------------- Assembly and Step
pavement.rootAssembly.DatumCsysByDefault(CARTESIAN)
pavement.rootAssembly.Instance(dependent=OFF, name=
    'Slabs-1', part=slabs)
pavement.StaticStep(name='Static', previous='Initial')

# -------------------- Boundary conditions
if not isImposeDisp:
    pavement.Pressure(amplitude=UNSET, createStepName='Static'
        , distributionType=UNIFORM, field='', magnitude=Pressure, name='Airplane_Load', 
        region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Gears'])
else:
    pavement.Pressure(amplitude=UNSET, createStepName='Static'
        , distributionType=UNIFORM, field='', magnitude=1.0, name='Disp_Joint_F1', 
        region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Joint_Face_1'])
    pavement.Pressure(amplitude=UNSET, createStepName='Static'
        , distributionType=UNIFORM, field='', magnitude=1.0, name='Disp_Joint_F2', 
        region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Joint_Face_2'])

pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Base_BC_0', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Base_BC_0'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Base_BC_1', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Base_BC_1'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Base_BC_2', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Base_BC_2'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Base_BC_3', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Base_BC_3'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Concrete_BC_0', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Concrete_BC_0'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Concrete_BC_1', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Concrete_BC_1'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Concrete_BC_2', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Concrete_BC_2'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Concrete_BC_3', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Concrete_BC_3'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='HMA_BC_0', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['HMA_BC_0'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='HMA_BC_1', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['HMA_BC_1'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='HMA_BC_2', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['HMA_BC_2'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='HMA_BC_3', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['HMA_BC_3'])
pavement.Pressure(amplitude=UNSET, createStepName='Static'
    , distributionType=UNIFORM, field='', magnitude=1.0, name='Bottom_face_BC', 
    region=pavement.rootAssembly.instances['Slabs-1'].surfaces['Bottom_face_BC'])

# -------------------- Mesh
pavement.rootAssembly.setMeshControls(elemShape=TET, 
    regions=
    pavement.rootAssembly.instances['Slabs-1'].cells, technique=FREE)
pavement.rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
    elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
    regions=(
    pavement.rootAssembly.instances['Slabs-1'].cells, ))
pavement.rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=C3D8R, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
    elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, distortionControl=DEFAULT)), regions=(
    pavement.rootAssembly.instances['Slabs-1'].cells, ))

pavement.rootAssembly.setMeshControls(allowMapped=True, 
    regions=pavement.rootAssembly.instances['Slabs-1'].cells)

pavement.rootAssembly.setMeshControls(allowMapped=False, 
    regions=
    pavement.rootAssembly.instances['Slabs-1'].cells.findAt(
    (( 0.5*Joint + 0.25*LocalZoneRefSize, 0.5*HMAThickness, 0.0), ), 
    ((-0.5*Joint - 0.25*LocalZoneRefSize, 0.5*HMAThickness, 0.0), ),
    (( 0.5*(Joint + LocalZoneRefSize) + 0.25*(TransitionZone - LocalZoneRefSize), 0.5*HMAThickness, 0.0), ), 
    ((-0.5*(Joint + LocalZoneRefSize) - 0.25*(TransitionZone - LocalZoneRefSize), 0.5*HMAThickness, 0.0), ), 
    (( 0.5*Joint + 0.25*LocalZoneRefSize, -0.5*ConcreteThickness, 0.0), ), 
    ((-0.5*Joint - 0.25*LocalZoneRefSize, -0.5*ConcreteThickness, 0.0), ),
    (( 0.5*(Joint + LocalZoneRefSize) + 0.25*(TransitionZone - LocalZoneRefSize), -0.5*ConcreteThickness, 0.0), ), 
    ((-0.5*(Joint + LocalZoneRefSize) - 0.25*(TransitionZone - LocalZoneRefSize), -0.5*ConcreteThickness, 0.0), ), ))

pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
    deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Ext_Z_Lines'].edges, size=maxElemSize)
pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
    deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Int_Z_Lines'].edges, size=transitionZoneMaxSize)
pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
    deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Ext_Y_Lines'].edges, size=thicknessRefSize)
pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
    deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Ext_Local_Ref_Lines'].edges, size=localZoneMaxSize)
pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
    deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Base_RefZone_Lines'].edges, size=transitionZoneMaxSize)
pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
    deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Joint_Lines'].edges, size=Joint)

if not isImposeDisp:
    if(len(wheels_lines) != 0):
        pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Gears_Lines'].edges, size=transitionZoneMaxSize)
    if(len(wheels_lines_TransZone) != 0):
        pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Gears_Lines_TransZone'].edges, size=localZoneMaxSize)
    if(len(wheels_lines_LocalProb) != 0):
        pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Gears_Lines_LocalProb'].edges, size=Joint)
    if(len(wheels_lines_Joint) != 0):
        pavement.rootAssembly.seedEdgeBySize(constraint=FINER, 
            deviationFactor=0.1, edges=pavement.rootAssembly.instances['Slabs-1'].sets['Gears_Lines_Joint'].edges, size=Joint)

pavement.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FINER, end1Edges=pavement.rootAssembly.instances['Slabs-1'].sets['Ext_X_Lines_dir1'].edges, end2Edges=
    pavement.rootAssembly.instances['Slabs-1'].sets['Ext_X_Lines_dir2'].edges, maxSize=maxElemSize, minSize=transitionZoneMaxSize)
pavement.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FINER, end1Edges=pavement.rootAssembly.instances['Slabs-1'].sets['LocalRefZone_XLines_dir1'].edges, end2Edges=
    pavement.rootAssembly.instances['Slabs-1'].sets['LocalRefZone_XLines_dir2'].edges, maxSize=localZoneMaxSize, minSize=Joint)
pavement.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FINER, end1Edges=pavement.rootAssembly.instances['Slabs-1'].sets['TransZone_XLines_dir1'].edges, end2Edges=
    pavement.rootAssembly.instances['Slabs-1'].sets['TransZone_XLines_dir2'].edges, maxSize=transitionZoneMaxSize, minSize=localZoneMaxSize)
pavement.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FINER, end1Edges=pavement.rootAssembly.instances['Slabs-1'].sets['Joint_Y_dir1'].edges, end2Edges=
    pavement.rootAssembly.instances['Slabs-1'].sets['Joint_Y_dir2'].edges, maxSize=transitionZoneMaxSize, minSize=Joint)
pavement.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FINER, end1Edges=pavement.rootAssembly.instances['Slabs-1'].sets['LocalZone_Y_dir1'].edges, end2Edges=
    pavement.rootAssembly.instances['Slabs-1'].sets['LocalZone_Y_dir2'].edges, maxSize=transitionZoneMaxSize, minSize=localZoneMaxSize)

pavement.rootAssembly.generateMesh(regions=(
    pavement.rootAssembly.instances['Slabs-1'], ))

# -------------------- Defining job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=ModelName, modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=JobName, nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

mdb.jobs[JobName].writeInput()
