# Project: FAA Reflective Cracking
# Development: UIUC team

# This Python script is used to run SIMPLIFIED pavement simulations in which
# the SIFs are assumed as constant along the crack front and therefore
# all the crack vertices move the same DeltaAMax, which is a user-specified parameter.

# First STEP:
# Solve a linear STATIC simulation to find the vertex
# in which the maximum Energy Release Rate (ERR) happens.

import os

print("******************************")
print("Starting 1st STEP")
print("****************************** \n")

print("Running the linear static simulation...")

if os.path.exists("./static.out"):
    os.remove("./static.out")

os.system("iset pavement_GFEMgl_static.tcl >> static.out")

print ("\nStatic simulation finished with success!\n")

# Second STEP:
# Reading .cm_sif file and finding vertex with Gmax

print("******************************")
print("Starting 2nd STEP")
print("****************************** \n")

print("Finding vertex with maximum G...")

file = open("pavement_GFEMgl_static.cm_sif", "r")

for iline in range(22):
    file.readline()

coords = dict()

nvertices = 0
end = True
while end:
    line = file.readline()
    if line == "\n":
        end = False
    else:
        # Storing also coordinates
        # That's needed when generating the bounding box
        split = line.split(",")
        coords[nvertices] = [float(split[k]) for k in (1, 2, 3)]

        nvertices += 1

nvertices_mls = ((nvertices - 1) * 3) + 1
for iline in range(3):
    file.readline()

Gs = dict()
count_id = 0
for iline in range(nvertices_mls):
    split = file.readline().split(",")

    index = int(split[0]); err = float(split[7])
    if ( index % 3 == 0 ):
        Gs[count_id] = err
        count_id += 1

file.close()
id_max = max(Gs, key=Gs.get)

xc = coords[id_max][0]
yc = coords[id_max][1]
zc = coords[id_max][2]

print ("\nCode found vertex with max. ERR!")
print ("Vertex ID where max. ERR happens = " + str(id_max) + "\n")

# Third STEP:
# Run the crack propagation simulation

print("******************************")
print("Starting 3rd STEP")
print("****************************** \n")

print("Running crack propagation simulation...")

cmd = "iset pavement_GFEMgl_propagation_constant_vel.tcl " + str(id_max) + " " + str(xc) + " " \
                                                                               + str(yc) + " " \
                                                                               + str(zc) + " >> propagation.out"

if os.path.exists("./propagation.out"):
    os.remove("./propagation.out")

os.system(cmd)
