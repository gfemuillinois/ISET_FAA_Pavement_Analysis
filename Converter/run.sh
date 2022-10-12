#!/bin/bash

while getopts "f:" flag
do
	case "${flag}" in
		f) filename=${OPTARG};;
	esac
done

sed -i 's/Slabs-1.//g' ./$filename

tr -d '\r' <./$filename> ./"linux_${filename}"

rm ./$filename

mv ./"linux_${filename}" $filename

./Converter 1 $filename

grfFilename="${filename%.*}.grf"

rm ./OutTet4_IGL.inp

mv ./OutTet4.grf ./$grfFilename

sed -i 's/Traction Base_BC_0	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring Base_BC_0	{ Springs = 60.0e06 0.0 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Base_BC_1	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring Base_BC_1	{ Springs =  0.0 0.0 60.0e06 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Base_BC_2	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring Base_BC_2	{ Springs = 60.0e06 0.0 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Base_BC_3	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring Base_BC_3	{ Springs =  0.0 0.0 60.0e06 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Concrete_BC_0	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/#Spring Concrete_BC_0	{ Springs = 5520.0e06 0.0 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Concrete_BC_1	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/#Spring Concrete_BC_1	{ Springs = 0.0 0.0 5520.0e06 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Concrete_BC_2	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/#Spring Concrete_BC_2	{ Springs = 5520.0e06 0.0 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Concrete_BC_3	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/#Spring Concrete_BC_3	{ Springs = 0.0 0.0 5520.0e06 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction HMA_BC_0	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring HMA_BC_0	{ Springs = 2000.0e06 0.0 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction HMA_BC_1	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring HMA_BC_1	{ Springs =  0.0 0.0 2000.0e06 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction HMA_BC_2	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring HMA_BC_2	{ Springs = 2000.0e06 0.0 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction HMA_BC_3	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring HMA_BC_3	{ Springs = 0.0 0.0 2000.0e06 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Bottom_face_BC	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Spring Bottom_face_BC	{ Springs = 0.0 55.0e06 0.0 Displacements= 0 0 0}/g' ./$grfFilename
sed -i 's/Traction Joint_Face_1	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Displacement Joint_Face_1	{ Components = -0.0003 0.0 0.0 Flags= 1 0 0}/g' ./$grfFilename
sed -i 's/Traction Joint_Face_2	{ Components = -1 0. 0. CoordinateFrame = FaceAligned }/Displacement Joint_Face_2	{ Components = 0.0003 0.0 0.0 Flags= 1 0 0}/g' ./$grfFilename

sed -i '/Spring HMA_BC_3/a \ \
#GFEM-gl BCs\
Displacement displ_local { Components = 0.00 0.00 0.00 Flags = 1 1 1 }\
Traction trac_local { Components = 0.00 0.00 0.00 }\
Spring spr_local { Springs = 0.00 0.00 0.00 Displacements = 0.00 0.00 0.00 } ' ./$grfFilename

