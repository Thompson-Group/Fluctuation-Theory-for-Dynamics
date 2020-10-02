#!/bin/bash

#------------------------------------
#- This is the read input file      -
#-   Direct Fluctions Calculations  -
#- Only change uncommented lines    -
#------------------------------------


# Program (Lammps or CP2K):
program=`sed -n '7p' ../input_file`
# Run Style (NVT or NPT):
run_style=`sed -n '9p' ../input_file`
# Number of Molecules
num_molecs=`sed -n '11p' ../input_file`
# Identifier of Molecules
molec1=`sed -n '13p' ../input_file`
molec2=`sed -n '14p' ../input_file`
molec3=`sed -n '15p' ../input_file`
molec4=`sed -n '16p' ../input_file`
# Starting Configuration (in steps)
start_config=`sed -n '18p' ../input_file`
# Ending Configuration (in steps)
end_config=`sed -n '20p' ../input_file`
# Separation of Configurations (in steps)
sep_config=`sed -n '22p' ../input_file`
# Timestep (in ps)
timestep=`sed -n '24p' ../input_file`
# Num Times
num_times=`sed -n '26p' ../input_file`
# Blocks
blocks=`sed -n '28p' ../input_file`
# Segments
segsplit=`sed -n '30p' ../input_file`
# Files
num_files=`sed -n '32p' ../input_file`
# NVE Length
nve_length=`sed -n '34p' ../input_file` 
# Num RPJ
num_rpj=`sed -n '36p' ../input_file`
# TYPE
calctype=`sed -n '38p' ../input_file`
# Constraint
constraint=`sed -n '40p' ../input_file`

echo "Input File Parameters"
echo "Program: "$program
echo "Run Style: "$run_style
echo "Num Molecs: "$num_molecs
echo "mol1: "$molec1
echo "mol2: "$molec2
echo "mol3: "$molec3
echo "mol4: "$molec4
echo "Start Config: "$start_config
echo "End Config: "$end_config
echo "Sep Config: "$sep_config
echo "Timestep: "$timestep
echo "Num Times: "$num_times
echo "Blocks: "$blocks
echo "Segments: "$segsplit
echo "Num Files: "$num_files
echo "NVE Length: "$nve_length
echo "NUM RPJ: "$num_rpj
echo "Calc Type: "$calctype
echo "Constraint: "$constraint
echo "End Input File Parameters"
