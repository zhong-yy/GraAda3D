mkdir -p model

# Construct a density model in a region of x=[0, 2000]m, y=[0, 2000]m, z=[0, 2000]m (-R 0/2000/0/2000/0/1000)
# The region is divided into 40 cells in the x direction, 40 cells in the y direction and 20 cells in the z directions (--nx 40 --ny 40 --nz 20)
# There are 40*40*20 cells in total.
# The anomalies are given in file anomalies.txt (-A anomalies.txt)
# The discretized model will be written to test_model1.xyz (-o test_model1)
makeModel -R 0/2000/0/2000/0/1000 --nx 40 --ny 40 --nz 20 -A anomalies.txt -o model/test_model1

# Construct a constraint model
# Input file: constraint_geometry.txt (-A constraint_geometry.txt)
# Output file: constraint.xyz (-o constraint)
makeModel -R 0/2000/0/2000/0/1000 --nx 40 --ny 40 --nz 20  -A constraint_geometry.txt -o model/constraint

# rearrange the columns of the constraint model to change it into a format required by GraAda3D
awk 'NR>2 {print $7, $8, $9, $10}' model/constraint.xyz > model/constraint_as_input_format.xyz
