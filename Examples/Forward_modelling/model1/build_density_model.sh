# Construct a density model in a region of x=[0, 2000]m, y=[0, 2000]m, z=[0, 2000]m (-R 0/2000/0/2000/0/1000)
# The region is divided into 40 cells in the x direction, 40 cells in the y direction and 20 cells in the z directions (--nx 40 --ny 40 --nz 20)
# There are 40*40*20 cells in total.
# Add 3 anomalous bodies (-a 600/900/550/850/200/500/500 -a 1200/1500/550/850/200/500/-500 -a 900/1500/1300/1500/200/500/500)
# The discretized model will be written to test_model1.xyz (-o test_model1)
mkdir -p model
makeModel -R 0/2000/0/2000/0/1000 --nx 40 --ny 40 --nz 20 -a 600/900/550/850/200/500/500 -a 1200/1500/550/850/200/500/-500 -a 900/1500/1300/1500/200/500/500 -o model/test_model1

# constraint model 1
# Add 1 anomalous bodies (-a 900/1500/1300/1500/200/500/200)
# Output file: constraint.xyz (-o constraint)
makeModel -R 0/2000/0/2000/0/1000 --nx 40 --ny 40 --nz 20  -a 900/1500/1300/1500/200/500/200 -o model/constraint

awk 'NR>2 {print $7, $8, $9, $10}' model/constraint.xyz > model/constraint_as_input_format.xyz
