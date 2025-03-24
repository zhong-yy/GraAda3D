../../makeModel -R 0/2000/0/2000/0/1000 --nx 40 --ny 40 --nz 20 -a 600/900/550/850/200/500/200 --a 1200/1500/550/850/200/500/-200 --a 900/1500/1300/1500/200/500/200 -o test_model1

# constraint model 1
../../makeModel -R 0/2000/0/2000/0/1000 --nx 40 --ny 40 --nz 20  --a --a 900/1500/1300/1500/200/500/10 -o constraint

awk 'NR>2 {print $7, $8, $9, $10}' constraint.xyz > constraint_as_input_format.xyz
