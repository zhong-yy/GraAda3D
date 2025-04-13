# compute
compute -m test_model1.xyz -p obs_points -c "gz/Txz/Tyz/Tzz" -n 0.02 -o result 
#The above command can also be written as:
#compute --model test_model1.xyz --observation obs_points --component "gz/Txz/Tyz/Tzz" --noise 0.02 --output result 

# plot 
python plot.py