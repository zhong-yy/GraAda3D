1. Construct a synthetic discretized model using `makeModel` (type `makeModel -h` to see the help message)

bash "./build_density_model.sh"

2. Generate observation points (type `generatePoints -h` to see the help message)

bash "./generate_observation_points.sh"

3. Compute gravity and gravity gradient, and plot the results ((type `compute -h` to see the help message))

bash "./compute.sh"

4. Plot the fields
 
```
python plot.py
python plot_true_model_slice.py
python plot_ref_model_slice.py
python plot_crg_ref_model_slice.py
```

Note: please open and read the .sh files before running them 
