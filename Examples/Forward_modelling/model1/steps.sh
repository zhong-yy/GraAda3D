bash "./build_density_model.sh"
bash "./generate_observation_points.sh"
bash "./compute.sh"
python plot.py
python plot_test_model_slice.py