1. Construct a synthetic discretized model using `makeModel` (type `makeModel -h` to see the help message)

There are two ways to to specify anomalous bodies in `makeModel`. One is specifying anomalies in command line using -a option. Option -a can be used multiple times to specify multiple anomalies. This is convenient if there are only a few simple anomalous bodies. However, when there are a lot of anomalous bodies, or the geometry of the anomalous body is too complex, you may need to repeat -a option many times, which will reduce the readability of your command. An example is given in "./build_density_model.sh".

The other way to specify anomalies is reading information of anomalous bodies from a file with option -A. -A can be used only once. -A and -a cannot be used at the same time. This is the recommended way when there are many anomalous bodies. An example is given in "build_density_model2.sh". 

To run the example, type bash "./build_density_model.sh" or bash "./build_density_model2.sh".

2. Generate observation points (type `generatePoints -h` to see the help message)

bash "./generate_observation_points.sh"

3. Compute gravity and gravity gradient, and plot the results ((type `compute -h` to see the help message))

bash "./compute.sh"
