

Type the following command in the command line window to plot data:

````bash
bash GMT_slice_node_reg.sh
#or
bash GMT_slice_node_reg.sh
````



**GMT_slice_node_reg.sh**: Slice 3D netcdf data, write out with node registration

**GMT_slice_node_reg.sh**: Slice 3D netcdf data, write out with pixel registration



Data slicing is implemented mainly via `project` and `grdtrack`. I have tried the new GMT feature `grdinterpolate`, but it only worked well for cutting horizontal slices and failed to extract vertical slices.



For detailed explation of the two ways of data organization, see https://docs.generic-mapping-tools.org/6.1/cookbook/options.html#option-nodereg

**Notes**
here is not much difference between the two ways of registration, but the pixel registration may better depict the inversion mesh, while the node registration may cause slight distortion. However, to use the pixel registration, extra caution is required for data slicing using `project` and `grdtrack`, because proper offsets should be added to the coordinates.