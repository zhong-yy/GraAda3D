1 GMT_slice_node_reg.sh: Slice 3D netcdf data, write out with node registration

2 GMT_slice_node_reg.sh: Slice 3D netcdf data, write out with pixel registration

For explation of the two ways of data organization, see https://docs.generic-mapping-tools.org/6.1/cookbook/options.html#option-nodereg

There is not much difference between the two ways of registration, but the pixel registration may better depict the inversion mesh, while the node registration may cause slight distortion. However, to use the pixel registration, extra caution is required for data slicing using `project` and `grdtrack`, because proper offsets should be added to the coordinates.
