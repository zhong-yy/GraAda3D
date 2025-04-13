# generate observation points
generatePoints -x 0:50:2000 -y 0:50:2000 -z -0.1 -o test_points

compute -m ../test_make_model/test_model1.xyz -p test_points -c "gz/Txz/Tyz/Tzz" -o result -n 0.02

# plot 
python plot.py