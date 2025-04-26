# compute
mkdir -p result
compute -m model/test_model2.xyz -p observation_points/obs_points -c "gz/Txz/Tyz/Tzz" -n 0.02 -o result/data 
awk 'NR>1 {print $1, $2, $4}' result/data > result/dobs_g_z
awk 'NR>1 {print $1, $2, $5, $6, $7}' result/data > result/dobs_Tzz_Txz_Tyz
awk 'NR>1 {print $1, $2, $5}' result/data > result/dobs_Tzz
awk 'NR>1 {print $1, $2, $4, $5}' result/data > result/dobs_gz_Tzz

# plot 
python plot.py