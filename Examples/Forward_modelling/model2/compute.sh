# compute
mkdir -p result
compute -m model2.xyz -p obs_points -c "gz/Txz/Tyz/Tzz" -n 0.02 -o result/data 
awk 'NR>1 {print $1, $2, $4}' result/data > result/dobs_g_z
awk 'NR>1 {print $1, $2, $5, $6, $7}' result/data > result/dobs_Tzz_Txz_Tyz

# plot 
python plot.py