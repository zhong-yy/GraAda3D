# compute
mkdir -p result
# use "compute -h" to see the help message
compute -m model/test_model.xyz -p observation_points/obs_points -c "gz/Txx/Txy/Txz/Tyy/Tyz/Tzz" -n 0.02 -o result/data

awk 'NR>1 {print $1, $2, $3, $4}' result/data > result/dobs_g_z
awk 'NR>1 {print $1, $2, $3, $5, $6, $7}' result/data > result/dobs_Tzz_Txz_Tyz
awk 'NR>1 {print $1, $2, $3, $5}' result/data > result/dobs_Tzz
awk 'NR>1 {print $1, $2, $3, $4,$8,$9,$6,$10,$7,$5}' result/data > result/dobs_gz_Txx_Tyx_Tzx_Tyy_Tyz_Tzz
