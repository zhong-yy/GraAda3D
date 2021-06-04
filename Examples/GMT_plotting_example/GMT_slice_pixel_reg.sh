function check_and_delete(){
    if [ -f "$1" ];then
        echo "delete $1"
        rm $1
    fi
}
# Using pixel registration


data_file=../Inv_gz/gz_result.nc

## Generate evenly spaced points along a profile. -N is necessary for the Cartesian coordinates
gmt project -C25/700 -E1975/700 -G50 -N -V > line1
gmt project -C25/1450 -E1975/1450 -G50 -N -V > line2


gmt math -o0 -T25/1975/50 --FORMAT_FLOAT_OUT=%10.2f T = zpoints.txt

## Slicing 
check_and_delete cross_section1.txt
check_and_delete cross_section2.txt
while read z; do
	gmt grdtrack line1 -G"${data_file}?density(${z})" | gawk -v PREC=100 -v dep=${z} '{print $3+25, dep, $4}' >> cross_section1.txt
	gmt grdtrack line2 -G"${data_file}?density(${z})" | gawk -v PREC=100 -v dep=${z} '{print $3+25, dep, $4}' >> cross_section2.txt
done < zpoints.txt

## xyz2grd convert ascii table to binary grid file
gmt xyz2grd cross_section1.txt -R0/2000/0/1000 -I40+n/20+n -rp -V -Gcross_section1.nc
gmt xyz2grd cross_section2.txt -R0/2000/0/1000 -I40+n/20+n -rp -V -Gcross_section2.nc

## Alternatively, use "gmt surface" to interpolate data
# gmt surface cross_section.txt `gmt gmtinfo cross_section.txt -I-` -I50/50 -V -Gcross_section.nc

gmt begin 1pix jpg E300
    gmt set FONT 9p
    gmt set MAP_FRAME_PEN 0p,black
    gmt basemap -Jx0.002/-0.002 -R0/2000/0/2000 -Bxa500fg+l"x (m)" -Bya500fg+l"y (m)" -BWSrt
    gmt grdimage "${data_file}?density(350)"  -nn -E300 -Q -Cjet
    echo "z=350 m" | gmt text -F+cTL+f9p, -D0.c/0.6c -N
    gmt colorbar -DJMR+v+w3c/0.18c+o0.9c/0c -Bxaf -By+L"kg/m@+3@+"
gmt end show

gmt begin 2pix jpg E300
    gmt set FONT 9p
    gmt set MAP_FRAME_PEN 0p,black
    gmt basemap -Jx0.002/-0.002 -R0/2000/0/1000 -Bxa500fg+l"Distance (m)" -Bya500fg+l"Depth (m)" -BWSrt
    gmt grdimage cross_section1.nc  -nn -E300 -Q -Cjet
    echo "y=700 m" | gmt text -F+cTL+f9p, -D0.c/0.6c -N
    gmt colorbar -DJMR+v+w3c/0.18c+o0.9c/0c -Bxaf -By+L"kg/m@+3@+"
gmt end show

gmt begin 3pix jpg E300
    gmt set FONT 9p
    gmt set MAP_FRAME_PEN 0p,black
    gmt basemap -Jx0.002/-0.002 -R0/2000/0/1000 -Bxa500fg+l"Distance (m)" -Bya500fg+l"Depth (m)" -BWSrt
    gmt grdimage cross_section2.nc  -nn -E300 -Q -Cjet
    echo "y=1450 m" | gmt text -F+cTL+f9p, -D0.c/0.6c -N
    gmt colorbar -DJMR+v+w3c/0.18c+o0.9c/0c -Bxaf -By+L"kg/m@+3@+"
gmt end show

gmt begin slices_pix jpg E300
    gmt set FONT 9p
    gmt set MAP_FRAME_PEN 1p,black
#    gmt grd2cpt "${data_file}?density(350)" cross_section1.nc cross_section2.nc -Crainbow
#    gmt grd2cpt "${data_file}?density(350)" cross_section1.nc cross_section2.nc -Crainbow
    gmt makecpt  -T-200/200 -D -Cjet
    gmt basemap -Jx0.002/-0.002 -R0/2000/0/2000 -Bxa500fg+l"x (m)" -Bya500fg+l"y (m)" -BWSrt
    gmt grdimage "${data_file}?density(350)"  -nn -E300 -Q
    echo "(a) z=350 m" | gmt text -F+cTL+f9p, -D0.c/0.6c -N
    gmt colorbar -DJBC+h+w3.9c/0.18c+o0.c/1.25c+e -Bxaf+L"Density contrast (kg/m@+3@+)" --MAP_FRAME_PEN=0.1p,black # -By+L"kg/m@+3@+"
    
    gmt basemap -Jx0.002/-0.002 -R0/2000/0/1000 -Bxa500fg+l"y (m)" -Bya500fg+l"Depth (m)" -BWSrt -X6.5c -Y2.8c
    gmt grdimage cross_section1.nc  -nn -E300 -Q
    echo "(b) y=700 m" | gmt text -F+cTL+f9p, -D0.c/0.6c -N
    #gmt colorbar -DJMR+v+w3c/0.18c+o0.9c/0c+e -Bxaf -By+L"kg/m@+3@+"    
    
    gmt basemap -Jx0.002/-0.002 -R0/2000/0/1000 -Bxa500fg+l"y (m)" -Bya500fg+l"Depth (m)" -BWSrt -Y-4.5c
    gmt grdimage cross_section2.nc  -nn -E300 -Q
    echo "(c) y=1450 m" | gmt text -F+cTL+f9p, -D0.c/0.6c -N
    #gmt colorbar -DJMR+v+w3c/0.18c+o0.9c/0c+e -Bxaf -By+L"kg/m@+3@+"
gmt end show
