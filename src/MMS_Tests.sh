#!./bin/bash
input="input.nml"
summary="../results/MMS_subsonic/vanLeer/summary.dat"
temp="temp.txt"
start=`date +%s`
flux_str=""
limiter_str=""
banner_str1="!========================================="
banner_str1+="=====================================!"
banner_str2="------------------------------------------"
banner_str2+="--------------------------------------"
touch $summary
if [ -f $summary ]; then
  rm -f "$summary"
fi

out_dir="MMS_subsonic/vanLeer/"
gamma=1.4

supersonic="F"
isMMS="T"
cart_grid="F"
flux_out="F"
cons="T"
flux_out="F"

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
Lmms=1.0

n_ghost=2

flux_scheme=1
limiter_scheme=0
beta_lim=2.0
eps_roe=0.01

epsM=1.0
kappaM=-1.0
limiter_freeze="T"

CFL=0.2
max_iter=50000

soln_save=$max_iter
res_save=1
res_out=100


for imax in 9 17 33 65 129 257
do
jmax=$imax
grid_dir="../grids/curvilinear-grids/"
grid_name="curv2d$imax"
grid_name+=".grd"
out_file="MMS_$((imax))x$((jmax))"

i_low=1
j_low=1
i_high=$(( imax - 1 ))
j_high=$(( jmax - 1 ))
ig_low=$(( i_low - n_ghost ))
jg_low=$(( j_low - n_ghost ))
ig_high=$(( i_high + n_ghost ))
jg_high=$(( j_high + n_ghost ))

num_BCs=4
bounds=()
if [ $supersonic = "T" ]; then
  if [ $cart_grid = "T" ]; then
    bounds+=(1,$ig_low,$((i_low-1)),$j_low,$j_high,1)
    bounds+=(4,$((i_high+1)),$ig_high,$j_low,$j_high,2)
    bounds+=(1,$i_low,$i_high,$jg_low,$((j_low-1)),3)
    bounds+=(4,$i_low,$i_high,$((j_high+1)),$jg_high,4)
  else
    bounds+=(4,$ig_low,$((i_low-1)),$j_low,$j_high,1)
    bounds+=(1,$((i_high+1)),$ig_high,$j_low,$j_high,2)
    bounds+=(4,$i_low,$i_high,$jg_low,$((j_low-1)),3)
    bounds+=(1,$i_low,$i_high,$((j_high+1)),$jg_high,4)
  fi
else
    bounds+=(1,$ig_low,$((i_low-1)),$j_low,$j_high,1)
    bounds+=(1,$((i_high+1)),$ig_high,$j_low,$j_high,2)
    bounds+=(1,$i_low,$i_high,$jg_low,$((j_low-1)),3)
    bounds+=(1,$i_low,$i_high,$((j_high+1)),$jg_high,4)
fi


if [ $flux_scheme -eq 1 ]; then
  flux_str="van Leer flux"
elif [ $flux_scheme -eq 2 ]; then
  flux_str="Roe's flux"
fi

if [ $limiter_scheme -eq 0 ]; then
  limiter_str="no limiter"
elif [ $limiter_scheme -eq 1 ]; then
  limiter_str="van Leer limiter"
elif [ $limiter_scheme -eq 2 ]; then
  limiter_str="van Albada limiter"
elif [ $limiter_scheme -eq 3 ]; then
  limiter_str="minmod limiter"
elif [ $limiter_scheme -eq 4 ]; then
  limiter_str="beta ($beta_lim) limiter"
fi
echo ""
echo "$banner_str1"
echo "| $out_file "
echo "$banner_str2"
echo "| $flux_str | $limiter_str | CFL=$CFL "
echo "$banner_str1"

echo "&grid" > $input
echo "  grid_dir = \"$grid_dir\"" >> $input
echo "  grid_name = \"$grid_name\"" >> $input
echo "  cart_grid = $cart_grid" >> $input
echo "  imax = $imax" >> $input
echo "  jmax = $jmax" >> $input
echo "  n_ghost = $n_ghost" >> $input
echo "/" >> $input
echo "" >> $input
echo "&geometry" >> $input
echo "  xmin =  $xmin" >> $input
echo "  xmax =  $xmax" >> $input
echo "  ymin =  $ymin" >> $input
echo "  ymax =  $ymax" >> $input
echo "  Lmms = $Lmms" >> $input
echo "/" >> $input
echo "" >> $input
echo "&constants" >> $input
echo "  gamma = $gamma" >> $input
echo "  num_BCs = $num_BCs" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  isMMS = $isMMS" >> $input
echo "/" >> $input
echo "" >> $input
echo "&numerical" >> $input
echo "  CFL = $CFL" >> $input
echo "  max_iter = $max_iter" >> $input
echo "/" >> $input
echo "" >> $input
echo "&boundary" >> $input
for ((i = 0; i < $num_BCs; i++)); do
echo "  bounds(:,$((i+1))) = ${bounds[$i]}" >> $input
done
echo "/" >> $input
echo "" >> $input
echo "&flux" >> $input
echo "  flux_scheme = $flux_scheme" >> $input
echo "  limiter_scheme = $limiter_scheme" >> $input
echo "  beta_lim = $beta_lim" >> $input
echo "/" >> $input
echo "" >> $input
echo "&output" >> $input
echo "  geometry_file = \"$geometry_file\"" >> $input
echo "  out_dir = \"$out_dir\"" >> $input
echo "  out_file = \"$out_file\"" >> $input
echo "  soln_save = $soln_save" >> $input
echo "  res_save = $res_save" >> $input
echo "  res_out = $res_out" >> $input
echo "  flux_out = $flux_out" >> $input
echo "/" >> $input
echo "" >> $input
echo "&reconstruction" >> $input
echo "  epsM = $epsM" >> $input
echo "  kappaM = $kappaM" >> $input
echo "  limiter_freeze = $limiter_freeze" >> $input
echo "/" >> $input

./../build/bin/test_program
cat "$temp" >> $summary
done
rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
