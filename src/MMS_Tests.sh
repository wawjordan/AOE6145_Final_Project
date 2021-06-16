#!./bin/bash
input="input.nml"
summary="../results/summary.dat"
temp="temp.txt"
start=`date +%s`
MMS_str=""
flux_str=""
limiter_str=""
test_str=""
touch $summary
if [ -f $summary ]; then
  rm -f "$summary"
fi
for imax in 9 17 33 65 129 257
do
grid_name="../grids/curvilinear-grids/"
grid_name+="curv2d$imax"
grid_name+=".grd"
cart_grid="T"
jmax=$imax
n_ghost=2

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
Lmms=1.0

gamma=1.4

CFL=0.1
max_iter=10000

flux_scheme=1
limiter_scheme=0
beta_lim=2.0
eps_roe=0.1

geometry_file="example.dat"
soln_save=$max_iter
res_save=1
res_out=100
cons="T"

epsM=1.0
kappaM=-1.0
limiter_freeze="T"
if [ "$isMMS" = "T" ]; then
  MMS_str="MMS"
else
  MMS_str="$test_str"
fi
if [ $flux_scheme -eq 1 ]; then
  flux_str="van Leer flux"
elif [ $flux_scheme -eq 2 ]; then
  flux_str="Roe's flux"
fi

if [ $limiter_scheme -eq 1 ]; then
  limiter_str="van Leer limiter"
elif [ $limiter_scheme -eq 2 ]; then
  limiter_str="van Albada limiter"
elif [ $limiter_scheme -eq 3 ]; then
  limiter_str="minmod limiter"
elif [ $limiter_scheme -eq 4 ]; then
  limiter_str="beta ($beta_lim) limiter"
fi
echo ""
echo "!==============================================================================!"
echo "| N=$imax | $MMS_str "
echo "--------------------------------------------------------------------------------"
echo "| $flux_str | $limiter_str | CFL=$CFL "
echo "!==============================================================================!"

echo "&grid" > $input
echo "  grid_name = \"$grid_name\"" >> $input
echo "  cart_grid = $cart_grid" >> $input
echo "  C_grid = F" >> $input
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
echo "  isAxi = F" >> $input
echo "  Lmms = $Lmms" >> $input
echo "/" >> $input
echo "" >> $input
echo "&constants" >> $input
echo "  gamma = $gamma" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  isMMS = T" >> $input
echo "/" >> $input
echo "" >> $input
echo "&numerical" >> $input
echo "  CFL = $CFL" >> $input
echo "  max_iter = $max_iter" >> $input
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
echo "  soln_save = $soln_save" >> $input
echo "  res_save = $res_save" >> $input
echo "  res_out = $res_out" >> $input
echo "  cons = $cons" >> $input
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
