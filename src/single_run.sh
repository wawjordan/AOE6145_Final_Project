#!./bin/bash
input="input.nml"
summary="../results/summary.dat"
temp="temp.txt"
start=`date +%s`
MMS_str=""
flux_str=""
limiter_str=""
test_str=""
#touch $summary
#if [ -f $summary ]; then
#  rm -f "$summary"
#fi
#grid_name="/home/grad3/wajordan/AOE_6145/AOE_6145_Final_project/grids"
#grid_name+="/curvilinear-grids/"
#grid_name+="curv2d17.grd"
grid_name="../grids/curvilinear-grids/"
grid_name+="curv2d33.grd"
cart_grid="T"
imax=33
jmax=33
n_ghost=2

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
isAxi="F"
Lmms=1.0

gamma=1.4

isMMS="T"
u0=-800.0
v0=0.0
u_inf=800.0
alpha=90.0
p_inf=100000.0
T_inf=500.0


CFL=0.01 #0.1 0.5 0.9
max_iter=800

flux_scheme=2
limiter_scheme=2
beta_lim=2.0
eps_roe=0.1

geometry_file="example.dat"
soln_save=8
res_save=1
res_out=1
cons="T"

epsM=0.0
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
echo "  isAxi = $isAxi" >> $input
echo "  Lmms = $Lmms" >> $input
echo "/" >> $input
echo "" >> $input
echo "&constants" >> $input
echo "  gamma = $gamma" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  isMMS = $isMMS" >> $input
echo "  u_inf = $u_inf" >> $input
echo "  alpha = $alpha" >> $input
echo "  p_inf = $p_inf" >> $input
echo "  T_inf = $T_inf" >> $input
echo "  u0 = $u0" >> $input
echo "  v0 = $v0" >> $input
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
#cat "$temp" >> $summary
#rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
