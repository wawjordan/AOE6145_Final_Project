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
imax=97
jmax=27
#grid_dir="../grids/curvilinear-grids/"
#grid_name="curv2d$imax"
#grid_name+=".grd"
#grid_name="../grids/inlet-grids/"
#grid_name+="Inlet.33x17.grd"
grid_dir="../grids/NACA64A006-grids/"
#grid_name="NACA64A006.extra-coarse.27x14.grd"
grid_name="NACA64A006.coarse.53x27.grd"
#grid_name="NACA64A006.fine.385x105.grd"
cart_grid="F"
C_grid="T"
index1=1
index2=16
#index2=64
n_ghost=2

i_low=1
j_low=1
i_high=$(( imax - 1 ))
j_high=$(( jmax - 1 ))
ig_low=$(( i_low - n_ghost ))
jg_low=$(( j_low - n_ghost ))
ig_high=$(( i_high + n_ghost ))
jg_high=$(( j_high + n_ghost ))

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
isAxi="F"
Lmms=1.0

gamma=1.4

isMMS="F"
u0=20.0
v0=0.0
u_inf=20.0 #1.0
alpha=0.0
#alpha=8.0
p_inf=65855.8
#p_inf=67243.5
T_inf=300.0
M_inf=0.84
#M_inf=0.75


CFL=0.1 #0.1 0.5 0.9
max_iter=200000
locTime="F"

flux_scheme=1
limiter_scheme=1
beta_lim=2.0
eps_roe=0.01

geometry_file="example.dat"
flux_out="F"
soln_save=50000
res_save=1
res_out=100
cons="T"

num_BCs=5
bounds=()
bounds+=(2,$i_low,$i_high,$((jg_high-n_ghost+1)),$jg_high,4)
bounds+=(3,$((ig_high-n_ghost+1)),$ig_high,$j_low,$j_high,2)
bounds+=(3,$ig_low,$((ig_low+n_ghost-1)),$j_low,$j_high,1)
bounds+=(5,$((index2+1)),$((imax-index2-1)),$jg_low,$((j_low-1)),3)
bounds+=(6,$ig_low,$index2,$jg_low,$((j_low-1)),3)

epsM=1.0
kappaM=-1.0
limiter_freeze="F"
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
echo "  grid_dir = \"$grid_dir\"" >> $input
echo "  grid_name = \"$grid_name\"" >> $input
echo "  cart_grid = $cart_grid" >> $input
echo "  C_grid = $C_grid" >> $input
echo "  index1 = $index1" >> $input
echo "  index2 = $index2" >> $input
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
echo "  num_BCs = $num_BCs" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  isMMS = $isMMS" >> $input
echo "  u_inf = $u_inf" >> $input
echo "  alpha = $alpha" >> $input
echo "  p_inf = $p_inf" >> $input
echo "  T_inf = $T_inf" >> $input
echo "  M_inf = $M_inf" >> $input
echo "  u0 = $u0" >> $input
echo "  v0 = $v0" >> $input
echo "/" >> $input
echo "" >> $input
echo "&numerical" >> $input
echo "  CFL = $CFL" >> $input
echo "  max_iter = $max_iter" >> $input
echo "  locTime = $locTime" >> $input
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
echo "  soln_save = $soln_save" >> $input
echo "  res_save = $res_save" >> $input
echo "  res_out = $res_out" >> $input
echo "  cons = $cons" >> $input
echo "  flux_out = $flux_out" >> $input
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
