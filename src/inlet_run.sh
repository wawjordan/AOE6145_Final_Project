#!./bin/bash
input="input.nml"
start=`date +%s`
flux_str=""
limiter_str=""
banner_str1="!========================================="
banner_str1+="=====================================!"
banner_str2="------------------------------------------"
banner_str2+="--------------------------------------"
#imax=53,105,209,417
#jmax=17,33,65,129
imax=417
jmax=129
grid_dir="../grids/inlet-grids/"
#grid_name="Inlet.33x17.grd"
#grid_name="Inlet.105x33.grd"
#grid_name="Inlet.209x65.grd"
grid_name="Inlet.417x129.grd"
inlet="T"
index1=1
#index2=20,40,80,160
index2=160
n_ghost=2

i_low=1
j_low=1
i_high=$(( imax - 1 ))
j_high=$(( jmax - 1 ))
ig_low=$(( i_low - n_ghost ))
jg_low=$(( j_low - n_ghost ))
ig_high=$(( i_high + n_ghost ))
jg_high=$(( j_high + n_ghost ))

gamma=1.4

p_inf=12270.0
T_inf=217.0
M_inf=4.0

CFL=0.4
max_iter=100000

flux_scheme=2
limiter_scheme=1
beta_lim=2.0
eps_roe=0.1

out_dir="INLET/Roe/"
out_file="inlet_$((imax))x$((jmax))"
flux_out="F"
soln_save=10000
res_save=1
res_out=100

# inlet
num_BCs=4
bounds=()
bounds+=(2,$i_low,$index2,$jg_low,$((jg_low+1)),4)
bounds+=(4,$((ig_high-n_ghost+1)),$ig_high,$j_low,$j_high,2)
bounds+=(5,$i_low,$i_high,$jg_low,$((jg_low+n_ghost-1)),3)
bounds+=(5,$((index2+1)),$i_high,$((jg_high-n_ghost+1)),$jg_high,4)

epsM=1.0
kappaM=-1.0
limiter_freeze="F"
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
echo "| $out_file | "
echo "$banner_str2"
echo "| $flux_str | $limiter_str | CFL=$CFL "
echo "$banner_str1"

echo "&grid" > $input
echo "  grid_dir = \"$grid_dir\"" >> $input
echo "  grid_name = \"$grid_name\"" >> $input
echo "  inlet = $inlet" >> $input
echo "  index1 = $index1" >> $input
echo "  index2 = $index2" >> $input
echo "  imax = $imax" >> $input
echo "  jmax = $jmax" >> $input
echo "  n_ghost = $n_ghost" >> $input
echo "/" >> $input
echo "" >> $input
echo "&geometry" >> $input
echo "/" >> $input
echo "" >> $input
echo "&constants" >> $input
echo "  gamma = $gamma" >> $input
echo "  num_BCs = $num_BCs" >> $input
echo "/" >> $input
echo "" >> $input
echo "&initial" >> $input
echo "  p_inf = $p_inf" >> $input
echo "  T_inf = $T_inf" >> $input
echo "  M_inf = $M_inf" >> $input
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
#rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
