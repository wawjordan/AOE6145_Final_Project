#!./bin/bash
input="input.nml"
start=`date +%s`
flux_str=""
limiter_str=""
banner_str1="!========================================="
banner_str1+="=====================================!"
banner_str2="------------------------------------------"
banner_str2+="--------------------------------------"

#imax=49 97 193 385
#jmax=14 27  53 105
imax=193
jmax=53
grid_dir="../grids/NACA64A006-grids/"
#grid_name="NACA64A006.extra-coarse.27x14.grd"
#grid_name="NACA64A006.coarse.53x27.grd"
grid_name="NACA64A006.medium.193x53.grd"
#grid_name="NACA64A006.fine.385x105.grd"
C_grid="T"
index1=1
#index2= 8 16 32 64
index2=32
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

T_inf=300.0
#alpha=0.0
#p_inf=65855.8
#M_inf=0.84

alpha=8.0
p_inf=67243.5
M_inf=0.75

CFL=0.5
max_iter=200000

flux_scheme=2
limiter_scheme=0
beta_lim=2.0
eps_roe=0.1

out_dir="AIRFOIL/alpha2/Roe/"
out_file="NACA64A006_$((imax))x$((jmax))"
flux_out="F"
soln_save=200000
res_save=1
res_out=100

# airfoil
num_BCs=5
bounds=()
bounds+=(2,$i_low,$i_high,$((jg_high-n_ghost+1)),$jg_high,4)
bounds+=(3,$((ig_high-n_ghost+1)),$ig_high,$j_low,$j_high,2)
bounds+=(3,$ig_low,$((ig_low+n_ghost-1)),$j_low,$j_high,1)
bounds+=(5,$((index2+1)),$((imax-index2-1)),$jg_low,$((j_low-1)),3)
bounds+=(6,$ig_low,$index2,$jg_low,$((j_low-1)),3)

epsM=1.0
kappaM=-1.0
limiter_freeze="T"
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
echo "  C_grid = $C_grid" >> $input
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
echo "  alpha = $alpha" >> $input
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
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
