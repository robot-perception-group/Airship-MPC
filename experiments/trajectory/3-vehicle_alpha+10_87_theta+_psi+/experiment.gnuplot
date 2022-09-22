load "../../trajectory_experiment.gnuplot"
alphastate="enabled"
phistate="enabled"
thetastate="enabled"
rho=87*pi/180


#set terminal qt enhanced size 1024,768
#set terminal qt enhanced size 800,300 font "Helvetica"
set terminal svg size 1600,600 font "Helvetica,22"
#set terminal qt size 1600,600 font "Helvetica,23"
set tmargin 0
set lmargin 0
set rmargin 0
set bmargin 0
set view equal xyz
set view 0,0,1.7
set grid vertical xtics ytics ztics
set xyplane at 5
set datafile separator ","
set isosamples 50,200
set parametric
set arrow 1 from 0,0,0 to 0,-6,0 lc 3 lw 2
set urange [0:2*pi]
set vrange [0:100]
set xrange [24:-24]
set yrange [24:-24]
set zrange [:5] reverse


#set size 0.33,0.5
#set origin 0,0.5
#set view 30,60,1.7
#splot f3x(u,v,0),f3y(u,v,0),f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", f4x(u,v,wind),f4y(u,v,wind),f4z(u,v,wind) with lines lc 1 title "valid trajectories moving subject", persx(u,v),persy(u,v),persz(u,v) with lines lc 1 title "subject", "trajectory.csv" using ($2):($1):(-$3) with lines  lc 2 lw 2 title "MPC trajectory", "" using 8:7:(-$9) lc 7 with lines title "sensor pointing error in m"
#pause -1

set multiplot
set key off
#set size 0.33,0.5
#set origin 0.022,0.5
#set view 60,120,1.4
##set xlabel "Y (m)" rotate parallel offset graph 0.0,0.1,0
#set xlabel "Y (m)" rotate parallel offset -9,1.0,0
##set ylabel "X (m)" rotate parallel offset graph -0.1,-0.1,0
#set ylabel "X (m)" rotate parallel offset 3,2.0,0 
##set zlabel "-Z (m)" offset graph 0.3,0,0.8
#set zlabel "Z (m)" offset 5,5,0
#splot f3x(u,v,0),-f3y(u,v,0),-f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", \
#f4x(u,v,wind),-f4y(u,v,wind),-f4z(u,v,wind) with lines lc 1 dt 2 title "valid trajectories moving subject", \
#persx(u,v),-persy(u,v),-persz(u,v) with lines lc 1 title "subject", \
#"trajectory.csv" using ($2):($1):($3) every 3::0 with lines lc 2 lw 0.5 title "MPC trajectory", \
#"" using ($2):($1):($3) every 3::1 with lines lc 8 lw 0.5 title "MPC trajectory", \
#"" using ($2):($1):($3) every 3::2 with lines lc 6 lw 0.5 title "MPC trajectory", \
#"" using 8:7:($9) every 3::0 lc 7 with lines title "sensor pointing error in m", \
#"" using 8:7:($9) every 3::1 lc 7 with lines title "sensor pointing error in m", \
#"" using 8:7:($9) every 3::2 lc 7 with lines title "sensor pointing error in m"
unset ztics
unset zlabel
set size 0.36,0.67
set origin 0.32,0.30
set view 359.999,0,2.5
set xlabel "Y (m)" rotate parallel offset graph 0,-2.04
set ylabel "X (m)" rotate parallel offset graph 0.95,0
set xtics offset graph 0,-1.24
set ytics offset graph 0.08,0
splot f3x(u,v,0),-f3y(u,v,0),-f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", \
f4x(u,v,wind),-f4y(u,v,wind),-f4z(u,v,wind) with lines lc 1 dt 2 title "valid trajectories moving subject", \
persx(u,v),-persy(u,v),-persz(u,v) with lines lc 3 title "subject", \
"trajectory.csv" using ($2):($1):($3) every 3::0 with lines lc 2 lw 0.5 title "MPC trajectory", \
"" using ($2):($1):($3) every 3::1 with lines lc 8 lw 0.5 title "MPC trajectory", \
"" using ($2):($1):($3) every 3::2 with lines lc 6 lw 0.5 title "MPC trajectory"

#, \
#"" using 8:7:($9) every 3::0 lc 7 with lines title "sensor pointing error in m", \
#"" using 8:7:($9) every 3::1 lc 7 with lines title "sensor pointing error in m", \
#"" using 8:7:($9) every 3::2 lc 7 with lines title "sensor pointing error in m"
#set ztics in
#unset ytics
#unset ylabel
#set xlabel "Y (m)" norotate offset graph 0.68,0,1.0
#set xtics offset graph 0.04,0,.28
#set zlabel "Z (m)" rotate parallel offset graph 0.13,0
#set ztics offset graph 0.05,0
#set size 0.33,0.8
#set origin 0.636,0.45
#set view 90,0,1.1
#splot f3x(u,v,0),-f3y(u,v,0),-f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", \
#f4x(u,v,wind),-f4y(u,v,wind),-f4z(u,v,wind) with lines lc 1 dt 2 title "valid trajectories moving subject", \
#persx(u,v),-persy(u,v),-persz(u,v) with lines lc 1 title "subject", \
#"trajectory.csv" using ($2):($1):($3) every 3::0 with lines lc 2 lw 0.5 title "MPC trajectory", \
#"" using ($2):($1):($3) every 3::1 with lines lc 8 lw 0.5 title "MPC trajectory", \
#"" using ($2):($1):($3) every 3::2 with lines lc 6 lw 0.5 title "MPC trajectory", \
#"" using 8:7:($9) every 3::0 lc 7 with lines title "sensor pointing error in m", \
#"" using 8:7:($9) every 3::1 lc 7 with lines title "sensor pointing error in m", \
#"" using 8:7:($9) every 3::2 lc 7 with lines title "sensor pointing error in m"
set ytics

set ylabel rotate parallel offset graph 0.1,0
set ytics offset graph 0.02,0
set xtics offset graph 0,0.1
set size 0.36,0.32
set origin -0.05,0.68
set autoscale x
set autoscale y
set yrange [-22.5:22.5]
set samples 500
set trange [0:166]
set format x ''
set grid xtics
unset arrow 1
set ylabel "X (m)"
unset xlabel
unset lmargin
set format x ''
plot "trajectory.csv" using ($0*0.25):($1) every 3::0 with lines  lc 2 lw 0.5 title "MPC trajectory X", \
"" using ($0*0.25):($1) every 3::1 with lines lc 8 lw 0.5 title "MPC trajectory X", \
"" using ($0*0.25):($1) every 3::2 with lines lc 6 lw 0.5 title "MPC trajectory X", \
t,-f3y(-pi/2+t*psidot,reference_radius,wind)  with lines lc 1 dt 2 title "valid trajectory moving subject X" axis x1y1

#set xlabel "Time (s)" offset 0,0.5
set size 0.36,0.32
set origin -0.05,0.345
set x2label "Time (s)" offset graph 0,-0.3
set ylabel "Y (m)"

plot "trajectory.csv" using ($0*0.25):($2) every 3::0 with lines  lc 2 lw 0.5 title "MPC trajectory Y", \
"" using ($0*0.25):($2) every 3::1 with lines  lc 8 lw 0.5 title "MPC trajectory Y", \
"" using ($0*0.25):($2) every 3::2 with lines  lc 6 lw 0.5 title "MPC trajectory Y", \
t,f3x(-pi/2+t*psidot,reference_radius,wind)  with lines lc 1 dt 2 title "valid trajectory moving subject Y" axes x1y1

#set grid vertical xtics y2tics
set xtics offset graph 0.04,0.3
set size 0.36,0.32
#set origin -0.05,0.05
set origin 0.60,0.68
set ylabel "Z (m)"
set ylabel  offset graph 0.06,0
#set key at screen 0.98,0.5
set yrange [0:-45]

plot "trajectory.csv" using ($0*0.25):($3) every 3::0 axes x1y2 with lines  lc 2 lw 0.5 title "MPC trajectory Z", \
"" using ($0*0.25):($3) every 3::1 axes x1y2 with lines  lc 8 lw 0.5 title "MPC trajectory Z", \
"" using ($0*0.25):($3) every 3::2 axes x1y2 with lines  lc 6 lw 0.5 title "MPC trajectory Z", \
t,-f3z(-pi/2+t*psidot,reference_radius,wind)  with lines lc 1 dt 2 title "valid trajectory moving subject Z" axes x2y2

set size 0.36,0.32
set origin 0.60,0.345
set yrange [0:15]
set ylabel "Cost - {𝜓፞} (°/s)"
set ylabel  offset graph 0.06,0
plot "trajectory.csv" using ($0*0.25):($13) every 3::0 with lines lc 2 lw 0.5 title "{𝜓̇}" axes x1y1, \
"" using ($0*0.25):($13) every 3::1 with lines lc 8 lw 0.5 title "{𝜓፞}" axes x1y1, \
"" using ($0*0.25):($13) every 3::2 with lines lc 6 lw 0.5 title "{𝜓፞}" axes x1y1, \
"" using ($0*0.25):(sqrt($7**2+$8**2+$9**2)) every 3::0 with lines lc 7 title "sensor pointing error in m" axes x1y1, \
"" using ($0*0.25):(sqrt($7**2+$8**2+$9**2)) every 3::1 with lines lc 7 title "sensor pointing error in m" axes x1y1, \
"" using ($0*0.25):(sqrt($7**2+$8**2+$9**2)) every 3::2 with lines lc 7 title "sensor pointing error in m" axes x1y1

set origin 0,0
set size 0.95,0.33
set key inside top right width -6 maxrows 2
set xrange [0:1]
set yrange [0:1]
unset xlabel
unset x2label
unset ylabel
unset y2label
unset xtics
unset x2tics
unset ytics
unset y2tics
set border 0
set encoding utf8
plot -1,-1 with lines notitle, \
keyentry with lines lc 1 dt 2 title "Analytic", \
keyentry with lines lc 7 title "Cost", \
keyentry with points lt 7 lc 3 title "Subject", \
keyentry with vectors lc 3 lw 2 title "Subject Motion", \
keyentry with lines lc 2 lw 0.5 title "Solver-generated  {{^1}p_{t}}", \
keyentry with lines lc 2 lw 0.5 title "Rotation  {{^1}𝜓̇_{t}}" , \
keyentry with lines lc 8 lw 0.5 title "Solver-generated  {^2}p_{t}}", \
keyentry with lines lc 8 lw 0.5 title "Rotation  {{^2}𝜓̇_{t}}" , \
keyentry with lines lc 6 lw 0.5 title "Solver-generated  {^3}p_{t}}", \
keyentry with lines lc 6 lw 0.5 title "Rotation  {{^3}𝜓̇_{t}}"

#keyentry with lines lc 0 title "Hyperplane", \
#keyentry with dots lc "#FFFFFFFF" title "Trajectory parametes:     ", \
#keyentry with dots lc "#FFFFFFFF" title "v_{min} = ".sprintf("%.1f",vmin)."\tm/s", \
#keyentry with dots lc "#FFFFFFFF" title "v_{max} = ".sprintf("%.1f",vmax)."\tm/s", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜓̇_{min}} = ".sprintf("%.1f",psidot*180.0/pi)."\t°/s", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜓̇_{max}} = ".sprintf("%.1f",14)."\t°/s", \
#keyentry with dots lc "#FFFFFFFF" title "~{{}_F}{{}^T}v = ".sprintf("%.1f",wind)."\tm/s", \
#keyentry with dots lc "#FFFFFFFF" title "{𝛾} = ".sprintf("%.1f",rho*180.0/pi)."\t°", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜙} = ".sprintf("%.1f",-beta*180.0/pi)."\t°", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜑}  is ".phistate, \
#keyentry with dots lc "#FFFFFFFF" title "{𝜃}  is ".thetastate, \
#keyentry with dots lc "#FFFFFFFF" title "{𝛽}  is ".alphastate
#pause -1
