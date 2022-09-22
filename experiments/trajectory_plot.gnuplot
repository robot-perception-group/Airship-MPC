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
set arrow 1 from 0,0,0 to 0,-6,0 lc 3 lw 2

set parametric
set urange [0:2*pi]
set vrange [0:100]
set xrange [24:-24]
set yrange [24:-24]
set zrange [:5] reverse
set view equal xyz


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
#set xlabel "Y (m)" rotate parallel offset graph 0.0,0.1,0
#set xlabel "Y (m)" rotate parallel offset 3,2.0,0
#set xlabel "Y (m)" rotate parallel offset -9,1.0,0
#set ylabel "X (m)" rotate parallel offset graph -0.1,-0.1,0
#set ylabel "X (m)" rotate parallel offset 3,2.0,0 
#set zlabel "-Z (m)" offset graph 0.3,0,0.8
#set zlabel "Z (m)" offset 5,5,0
#splot "trajectory.csv" using ($2):($1):($3) with lines  lc 2 lw 2 title "MPC trajectory", \
#"" using 8:7:($9) lc 7 with lines title "sensor pointing error in m", \
#f3x(u,v,0),-f3y(u,v,0),-f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", \
#f4x(u,v,wind),-f4y(u,v,wind),-f4z(u,v,wind) with lines lc 1 title "valid trajectories moving subject", \
#persx(u,v),-persy(u,v),-persz(u,v) with lines lc 1 title "subject", \
#"trajectory.csv" using ($2):($1):($3) with lines  lc 2 lw 0.5 title "MPC trajectory", \
#"" using 8:7:($9) lc 7 with lines title "sensor pointing error in m"
unset ztics
unset zlabel
set size 0.36,0.67
set origin 0.32,0.30
set view 359.999,0,2.5
set xlabel "Y (m)" rotate parallel offset graph 0,-2.16
set ylabel "X (m)" rotate parallel offset graph 0.95,0
set xtics offset graph 0,-1.24
set ytics offset graph 0.08,0
splot f3x(u,v,0),-f3y(u,v,0),-f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", \
f4x(u,v,wind),-f4y(u,v,wind),-f4z(u,v,wind) with lines lc 1 title "valid trajectories moving subject", \
persx(u,v),-persy(u,v),-persz(u,v) with lines lc 3 title "subject", \
"trajectory.csv" using ($2):($1):($3) with lines  lc 2 lw 2 title "MPC trajectory"

#, \
#"" using 8:7:($9) lc 7 with lines title "sensor pointing error in m"
#set ztics in
#unset ytics
#unset ylabel
#set xlabel "Y (m)" norotate offset graph 0.68,0,1.0
#set xtics offset graph 0.04,0,.28
#set zlabel "Z (m)" rotate parallel offset graph 0.15,0
#set ztics offset graph 0.05,0
#set size 0.33,0.8
#set origin 0.66,0.45
#set view 90,0,1.9
#splot f3x(u,v,0),-f3y(u,v,0),-f3z(u,v,0) with lines lc 0 title "hyperplane stationary subject", \
#f4x(u,v,wind),-f4y(u,v,wind),-f4z(u,v,wind) with lines lc 1 title "valid trajectories moving subject", \
#persx(u,v),-persy(u,v),-persz(u,v) with lines lc 1 title "subject", \
#"trajectory.csv" using ($2):($1):($3) with lines  lc 2 lw 2 title "MPC trajectory", \
#"" using 8:7:($9) lc 7 with lines title "sensor pointing error in m"
set ytics nomirror

set ylabel rotate parallel offset graph 0.08,0
set ytics offset graph 0.02,0
set xtics offset graph 0,0.1
set size 0.36,0.32
set origin -0.05,0.68
set autoscale x
set autoscale y
set yrange [-22.5:22.5]
set samples 500
set trange [0:500]
set format x ''
set grid xtics
unset arrow 1
set ylabel "X (m)"
unset xlabel
unset lmargin
plot "trajectory.csv" using ($0*0.25):($1) with lines  lc 2 lw 2 title "MPC trajectory X", t,-f3y(-pi/2+t*psidot,reference_radius,wind)  with lines lc 1 title "valid trajectory moving subject X" axis x1y1

#set xlabel "Time (s)" offset 0,0.5
set x2label "Time" offset graph 0,-0.3
set size 0.36,0.32
set origin -0.05,0.345
set ylabel "Y (m)"
plot "trajectory.csv" using ($0*0.25):($2) with lines  lc 2 lw 2 title "MPC trajectory Y", t,f3x(-pi/2+t*psidot,reference_radius,wind)  with lines lc 1 title "valid trajectory moving subject Y" axes x1y1
#set grid vertical xtics y2tics
set xtics offset graph 0.04,0.3
set size 0.36,0.32
#set origin -0.05,0.05
set origin 0.60,0.68
set x2label "Time" offset graph 0,-0.3
#set format x '% h'
set ylabel "Z (m)"
set yrange [0:-45]
set y2tics
#set format y2 '% h'
set y2range [0:8.5]
set y2label "Cost - AoA 𝛽 (°)" offset graph -0.05,0 
#set key at screen 0.98,0.5
plot "trajectory.csv" using ($0*0.25):($3) axes x1y1 with lines  lc 2 lw 2 title "MPC trajectory Z", \
t,-f3z(-pi/2+t*psidot,reference_radius,wind)  with lines lc 1 title "valid trajectory moving subject Z" axes x1y1, \
"trajectory.csv" using ($0*0.25):($12) with lines lc 4 title "AoA Alpha in °" axes x1y2, "" using ($0*0.25):(sqrt($7**2+$8**2+$9**2)) with lines lc 7 title "sensor pointing error in m" axes x1y2
#set object 5 rect from -1000,-1000 to 1000,1000
#set object 5 rect behind fc rgb 'white' fillstyle solid 0.0 noborder
set size 0.35,0.31
set origin 0.30,0
#set lmargin 0
#unset rmargin
#set grid x2tics
#set x2tics offset graph 0.04,-0.3
#set x2label "Time (s)" offset graph 0,-0.5
#set x2tics mirror
#unset xtics
#set format x2 '% h'
#set format x ''
#set yrange [0:8.5]
#set ylabel "Sensor error (m),\nAoA 𝛽 (°)"
#set ylabel  offset graph 0.03,0
#plot "trajectory.csv" using ($0*0.25):($12) with lines lc 4 title "AoA Alpha in °" axes x1y1, "" using ($0*0.25):(sqrt($7**2+$8**2+$9**2)) with lines lc 7 title "sensor pointing error in m" axes x1y1
#unset object 5

set origin 0.63,0.33
set size 0.32,0.33
set key inside top center width -8
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
keyentry with lines lc 1 title "Analytic", \
keyentry with points lt 7 lc 3 title "Subject", \
keyentry with vectors lc 3 lw 2 title "Subject Motion", \
keyentry with lines lc 2 lw 2 title "Solver-generated", \
keyentry with lines lc 4 title "Angle of Attack {𝛽} (°)", \
keyentry with lines lc 7 title "Cost"


#keyentry with lines lc 0 title "Hyperplane", \
#keyentry with dots lc "#FFFFFFFF" title "   ", \
#keyentry with dots lc "#FFFFFFFF" title "Trajectory parametes:     ", \
#keyentry with dots lc "#FFFFFFFF" title "v_{min} = ".sprintf("%.1f",vmin)."\tm/s", \
#keyentry with dots lc "#FFFFFFFF" title "v_{max} = ".sprintf("%.1f",vmax)."\tm/s", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜓̇} = ".sprintf("%.1f",psidot*180.0/pi)."\t°/s", \
#keyentry with dots lc "#FFFFFFFF" title "~{{}_F}{{}^T}v = ".sprintf("%.1f",wind)."\tm/s", \
#keyentry with dots lc "#FFFFFFFF" title "{𝛾} = ".sprintf("%.1f",rho*180.0/pi)."\t°", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜙} = ".sprintf("%.1f",-beta*180.0/pi)."\t°", \
#keyentry with dots lc "#FFFFFFFF" title "{𝜑}  is ".phistate, \
#keyentry with dots lc "#FFFFFFFF" title "{𝜃}  is ".thetastate, \
#keyentry with dots lc "#FFFFFFFF" title "{𝛽}  is ".alphastate
#pause -1

