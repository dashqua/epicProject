#!/bin/sh
echo "########### Running stativ and moving mesh ##########"
start=`date +%s`
cd eulerVortex.cyclic.twitch && bash script.sh && cd ..
time1=`date +%s`
cd eulerVortex.cyclic.twitch.moving && bash script.sh
time2=`date +%s`

echo "########### post-processing #########################"
rm *.eps *.tex
cat eulerVortex.cyclic.twitch/postProcessing/probes/0/U | sed -e "s/(//" -e "s/)//" > tmpU
cat eulerVortex.cyclic.twitch.moving/postProcessing/probes/0/U | sed -e "s/(//" -e "s/)//" > tmpU2
paste tmpU tmpU2 > tmpU3

# PLOT PRESSURE
gnuplot -p -e " \
set term epslatex color; \
set tics front ; \
set grid lt 0 lw 2 lc rgb 'black' ; \
\
\
set out 'compare_p.tex'; \
set size 1.0,0.5; \
set key bottom right; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 0.4; \
set xlabel 't\;[s]' offset 0,1,0; set ylabel '\$p\;[kg.m^{-1}.s^{-2}]\$'; \
p 'eulerVortex.cyclic.twitch/postProcessing/probes/0/p'        u 1:2 with lines lw 18 lc rgb 'orange' t 'static domain', \
  'eulerVortex.cyclic.twitch.moving/postProcessing/probes/0/p' u 1:2 with lines lw 6  lc rgb 'blue'   t 'moving domain' smooth sbezier; "

# PLOT VELOCITY
gnuplot -p -e " \
set term epslatex color; \
set tics front ; \
set grid lt 0 lw 2 lc rgb 'black' ; \
\
\
set out 'compare_u.tex'; \
set multiplot layout 2,1; \
set key bottom right ; set size 1.0,0.5; set origin 0.0,0.5 ; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 0.4; \
set xlabel 't\;[s]' offset 0,1,0; set ylabel '\$||U||\$'; \
p 'tmpU'  u 1:(sqrt(\$2*\$2+\$3*\$3+\$4*\$4)) with lines lw 18 lc rgb 'orange' title 'static domain', \
  'tmpU2' u 1:(sqrt(\$2*\$2+\$3*\$3+\$4*\$4)) with lines lw 6  lc rgb 'blue'   title 'moving domain' smooth sbezier; \
set key top right ;  set size 1.0,0.5;  set origin 0.0,0.0 ; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 0.1; \
set xlabel 't\;[s]' offset 0,1,0; set ylabel '\$L^2\$-error'; \
p 'tmpU3'  u 1:(sqrt((\$2-\$6)*(\$2-\$6))) with lines lw 6 lc rgb 'red' title '\$|| u-u_{ALE}||\$' smooth sbezier,   \
  'tmpU3'  u 1:(sqrt((\$3-\$7)*(\$3-\$7))) with lines lw 6 lc rgb 'blue' title '\$|| v-v_{ALE}||\$' smooth sbezier,  \
  'tmpU3'  u 1:(sqrt((\$4-\$8)*(\$4-\$8))) with lines lw 6 lc rgb 'green' title '\$|| w-w_{ALE}||\$' smooth sbezier; \
unset multiplot;"

#PLOT ENERGY
gnuplot -p -e " \
set term epslatex color; \
set tics front ; \
set grid lt 0 lw 2 lc rgb 'black' ; \
\
\
set out 'compare_e.tex'; \
set multiplot; set size 1.0,0.5; \
set key bottom right; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 100;\
set xlabel 't\;[s]' offset 0,1,0; set ylabel '\$E\;[m^2.s^{-2}]\$'; \
set origin 0.0,0.5;
p 'eulerVortex.cyclic.twitch/postProcessing/probes/0/e'        u 1:2 with lines lw 18 lc rgb 'orange' t 'static domain', \
  'eulerVortex.cyclic.twitch.moving/postProcessing/probes/0/e' u 1:2 with lines lw 6  lc rgb 'blue'   t 'moving domain' smooth sbezier; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 100;\
set format y '\$%.3t.10^{%L}\$'; \
set xlabel 't\;[s]' offset 0,1,0; set ylabel '\$E\$'; \
set key top right; \
set origin 0.0,0.0;
p 'eulerVortex.cyclic.twitch/postProcessing/probes/global/e'        u 1:2 with lines lw 18 lc rgb 'orange' t 'static domain', \
  'eulerVortex.cyclic.twitch.moving/postProcessing/probes/global/e' u 1:2 with lines lw 6  lc rgb 'blue'   t 'moving domain' smooth sbezier; \
unset multiplot;\
set out;"
#set out 'integral_e.tex'; \




sed -i -e "s/graphics{/graphics{images\/gnuplot\//" compare_p.tex
sed -i -e "s/graphics{/graphics{images\/gnuplot\//" compare_u.tex
sed -i -e "s/graphics{/graphics{images\/gnuplot\//" compare_e.tex
rm tmpU tmpU2 tmpU3


# conversion to pdf






runtime1=$((time1-start))
runtime2=$((time2-time1))
runtimetot=$((runtime1+runtime2))
echo "time(s): $runtime1 - $runtime2 - $runtimetot"
