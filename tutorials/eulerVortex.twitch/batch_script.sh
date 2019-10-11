#!/bin/sh
echo "########### Running stativ and moving mesh ##########"
start=`date +%s`
#cd eulerVortex.cyclic.twitch && bash script.sh && cd ..
time1=`date +%s`
#cd eulerVortex.cyclic.twitch.moving && bash script.sh
time2=`date +%s`

echo "########### post-processing #########################"
rm *.eps *.tex
cat eulerVortex.cyclic.twitch/postProcessing/probes/0/U | sed -e "s/(//" -e "s/)//" > tmpU
cat eulerVortex.cyclic.twitch.moving/postProcessing/probes/0/U | sed -e "s/(//" -e "s/)//" > tmpU2

gnuplot -p -e " \
set term epslatex color; \
set out 'compare_p.tex'; \
set multiplot ; set size 1,0.5; \
set format '$%g$' ; \
set tics front ; \
set key bottom right; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 0.4;\
set grid lt 0 lw 2 lc rgb 'black' ; \
set origin 0.0,0.5; \
set xlabel 't\;[s]'; set ylabel '\$p\;[kg.m^{-1}.s^{-2}]\$'; \
p 'eulerVortex.cyclic.twitch/postProcessing/probes/0/p'        u 1:2 with lines lw 18 lc rgb 'orange' t 'static domain', \
  'eulerVortex.cyclic.twitch.moving/postProcessing/probes/0/p' u 1:2 with lines lw 6  lc rgb 'blue'   t 'moving domain' smooth sbezier; \
set origin 0.0,0.0; \
set xlabel 't\;[s]'; set ylabel '\$||U||\$'; \
p 'tmpU'  u 1:(sqrt(\$2*\$2+\$3*\$3+\$4*\$4)) with lines lw 18 lc rgb 'orange' title 'static domain', \
  'tmpU2' u 1:(sqrt(\$2*\$2+\$3*\$3+\$4*\$4)) with lines lw 6  lc rgb 'blue'   title 'moving domain' smooth sbezier; \
unset multiplot; \
set out 'compare_e.tex'; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 100;\
set xlabel 't\;[s]'; set ylabel '\$E\;[m^2.s^{-2}]\$'; \
p 'eulerVortex.cyclic.twitch/postProcessing/probes/0/e'        u 1:2 with lines lw 18 lc rgb 'orange' t 'static domain', \
  'eulerVortex.cyclic.twitch.moving/postProcessing/probes/0/e' u 1:2 with lines lw 6  lc rgb 'blue'   t 'moving domain' smooth sbezier; \
set out 'integral_e.tex'; \
set mxtics 4; set mytics 4; set xtics autofreq 1.; set ytics autofreq 100;\
set format y '\$%.3t.10^{%L}\$'; \
set xlabel 't\;[s]'; set ylabel '\$E\$'; \
set key top right; \
p 'eulerVortex.cyclic.twitch/postProcessing/probes/global/e'        u 1:2 with lines lw 18 lc rgb 'orange' t 'static domain', \
  'eulerVortex.cyclic.twitch.moving/postProcessing/probes/global/e' u 1:2 with lines lw 6  lc rgb 'blue'   t 'moving domain' smooth sbezier; \

set out;"

sed -i -e "s/graphics{/graphics{images\/gnuplot\//" compare_p.tex
sed -i -e "s/graphics{/graphics{images\/gnuplot\//" compare_e.tex
sed -i -e "s/graphics{/graphics{images\/gnuplot\//" integral_e.tex
rm tmpU tmpU2

runtime1=$((time1-start))
runtime2=$((time2-time1))
runtimetot=$((runtime1+runtime2))
echo "time(s): $runtime1 - $runtime2 - $runtimetot"
