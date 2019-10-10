
while (1){
      plot "< cat postProcessing/probes/0/U | sed 's/(//g' | sed 's/)//g'" us 1:(sqrt($2*$2+$3*$3+$4*$4))
      pause 2
      }