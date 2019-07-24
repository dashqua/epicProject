import os,sys,subprocess
import matplotlib.pyplot as plt
import numpy as np

plt.yscale("symlog")
plt.figure( figsize=(16, 9) )
cwd = os.getcwd()
command = "grep -rnw "+str(cwd)+" -e 'bounding h' | cut -f2 -d ','"
output = subprocess.check_output(command, shell=True)
f_out = (output.decode("utf-8")).split("\n")
hmin = []; hmax = []; have = [];
for line in f_out:
    if "average"  in line: # ' min: -3448.44 max: 687.601 average: 650.152'
        min = float(line.split("min:")[1].split("max:")[0])
        max = float(line.split("max:")[1].split("average:")[0])
        ave = float(line.split("average:")[1])
        hmin.append(min);hmax.append(max);have.append(ave); 

command = "grep -rnw "+str(cwd)+" -e 'bounding rho' | cut -f2 -d ','"
output = subprocess.check_output(command, shell=True)
f_out = (output.decode("utf-8")).split("\n")
rhomin = []; rhomax = []; rhoave = [];
for line in f_out:
    if "average"  in line: # ' min: -3448.44 max: 687.601 average: 650.152'
        min = float(line.split("min:")[1].split("max:")[0])
        max = float(line.split("max:")[1].split("average:")[0])
        ave = float(line.split("average:")[1])
        rhomin.append(min);rhomax.append(max);rhoave.append(ave); 

        
command = "grep -rnw "+str(cwd)+" -e 'new min' | cut -f4 -d ':'"
output = subprocess.check_output(command, shell=True)
f_out = (output.decode("utf-8")).split("\n")

new_min = []
for x in f_out:
    new_min.append(x)
plt.subplot(131)
plt.plot(new_min, label="new_min")
plt.legend(loc='best')

plt.subplot(132)
plt.plot(hmin, 'r', label="hmin")
plt.plot(hmax, 'b', label="hmax")
plt.plot(have, 'g', label="have")
plt.legend(loc='best')

plt.subplot(133)
plt.plot(rhomin, 'r', label="rhomin")
plt.plot(rhomax, 'b', label="rhomax")
plt.plot(rhoave, 'g', label="rhoave")
plt.legend(loc='best')


plt.savefig("test.png")
