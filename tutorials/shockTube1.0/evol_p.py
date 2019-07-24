import sys,os
import numpy as np
import matplotlib.pyplot as plt

cwd = os.getcwd()
lstdir = os.listdir(cwd)
x=[]
p = []
for edir in lstdir:
    if os.path.isdir(cwd+"/"+edir):
        try:
            with open(cwd+"/"+edir+"/p", "r") as f:
                p.append(float(f.readlines()[200]))
                x.append(float(edir))
        except:
            continue

with open("evolp_"+sys.argv[1], 'w') as f:
    for i in range(len(x)):
        f.write(str(x[i])+"\t"+str(p[i])+"\n")
#print(p)
#plt.plot(p, label="p")
#plt.legend(loc='best')
#plt.savefig("testos.png")
