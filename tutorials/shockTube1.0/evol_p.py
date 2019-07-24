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

print(p)
plt.plot(p, label="p")
plt.legend(loc='best')
plt.savefig("testos.png")
