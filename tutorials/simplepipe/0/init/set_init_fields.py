import sys
import numpy as np
nbe = int(sys.argv[1])
nbeq =  int(nbe/4)
nbeqq = int(nbe/nbeq)
nbeqq2 = int(nbeqq/2)
print(nbe)
print(nbeq)
print(nbeqq)
print(nbeqq2)
print(np.sqrt(nbe))

strval=""
for z in range(nbeq):
    for i in range(nbeqq2):
        strval+="350"
        for j in range(nbeqq2 -1):
            strval+="\n273"

content ="// internal Field uniform 273.15;  \n\
internalField nonuniform List<scalar> \n\
%d                                    \n\
(                                     \n\
%s                                    \n\
)                                     \n\
;" % (nbe,strval)
f = open(sys.argv[0] + "internalT", "w")
f.write(content)
f.close()


#for j in range(128):
#    print("500000")
#    for i in range(95):
#        print("100000")
        
    #303975
#101325



strval=""
for z in range(nbeq):
    for i in range(nbeqq2):
        strval+="1"
        for j in range(nbeqq2 -1):
            strval+="3"

content ="// internal Field uniform 273.15;  \n\
internalField nonuniform List<scalar> \n\
%d                                    \n\
(                                     \n\
%s                                    \n\
)                                     \n\
;" % (nbe,strval)
f = open(sys.argv[0] + "internalp", "w")
f.write(content)
f.close()
