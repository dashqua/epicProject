import os, sys

nbelem=206
with open("internalp", "w") as f:
    f.write("internalField nonuniform List<scalar>"+"\n")
    f.write(str(nbelem)+"\n")
    f.write("("+"\n")

    for i in range(int(nbelem/2)):
        f.write(str(3)+"\n")
        f.write(str(1)+"\n")

    f.write(")"+"\n")
    f.write(";"+"\n")
