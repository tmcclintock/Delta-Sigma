import numpy as np

text = "# M in Msun/h"

lmlo = np.log10(1e12)
lmhi = np.log10(7e15)
M = np.logspace(lmlo,lmhi,1000,base=10)

outfile = open("mass_list.dat","w")
outfile.write(text)
outfile.write("%e\n"%1e12)
for i in range(len(M)):
    outfile.write("%e\n"%M[i])
outfile.close()
