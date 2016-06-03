import numpy as np

paths = "matter_power_%s"

pks = ["lin","nl"]

for pk in pks:
    k = np.genfromtxt(paths%pk+"/k_h.txt")
    p = np.genfromtxt(paths%pk+"/p_k.txt")
    z = np.genfromtxt(paths%pk+"/z.txt")
    
    index = np.where(z!=0)[0]
    print p.shape
    print index
    p = p[index][0]
    z = z[index]
    print p.shape, z.shape, z

    pfile = open(paths%pk+"/p_k.txt.new","w")
    pfile.write("# p_k\n")
    zfile = open(paths%pk+"/z.txt.new","w")
    zfile.write("# z\n")
    for i in range(len(p)):
        pfile.write("%.18e "%p[i])
    pfile.close()
    zfile.write("%.18e "%z)
