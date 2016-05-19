"""
This is a script to run the weak lensing suite of code.
"""
import numpy as np
import os

#This is the base path to the data
cosmosis_path = "/home/tmcclintock/Cosmosis/"
ds_path = cosmosis_path+"my_tests/my_output/"

def make_start_cosmology(lam,i,j):
    lM = np.log10(1.56e14*(lam/30.0)**1.31) #log10 M200_c, but w/e
    Rmis = np.exp(-1.13)*R_lams[i,j]
    print "\tMCMC cosmology start mass: ",lM
    inpath  = cosmosis_path+"my_tests/WL_ANALYSIS/cosmology.ini"
    outpath = cosmosis_path+"my_tests/WL_ANALYSIS/temp.ini"
    in_file = open(inpath,"r")
    out_file = open(outpath,"w")
    for line in in_file:
        parts,outline = line.split(), line
        if len(parts)>0:
            if parts[0]=="log10_cluster_mass":
                outline = parts[0]+" "+parts[1]+" "+parts[2]+" "+str(lM)+" "+parts[4]+"\n"
            if parts[0]=="sigma_miscentering":
                outline = "%s = %f\n"%(parts[0],Rmis)
            out_file.write(outline)
    in_file.close(),out_file.close()
    os.system("mv "+outpath+" "+inpath)
    return

z = 0.5
cosmoname = "assumed"
command = "z=%f cosmoname=%s cosmosis my_tests/wl_test/pipeline.ini"%(z,cosmoname)

os.system(command)

