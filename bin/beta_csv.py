#!/usr/bin/env python
from string import *
import sys

# Use: beta_csv.py file.in file.pdb new.pdb
# occup - WC-WC; beat = others

f = open(sys.argv[1], "r")
d = {}
for i in f.readlines():
    tmp = i.split(",")
    d[int(tmp[1])] = tmp[-2:]
f.close()

f1 = open(sys.argv[2], "r")
f2 = open(sys.argv[3], "w")

nuc = ["RA", "RU", "RG", "RC", "A", "U",
       "G", "C", "AD", "UR", "GU", "CY"]

for linia in f1.readlines():
    if linia[0:4] == ("ATOM") and linia[17:20].replace(" ", "")[:2] in nuc:
        resid = int(linia[23:26])
        occup = '{0:.2f}'.format(round(float(d[resid][0]), 2))
        beta = '{0:.2f}'.format(round(float(d[resid][1]), 2))
        f2.write(linia[:54] + " "*(6-len(occup)) + occup +
                 " "*(6-len(beta)) + beta + linia[66:])
f2.close()
f1.close()
