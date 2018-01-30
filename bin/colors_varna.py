#/****************************************************************************
#*                                                                          *
#*   Copyright (C) 2013: University of Warsaw                               *
#*   Author: Anna Gorska <agorska@cent.uw.edu.pl>                           *
#*   Author: Maciej Jasinski <maciejj@cent.uw.edu.pl>                       *
#*                                                                          *
#*   This program is free software; you can redistribute it and/or modify   *
#*   it under the terms of the GNU General Public License version 2,        *
#*   as published by the Free Software Foundation.                          *
#*                                                                          *
#*   This program is distributed in the hope that it will be useful,        *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
#*   GNU General Public License for more details.                           *
#*                                                                          *
#*   You should have received a copy of the GNU General Public License      *
#*   along with this program; if not, write to the                          *
#*   Free Software Foundation, Inc.,                                        *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              *
#*                                                                          *
#****************************************************************************/

# !/usr/bin/env python
from string import *
import os

# Use: colors_varna.py nucleotides_evaluate.csv
# javaindex.html sequence dotbracket


def avg(max1, min2):
#    return max1-((abs(max1)+abs(min2))/2)
    return max1- ((abs(max1 -min2))/2)


def run(parms):
    if parms["mode"] == "Traj":
        indeces = ["num", "WC", "WC-hbonds sd", "D", "non-WC-hbonds sd",
                   "T", "WC-total sd", "C", "Coulomb sd", "V", "VDW sd",
                   "S", "sum sd"]
    else:
        indeces = ["num", "WC", "D", "T", "C", "V", "S"]

    ll = []
    for k in parms["nucleotides_eval"]:
        tmp = [int(k.split(":")[1])]
        tmp.extend(parms["nucleotides_eval"][k])
        ll.append(tmp)
    ll.sort()

    dd = {"WC": [], "D": [], "T": [], "C": [], "V": [], "S": []}
    keys = ["WC", "D", "T", "C", "V", "S"]

    for i, k in enumerate(keys):
        for j in ll:
            dd[k].append(j[indeces.index(k)])

    c1 = "mn:#00FFFF,av:#FFFF00,mx:#FF0000"
    c2 = "mn:#FF0000,av:#FFFF00,mx:#00FFFF"
    cc = []

    for k in ["WC", "D", "T"]:
        cc.append(c1.replace("mn", str(min(dd[k]))).
                  replace("mx", str(max(dd[k]))).
                  replace("av", str(avg(max(dd[k]), min(dd[k])))))

    for k in ["S", "C", "V"]:
        cc.append(c2.replace("mn", str(min(dd[k]))).
                  replace("mx", str(max(dd[k]))).
                  replace("av", str(avg(max(dd[k]), min(dd[k])))))

    dd_str = {}
    for k in keys:
        dd_str[k] = ",".join([str(i) for i in dd[k]])

    sequence = parms["sequence"]
    dotbracket = parms["dot_bracket"]
    null = " 2>/dev/null &"
    captions = ["Hbonds WC", "Hbonds non-WC", "Hbonds total",
                "Stacking Coulomb", "Stacking VDW", "Stacking total"]
    f_names = ["Hbonds-WC", "Hbonds-non-WC", "Hbonds-total",
                "Stacking-Coulomb", "Stacking-VDW", "Stacking-total"]
    java_start = "java -cp "+parms["home"] + \
        "/bin/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd"\
        " -sequenceDBN \""+sequence+"\" -structureDBN \"" +\
        dotbracket + "\" -resolution \"4.0\" "

    for i, k in enumerate(keys):
        part = "-colorMap \"" +\
            dd_str[k] + "\" -colorMapCaption \"" + captions[i] +\
            "\" -colorMapStyle \"" + cc[i]+"\" -o \"" + \
            parms["output_pictures"] + f_names[i]+".png\""
        os.system(java_start+part+null)
