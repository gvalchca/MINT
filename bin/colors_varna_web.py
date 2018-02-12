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
import time
from wand.image import Image
#import Image, sys

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
    null = " 2>/dev/null "
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

    for i,k in enumerate(keys):
	#with Image(filename=parms["output_pictures"]+f_names[i]+".png") as img:
        #    img.resize(260)
        #    img.save(filename=parms["output_pictures"]+f_names[i]+"_small.png")
	while not os.path.exists(parms["output_pictures"]+f_names[i]+".png"):
	    print "Waiting for PNG nr  "+str(i)
	    time.sleep(1)
	command="convert "+parms["output_pictures"]+f_names[i]+".png -resize 260 "+parms["output_pictures"]+f_names[i]+"_small.png"
        #print command
	os.system(command)

    #changing absolute path to relative:
    parms["output_pictures"] =  parms["output_pictures"][22:]
    print   parms["output_pictures"] 

    javainput = """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"\
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<link rel="icon" type="image/png" href="http://mint.cent.uw.edu.pl/favicon.ico">
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<META NAME="ROBOTS" CONTENT="NOINDEX, NOFOLLOW">
<title>Secondary structure visualization</title>
<style type="text/css">
    body{
        font-family: Helvetica, sans-serif;
        color:#000033;
                float:center;
        }

    .menu{
        text-align:center;
        float:center;
width:800px;
        display: inline;
        font-family:Helvetica,sans-serif;
        height:16px;
        font-size: 12px;
        color:#FFFFFF;
        display:block;
        background-color:#707070;
        }

        .tresc{
                width:800px;
                position: relative;
                top:30px;
                color:#707070 ;
         }

    ul {
        margin: 0;
        padding: 0;
        list-style: none;
        color: #FFFFFF;
        background-color:#707070;
    }
    ul li {
        font-weight:bold;
        position: relative;
        display:inline-block;
        padding-right:25px;
        padding-left:25px;
        }

    li ul {
        position: absolute;
        left: 200px;
        top: 0;
        display: none;
        }

    ul li a {
        display: block;
        text-decoration: none;
        color: #FFFFFF;
        background-color:#707070;
        padding: 5px;
        border-bottom: 0;
        }

    fieldset{
        border:0px;
        padding:1em;
        float:center;
        font-family:Arial;
        }
    .h1{
        color: #606060;
        text-align:center;
        font-family: Arial;
        font-size:x-large;
        }

    .h2{
        color: #606060;
        background-size:750px 200px;
        text-align:right;
        font-family: Courier;
        font-size:50px;
        color:#009966;
        margin-left:10px;
        width:740px;
        height: 205px;
        }

    .h3{
        color: #606060;
        text-align:left;
        font-family: Arial;
        font-size:large;
        margin-left:10px;
        }

    div.p {
        margin-top: 10pt;
        text-align: justify;
        font-family: Arial;
        color: #606060;
        }

    div.a {
        margin-bottom: 12pt;
        text-align: justify;
        font-size:small;
        margin-right:0.5cm;
        }

    div.c {
        text-align: left;
        font-family:Courier New;
        font-size:small;
        }

    div.b {
        text-align: left;
        font-family:Arial;
        font-size:normal;
        font-weight:bold;
        }

    div.d {
        margin-left: 2cm;
        text-align: left;
        text-align: justify;
        font-family: Arial;
        margin-right:0.5cm;
        color: #606060;
        padding: 0.5cm;
       }

    li a:hover{
        background-color:#009966;
        }


</style>
</head>

<body>
<center>

<div class="h2">
<table align="center">
<tr><td>
<a href="http://mint.cent.uw.edu.pl/index.php?strona=MintInt"><img src="http://mint.cent.uw.edu.pl/obrazki/Mint.png"></a>
</td><td><br><br><a href="http://mint.cent.uw.edu.pl/index.php?strona=MintInt" style="text-decoration: none; color:#009966">Bionano MINT server</a></td>
</tr>
</table>
</div>

<div class="menu">
<ul>
<li><a href="http://mint.cent.uw.edu.pl/index.php?strona=MintInt" style="background-color:#009966;">About MINT</a></li><li><a href="http://mint.cent.uw.edu.pl/index.php?strona=MintDocs">Docs</a></li><li><a href="http://mint.cent.uw.edu.pl/index.php?strona=MintRun">Run</a></li><li><a href="http://mint.cent.uw.edu.pl/index.php?strona=MintDownload">Download</a></li><li><a href="http://mint.cent.uw.edu.pl/index.php?strona=MintResults">Results</a></li><li><a href="http://mint.cent.uw.edu.pl/index.php?strona=FAQ">FAQ</a></li>
</ul>
</div>


</div>

<div class="tresc">
<div class="p">

Run title: <b> """+parms["run_title"]+""" </b>
<br>
<br>

Secondary structure visualizations of computed parameters. The pictures were created using <a href="http://varna.lri.fr/">VARNA</a>.
<br>
Lowercase letters (a,u,g,c) indicate modified nucleotides for which MINT has no parameters (see description file). For these nucleotides the values of stacking energy were set to 0 and all recognized hydrogen bonds were classified to non Watson-Crick basepairs. 
<br><br>
For details or tertiary visualizations download the zipped output and refer to the manual.

<br><br>
</div>

<table border="1" align="center"><tr style = "background-color:#606060; color:white; font-family: Courier;  text-align:center;">
<tr> <td align = 'center' colspan='3'> <b> Basepairs parameters </b> </td> </tr>
<tr>
    <td> <a href = """+parms["output_pictures"]+f_names[0]+""".png > <img src =  """+parms["output_pictures"]+f_names[0]+"""_small.png> </a> <br> 
    Hbonds WC: <br>Number of hydrogen bonds in Watson-Crick base pairs. Pseudo-knots are shown as nucleotides with high number of hydrogen bonds but not in the secondary structure.<br>
    </td>
    <td> <a href = """+parms["output_pictures"]+f_names[1]+""".png > <img src =  """+parms["output_pictures"]+f_names[1]+"""_small.png> </a> 
    Hbonds non-WC: <br>Number of hydrogen bonds in non-Watson-Crick base pairs. 
    In this representation bases forming typical Watson-Crick pairs (e.g. in helices) will be visible as having low number of bonds. <br>
    </td>
    <td> <a href = """+parms["output_pictures"]+f_names[2]+""".png > <img src =  """+parms["output_pictures"]+f_names[2]+"""_small.png> </a> 
    <br> Hbonds total: <br>Number of hydrogen bonds in both Watson-Crick and non Watson-Crick basepairs.<br> <br> <br> <br>
    </td>
    </tr>


<tr> <td align = 'center' colspan='3'> <b> Stacking parameters </b> </td> </tr>
<tr>
    <td> <a href = """+parms["output_pictures"]+f_names[3]+""".png > <img src =  """+parms["output_pictures"]+f_names[3]+"""_small.png> </a> 
    Stacking Coulomb: <br> Energy of Coulomb interactions in kcal/mol.</td>
    <td> <a href = """+parms["output_pictures"]+f_names[4]+""".png > <img src =  """+parms["output_pictures"]+f_names[4]+"""_small.png> </a> 
    Stacking VDW: <br> Energy of Van der Waals interactions in kcal/mol.</td>
    <td> <a href = """+parms["output_pictures"]+f_names[5]+""".png > <img src =  """+parms["output_pictures"]+f_names[5]+"""_small.png> </a> 
    Stacking total: <br>  Sum of van der Waals and Coulomb interactions in kcal/mol.
    </td>
</tr>
</table>
<br>
<br>
<br>
<br>
</div>
</center>
</body>
</html>
"""

    f = open(parms["output"], "w")
    print parms["output"] #tu
    f.write(javainput)
    f.close()





