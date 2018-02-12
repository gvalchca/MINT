
#****************************************************************************
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

# !/usr/bin/python
import sys
import os

floats = ["time_cutoff", "max_memory_GB", "cutoff_stacking",
          "cutoff", "vdw_cutoff_stacking", "OP_stacking_distance_cutoff",
          "h_bond_l", "h_bond_angle", "margin", "vdw_cutoff_stacking",
          "pairing_cutoff"]
ints = ["first_frame", "last_frame", "vmd", "only_analysis",
        "threads", "stride", "create_csvs", "MerLength", "all_nucs_DNA"]
strings = ["out_name", "file_name", "file_dcd", "h_bond_atom",
           "table_nucleotides", "table_charges", "force_field",
           "list_of_modified_nucs", "Mode", "working_dir",
           "energy_library_directory", "nucleotides"]
lists = ["chains_names", "pdb_list", "files_dcd"]

configuratations = {"obligatory": ["Mode", "table_nucleotides",
                                   "table_charges", "list_of_modified_nucs",
                                   "force_field", "margin", "threads", "cutoff",
                                   "cutoff_stacking", "vdw_cutoff_stacking",
                                   "OP_stacking_distance_cutoff", "h_bond_atom",
                                   "h_bond_l", "h_bond_angle"],
                    "Single": ["file_name", "out_name"],
                    "Traj": ["file_name", "file_dcd", "first_frame",
                             "last_frame", "time_cutoff", "stride", "out_name",
                             "pairing_cutoff"],
                    "Download": ["pdb_list"],
                    "NNM": ["MerLength", "pdb_list", "energy_library_directory"],
                    "TrajHbonds": ["files_dcd", "working_dir", "nucleotides"]}
modes = ["Traj", "Single", "Download", "NNM", "TrajHbonds"]


def check(parms, defined_parameters):
    if "Mode" not in parms.keys() or parms["Mode"] not in modes:
        print "Please specify proper Mode parameter, ",
        print "you can run ", " ".join(modes), " modes."
        sys.exit()

    missing = [i for i in configuratations["obligatory"]
               if i not in defined_parameters]
    if len(missing) > 1:
        print "You miss parameters in your CONF file: ", "\n".join(missing)
        sys.exit()

    missing_spec = set(configuratations[parms["Mode"]]).difference(set(defined_parameters))
    if len(list(missing_spec)) >= 1:
        print "While running MINT in ", parms["Mode"], "Mode,",
        print "you need to specify aditionally parameters: ", " ".join(missing_spec)
        sys.exit()

    all_there_can_be = []
    for i in [floats, ints, strings, lists]:
        all_there_can_be.extend(i)

    aditional = [i for i in defined_parameters if i not in all_there_can_be]
    if len(aditional) > 0:
        print "I do not understand the parameters: "+" ".join(aditional)
        sys.exit()

    if parms["Mode"] == "Traj" and parms["stride"] == 0:
        print "Cant run with the 0 stride, running with stride=1 "
        parms["stride"] = 1

    if "create_csvs" not in parms.keys() or parms["create_csvs"] not in [0, 1]:
        print "You did not specify create_csvs parameter, only xls is going " +\
            "to be created."
        parms["create_csvs"] = 0

    if "all_nucs_DNA" not in parms.keys() or parms["all_nucs_DNA"] not in [0, 1]:
        parms["all_nucs_DNA"] = 0
    if  parms["all_nucs_DNA"] == 1:
	print "all_nucs_DNA set to 1. Nucleotides with standard names: A,C,G,RA,RC,RG,ADE,CYT,GUA, etc. wil be treated as DNA not RNA as usually."


    if (parms["Mode"] == "NNM" and
       (parms["MerLength"] == 0 or parms["MerLength"] % 2 != 0)):
        print "MerLength has to be divisible by 2, and grater then 0"
        print "Setting to 4.."
        parms["MerLength"] = 4

    if parms["Mode"] in ["Single", "Traj"]:
        if "/" not in parms["out_name"]:
            parms["out_name"] = "/".join(parms["file_name"].split("/")[:-1]) + \
                                "/" + parms["out_name"]
        else:
            catalogues = parms["out_name"].split("/")[:-1]
            directory = "/".join(catalogues[:-1])
            i = 1
            print 'Checking if directory:', directory, 'exists'
            while True:
                if os.path.exists(directory):
                    pass
                else:
                    print directory, "doesnt exist, creating "
                    os.mkdir(directory)
                try:
                    directory += "/" + catalogues[i]
                    i += 1
                except Exception:
                    break
    return True


def read_pdbs_list(PARMS):
    if len(PARMS["pdb_list"]) == 1 and PARMS["pdb_list"][0] in os.listdir("./"):
        f = open(PARMS["pdb_list"][0], "r")
        delim = ""
        out = []
        ll = f.readlines()
        for i in ll:
            delim = [de for de in [",", ";"] if de in i][0]
            out.extend([j.strip() for j in i.split(delim)])
        PARMS["pdb_list"] = out


def read_in_parameters_from_file(file_in):
    f = open(file_in, "r")
    p = [i.split("#")[0].replace("\n", "").replace(",", " ")
         for i in f.readlines() if not i.startswith("#")]
    f.close()

    params = {"only_analysis": 0}
    defined_parameters = []

    for i in p:
        tmp = i.split(":")
        if tmp[0].strip():
            tmp[0] = tmp[0].strip()
            try:
                if tmp[0] in floats:
                    params[tmp[0]] = float(tmp[1])
                elif tmp[0] in ints:
                    params[tmp[0]] = int(tmp[1])
                elif tmp[0] in lists:
                    params[tmp[0]] = tmp[1].split()
                elif tmp[0] in strings:
                    params[tmp[0]] = tmp[1].strip() #  tmp[1].replace(" ", "")
                defined_parameters.append(tmp[0])
            except Exception:
                params[tmp[0]] = None
    if "pdb_list" in params.keys():
        read_pdbs_list(params)
    if check(params, defined_parameters):
        return params
