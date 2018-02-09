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

#3

# !/usr/bin/python
import gzip
import math
import Bio
import sys
import os
from collections import Counter
import multiprocessing
import MDAnalysis
import pickle
import numpy as np

from operator import itemgetter
#from pympler.asizeof import asizeof
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/bin")

import rna_structure as RS
import RNAStructML as XML
import motif_cluster as MC
import csvToxls as cX
import colors_varna as Varna
import download_and_protonate as Download
import param_reader as ParamReader
import time

from Bio.PDB.MMCIFParser import MMCIFParser

if len(sys.argv) < 2:
    print "Specify CONFIG file!"
    sys.exit()

elif len(sys.argv) > 2:
    print "Too many arguments! Specify single configuration file!"
    sys.exit()


def timer(f):
    def g(*args, **kwargs):
        start = time.clock()
        result = f(*args, **kwargs)
        open("../time.txt", 'a').write(f.__name__ + ": " +
                                       str(time.clock() - start) + "\n")
        return result
    return

# warnings.simplefilter("error", RuntimeWarning)
aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "HSD",
      "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
      "VAL"]
NUCLEOTIDES = ["A", "U", "C", "G", "ADE", "URA", "CYT", "GUA", "RA", "RG",
               "RC", "RU", "T", "DA", "DC", "DT", "DG", "G5", "G3", "A5",
               "A3", "C3", "C5", "U5", "U3", "T5", "T3", "DG5", "DG3", "DA5",
               "DA3", "DC3", "DC5", "DU5", "DU3", "DT5", "DT3", "LG5", "LG3",
               "LA5", "LA3", "LC3", "LC5", "LU5", "LU3", "LT5", "LT3", "LA",
               "LC", "LT", "LG", "LU", "TPN", "GPN", "CPN", "APN"]

a = [NUCLEOTIDES.extend([i+"5", i+"3"]) for i in ["RA", "RG", "RC", "RU"]]

ions = ["MG", "ZN", "NA", "K", "HOH", "SOD", "CLA", "TIP", "MN", "SO4", "WAT"]
possible_ps = ["O1P", "O2P", "OP1", "OP2"]
d = ["WC donor", "Hoogsteen donor", "Sugar donor"]
a = ["WC acceptor", "Hoogsteen acceptor", "Sugar acceptor"]
IF_MODIFIED = {}
UNKNOWN_NUC = []
HAVE_FORCE_FIELD = []
masses = {'C': 12.0, 'N': 14.0, 'O': 16.0, 'P': 15.0, 'H': 1.0}


# read csv file and create the table
def read_from_csv(filename):
    table = []
    f = open(filename, "r")
    for i in f.readlines():
        table.append(i.replace("\n", "").split(";"))
    f.close()
    return table


# function uncompressing the file
def uncompress(filename):
    f_in = gzip.open(filename, 'rb')
    pdb = filename.split(".")[0].replace("pdb", "") + ".ent"
    f_out = open(pdb, "w")
    f_out.write(f_in.read())
    f_in.close()
    f_out.close()
    print "File with the strucuture: " + pdb
    return pdb


def if_modified_list(nucs, PARMS, if_mod=False):
    if if_mod is False:
        if_mod = IF_MODIFIED
    for i in nucs:
        if_mod[i.get_id()] = if_modified(i, PARMS)


def if_modified(res, PARMS):
    if not isinstance(res, str):
        resname = res.get_resname()
    else:
        resname = res
    resname = resname.strip()
    if resname in NUCLEOTIDES:
        return False
    origin = check_origin(resname, PARMS)
    if origin:
        return origin
    else:
        return "Completly unknown Nucleotide ", resname


# #@timer
def get_nucleic_from_pdb(PARMS):
    # Using BioPython to readin the PDB file
    filename = PARMS["file_name"]
    if filename.endswith(".cif"):
        parser = MMCIFParser()
        structure = parser.get_structure("1", filename)
    else:
        p = Bio.PDB.PDBParser(PERMISSIVE=1)
        structure = p.get_structure("1", filename)
    nucleotides = []
    seq = ""
    # searching for RNA
    chains_names = PARMS["chains_names"]
    if not chains_names:
        for i in structure[0]:
            chains_names.append(i.get_id())
    for chain_name in chains_names:
        for residue in structure[0][chain_name]:
                name = residue.get_resname().replace(" ", "")
                if name in NUCLEOTIDES:
                    nucleotides.append(residue)
                    seq += name + " "
                elif name in UNKNOWN_NUC:
                    nucleotides.append(residue)
                    seq += name + " "
                elif name not in aa and name not in ions:
                    w = "\nWARNING RESIDUE " + str(name) +\
                        "\ndoesn't seem to be nucleotide, amioacid or an ion."+\
                        "\nRunning the analysis without it."+\
                        "\nIf you wish for it to be analysed add parameters.\n"
                    print w
                    PARMS["OUT_FILE"].write(w)

    tt = ("\n  " + str(len(nucleotides)) +
          " nucleotides, sequence: \n" + seq + "\n")
    PARMS["OUT_FILE"].write(tt)
    c = Counter(seq.split(" ")[:-1])
    print "\nNucleotides statistics:"
    PARMS["OUT_FILE"].write("Nucleotides statistics:\n")
    k = c.keys()
    if_modified_list(nucleotides, PARMS)
    for i in k:
        if i in NUCLEOTIDES:
            tmp = str(i) + " "*(3-len(i)) + "-> " +\
                str(c[i]) + " known parameters"
            print tmp
            PARMS["OUT_FILE"].write(tmp+"\n")
        elif i in UNKNOWN_NUC:
            tmp = str(i) + " "*(3-len(i)) +\
                "-> " + str(c[i]) + " Unknown parameters"
            PARMS["OUT_FILE"].write(tmp+"\n")
            print tmp
    return nucleotides


# read pdb, iter through its nucleotides:
# @timer
def read_in_charges(nucleotides, PARMS):
    f = open(PARMS["table_charges"])
    tmp = f.readlines()
    f.close()
    l = [i.replace("\n", "").split(";") for i in tmp[1:]]
    column = 0
    if PARMS["force_field"] == "AMBER":
        column = [2, 5, 6]
    elif PARMS["force_field"] == "CHARMM":
        column = [3, 7, 8]
    elif PARMS["force_field"] == "MY_OWN":
        column = [4, 9, 10]
    d = {}
    nuc = list(set([names(n.get_resname()) for n in nucleotides]))
    for i in l:
        # charge, vdw_radius, vdw_depth
        try:
            if i[0] in nuc:   # and i[0] not in UNKNOWN_NUC:
                d[i[0]+"_"+i[1]] = [float(i[column[0]]), float(i[column[1]]),
                                    float(i[column[2]])]
        except Exception:
            
            print "We do not have data for ", i[0], " in ", i, " in"
            print PARMS["force_field"], " force field."
            suggested = [i for i in ["AMBER", "CHARMM"]
                         if i not in [PARMS["force_field"]]]
            print "Try running with ", " ".join(suggested),
            print " or defining your own field."
            sys.exit()
    return d


def coords_and_charges_and_vdw(nuc, charges, TimeTable, N):
    out = {}
    name = names(nuc.get_resname())
    atoms_names = [i.split("_")[1] for i in charges.keys()
                   if (i.split("_")[0] == name)]
    for at in atoms_names:
        try:
            out[at] = [charges[name+"_"+at],
                       TimeTable[nuc[at].get_full_id()][N], at]
        except Exception:
            try:
                tmp = at[-1] + at[0:-1]
                out[at] = [charges[name + "_" + at],
                           TimeTable[nuc[tmp].get_full_id()][N], at]
            except Exception:
                # print name+"_"+at, " ",charges[name+"_"+at]
                pass
    return out


def vdw_and_coulomb(base1, base2):
    energy_c = 0.0
    energy_vdw = 0.0

    for i in base1:
        for t in base2:
            atom_distance = np.linalg.norm(base1[i][1] - base2[t][1])
            # kcal/mol
            energy_c += 332*base1[i][0][0]*base2[t][0][0]/atom_distance
            #  If you measure charge in units of the electron charge,
            #  and distance in Angstroms,
            #  then an electrostatic energy looks like:
            #  E (kcal/mol) =  332 * q1*q2/r
            vdw_r = (base1[i][0][1] + base2[t][0][1])
            energy_vdw += (np.sqrt(base1[i][0][2]*base2[t][0][2]) *
                           (((vdw_r/atom_distance)**12) -
                            2*((vdw_r/atom_distance)**6)))  # kcal/mol
    return (round(energy_c, 2),
            round(energy_vdw, 2),
            round(energy_c + energy_vdw, 2))


# measure distances
def pi_stacking(n1, n2, base1, base2, PARMS):
    ps = {}
    t1 = False
    t2 = False
    for t in possible_ps:
        # bez wodorow
        try:
            t1 = atom_to_base_mas_center_distance(base2, base1[t])
            t2 = atom_to_base_mas_center_distance(base1, base2[t])
            if t1 and t1 < PARMS["OP_stacking_distance_cutoff"]:
                ps[RS.resid_description(n2) + "_" +
                   RS.resid_description(n1) + ":" +
                   t] = [t1, vdw_and_coulomb(only_base(base2),
                         phosphate_group(base1))]

            if t2 and t2 < PARMS["OP_stacking_distance_cutoff"]:
                ps[RS.resid_description(n1) + "_" +
                   RS.resid_description(n2) + ":" +
                   t] = [t2, vdw_and_coulomb(only_base(base1),
                         phosphate_group(base2))]
        except Exception:
            pass
    return ps


def phosphate_group(base):
    phosphate = {}
    ps = [i for i in base.keys() if "P" in i]
    for i in ps:
        phosphate[i] = base[i]
    return phosphate


def only_base(base):
    not_phosphate = {}
    ps = [i for i in base.keys() if "P" not in i]
    for i in ps:
        not_phosphate[i] = base[i]
    return not_phosphate


def atom_to_base_mas_center_distance(base, op):
    mass_list = []
    dist_list = []
    for i in base:
        if i not in possible_ps:
            m = masses[i[0]]
            dist_list.append(m*dist(base[i][1], op[1]))
            mass_list.append(m)
    return round(sum(dist_list)/float(sum(mass_list)), 2)


def atom_to_base_distance(base, op):
    dist_list = []
    for i in base:
        if i not in possible_ps:
            dist_list.append(dist(base[i][1], op[1]))
    return round(sum(dist_list)/float(len(dist_list)), 2)


def vmd_for_dipole_selection(vdw_and_charges):
    t = ["A", "C", "G", "U"]
    o = {}
    for m in t:
        o[m] = [i.split("_")[1] for i in vdw_and_charges
                if i.split("_")[0] == m]
    return o


def measure_for_all(nucleotides, charges, PARMS, TimeTable, N):
    table = PARMS["table_nuc_read"]
    margin = int(PARMS["margin"])
    h_bond_l = float(PARMS["h_bond_l"])
    h_bond_angle = float(PARMS["h_bond_angle"])

    # WCWC
    WCWC_h_bonds = []
    WCWC_clas = []
    WCWC_conf = []
    WCWC_nucs = []

    # NON_SPEC
    non_spec_h_bonds = []
    non_spec_clas = []
    non_spec_spec = []
    non_spec_nucs = []

    # STACKING
    stacking_nucs = []
    stacking_energies = []
    stacking_pi = {}

    # NUC CHAR
    nucleotides_characteristics = {}
    for i in nucleotides:
        nucleotides_characteristics[i] = [[0], [None]]

    # ALL TO ALL
    for a in range(0, len(nucleotides)):
        for b in range(a+1, len(nucleotides)):
            i = nucleotides[a]
            j = nucleotides[b]
            if_a_WCWC_pair = False
            if whether_not_too_far_simple(i, j, PARMS["cutoff"], TimeTable, N):
                potential_h_bonds = detect_H_bonds(i, j, table,
                                                   h_bond_l, h_bond_angle,
                                                   PARMS, TimeTable, N)
                h_bonds = [t for t in potential_h_bonds
                           if (t[3][0] != -1 and t[3][1] != -1) and
                           t[3][0] >= h_bond_angle and t[3][1] <= h_bond_l]
                if len(h_bonds) > 0:
                    clas = classify_H_bonds(i, j, h_bonds, table, PARMS)
                    tmp = clas.split("/")
                    nucleotides_characteristics[i][0].append(len(h_bonds))
                    nucleotides_characteristics[i][1].append(j)
                    nucleotides_characteristics[j][0].append(len(h_bonds))
                    nucleotides_characteristics[j][1].append(i)
                    if_a_WCWC_pair = True
                    if tmp[0] == tmp[1] == "WC":
                        WCWC_h_bonds.append(h_bonds)
                        WCWC_clas.append(clas)
                        WCWC_conf.append(nucleotides_conformation(i, j,
                                                                  TimeTable, N))
                        WCWC_nucs.append([a, b])
                    else:
                        non_spec_nucs.append([a, b])
                        non_spec_h_bonds.append(h_bonds)
                        non_spec_clas.append(clas)
                        non_spec_spec.append(nucleotides_conformation(i, j,
                                                                      TimeTable, N))
                else:
                    h_bonds_margin = [t for t in potential_h_bonds if
                                      (t[3][0] >= h_bond_angle *
                                       (1.0 - margin) and
                                       t[3][1] <= h_bond_l*(1.0 + margin))]

                    if len(h_bonds_margin) > 0:
                        non_spec_nucs.append([a, b])
                        non_spec_h_bonds.append(h_bonds_margin)
                        non_spec_clas.append("")
                        non_spec_spec.append("Large margin")
                # stacking
                ccvi = coords_and_charges_and_vdw(i, charges, TimeTable, N)
                ccvj = coords_and_charges_and_vdw(j, charges, TimeTable, N)
                if (i != j and not if_a_WCWC_pair and
                    whether_not_too_far_simple(i, j, PARMS["cutoff_stacking"],
                                               TimeTable, N)):
                    base1 = only_base(ccvi)
                    base2 = only_base(ccvj)
                    stacking = vdw_and_coulomb(base1, base2)
                    if stacking[1] <= PARMS["vdw_cutoff_stacking"]:
                        stacking_nucs.append([i, j])
                        stacking_energies.append(stacking)
                # distance from OPs
                stacking_pi.update(pi_stacking(i, j, ccvi, ccvj, PARMS))

    out_dicti = {
        "WCWC_h_bonds": WCWC_h_bonds,
        "WCWC_clas":  WCWC_clas,
        "WCWC_conf": WCWC_conf,
        "WCWC_nucs": WCWC_nucs,
        "non_spec_h_bonds": non_spec_h_bonds,
        "non_spec_clas": non_spec_clas,
        "non_spec_spec": non_spec_spec,
        "non_spec_nucs":  non_spec_nucs,
        "translated_nucleotides":
        translate_nucleotides_characteristics(nucleotides_characteristics),
        "stacking_nucs": stacking_nucs,
        "stacking_energies": stacking_energies,
        "stacking_pi": stacking_pi}

    return out_dicti


def base_atoms(residue):
    # print [at for at in residue if at.get_name]
    out = []
    for at in residue:
        name = at.get_name()
        if "'" not in name and "P" not in name:
            out.append(at)
    return out


def check_origin(resname, PARMS):
    try:
        return PARMS["dicti_modified_nucs"][resname]
    except Exception:
        return False


def translate_nucleotides_characteristics(nucleotides_characteristics):
    out = {}
    for i in nucleotides_characteristics.keys():
        out[RS.resid_description(i)] = [sum(nucleotides_characteristics[i][0]),
                                        " ".join([RS.resid_description(nuc)
                                                  for nuc in
                                                  nucleotides_characteristics[i][1]])]
    return out


# return cis/trans - nucleotide conformation for a bond
def nucleotides_conformation(i, j, TimeTable, N):
    atoms_coords = []
    for n in [i, j]:
        if if_exists(n, "C1'"):
            atoms_coords.append(TimeTable[n["C1'"].get_full_id()][N])
        elif if_exists(n, "C1*"):
            atoms_coords.append(TimeTable[n["C1*"].get_full_id()][N])
        elif if_exists(n, "C8*"): #PNA in Amber
            atoms_coords.append(TimeTable[n["C8*"].get_full_id()][N])
#        elif if_exists(n, "N1*"):
#            atoms_coords.append(TimeTable[n["N1*"].get_full_id()][N])

        elif if_exists(n, "C4'"): #PNA in CHARMM
            atoms_coords.append(TimeTable[n["C4'"].get_full_id()][N])
        else:
            print "No c1/n1 in", n
        if if_exists(n, "N9"):
            atoms_coords.append(TimeTable[n["N9"].get_full_id()][N])
        elif if_exists(n, "N1"):
            atoms_coords.append(TimeTable[n["N1"].get_full_id()][N])
        else:
            print "No n9/n1 in", n

    angle = measure_torsion_angle(atoms_coords)
    if angle > -90 and angle < 90:
        return "Cis"
    else:
        return "Trans"


def get_atoms_of_the_rest(res):
    out = []
    for i in res:
        name = i.get_name()
        if "'" not in name and "*" not in name and "H" not in name:
            out.append(i)
    out.sort()
    return out


def measure_torsion_angle(list_of_arrays):
    v1_2 = list_of_arrays[0] - list_of_arrays[1]
    v3_2 = list_of_arrays[3] - list_of_arrays[1]
    v3_4 = list_of_arrays[3] - list_of_arrays[2]
    il_skal_1_2_3_2 = np.dot(v1_2, v3_2)
    il_skal_3_4_3_2 = np.dot(v3_4, v3_2)
    mian_1 = np.dot(v3_2, v3_2)
    v1 = v1_2 - (il_skal_1_2_3_2/mian_1)*v3_2
    v2 = -v3_4 + (il_skal_3_4_3_2/mian_1)*v3_2
    cos_dih = np.dot(v1, v2/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    #added because of python inaccuracy in the floating point operations
    if cos_dih > 1:
        cos_dih = 1
    else:
        if cos_dih < -1:
            cos_dih = -1
    dih_rad = np.arccos(cos_dih)
    znak_dih = math.copysign(1, np.dot(v1_2, (np.cross(v3_2, v3_4))))
    dih = math.degrees(dih_rad)*int(znak_dih)
    return dih


# Cutoff - simple
def whether_not_too_far_simple(i, j, cutoff, TimeTable, N):
    coords = []
    for n in [i, j]:
        if if_exists(n, "C1'"):
            coords.append(list(TimeTable[n["C1'"].get_full_id()][N]))
        elif if_exists(n, "C1*"):
            coords.append(list(TimeTable[n["C1*"].get_full_id()][N]))
        elif if_exists(n, "C8*"): #PNA in Amber
            coords.append(list(TimeTable[n["C8*"].get_full_id()][N]))
        elif if_exists(n, "C4'"): #PNA in CHARMM
            coords.append(list(TimeTable[n["C4'"].get_full_id()][N]))
#        elif if_exists(n, "N1*"):
#            coords.append(list(TimeTable[n["N1*"].get_full_id()][N]))
#        elif if_exists(n, "N4'"):
#            coords.append(list(TimeTable[n["N4'"].get_full_id()][N]))
        else:
            return False
    if dist(coords[0], coords[1]) < cutoff:
        return True
    return False


def purine_or_pyrimidyne(nuc):
    if names(nuc.get_resname()) in ["A", "G"]:
        return "purine"
    else:
        return "pyrimidyne"


def classify_H_bonds(i, j, h_bonds, table, PARMS):
    atoms_i = [t[0] for t in h_bonds]
    atoms_j = [t[2] for t in h_bonds]
    return "/".join([classify_list_of_atoms(i.get_resname(),
                                            atoms_i, table,
                                            PARMS),
                     classify_list_of_atoms(j.get_resname(),
                                            atoms_j, table,
                                            PARMS),
                     str(len(h_bonds))])


#  * means corner
def classify_list_of_atoms(res_name, atoms, table, PARMS):
    tmp = []
    modi = False
    if if_modified(res_name, PARMS):
        modi = True
        res_name = check_origin(res_name.strip(), PARMS)
    for i in atoms:
        a = what_site_atom(res_name, i, table)
        if a == [] and modi:
            a = ["Modified"]
        tmp.append(a)
    if len(tmp) > 1:
        t = list(reduce(lambda x, y: x & y, (set(i) for i in tmp)))
        if not t:
            return " or ".join(a)
        return t[0]
    else:
        return "*".join(tmp[0])


def what_site_atom(res_name, atom, table):
    row = table[[i[0] for i in table].index(names(res_name))]
    sites = []
    for i in range(0, len(row)):
        if atom in row[i].split():
            tmp = table[0][i].split()[0]
            if tmp not in sites:
                sites.append(tmp)
    return sites


def atom_coords(res, at, TimeTable, N, PARMS):
    at = at.replace(" ", "")
    try:
        return TimeTable[res[at].get_full_id()][N]
    except Exception:
        try:
            at = at.replace("'", "*")
            return TimeTable[res[at].get_full_id()][N]
        except Exception:
            try:
                at = at + "'"
                return TimeTable[res[at].get_full_id()][N]
            except Exception:
		#If and elif below are an ugly solution for two different names for the methyl C in PNA TPN residue in Amber and CHARMM.
		#Will work probably also for the DNA Thymine
    		if PARMS["force_field"] == "AMBER" and "C5M" in at:
			return False
    	        elif PARMS["force_field"] == "CHARMM" and "C7" in at:
                        return False
                else:
                        print "There is no atom ", at, " in resid :",
                        print res.get_resname(), RS.resid_description(res)
                        return False


# Function detecting all possible H-Bonds between two residues
def detect_H_bonds(res1, res2, table, h_bond_l, h_bond_angle,
                   PARMS, TimeTable, N):
    Hbonds = []
    for at_donor in donor_acceptor_list(res1, "d", table, if_modified(res1,
                                                                      PARMS)):
        for at_acceptor in donor_acceptor_list(res2, "a",
                                               table, if_modified(res2, PARMS)):
            for hydrogen in hydrogens_names_for_donor_name(res1, at_donor):
                tmp = if_hbond(atom_coords(res1, at_donor, TimeTable, N, PARMS),
                               TimeTable[res1[hydrogen].get_full_id()][N],
                               atom_coords(res2, at_acceptor, TimeTable, N, PARMS),
                               h_bond_l, h_bond_angle, PARMS)
                Hbonds.append([at_donor, hydrogen, at_acceptor, tmp])

    for at_donor in donor_acceptor_list(res2, "d", table, if_modified(res2,
                                                                      PARMS)):
        for at_acceptor in donor_acceptor_list(res1, "a",
                                               table, if_modified(res1, PARMS)):
            for hydrogen in hydrogens_names_for_donor_name(res2, at_donor):
                tmp = if_hbond(atom_coords(res2, at_donor, TimeTable, N, PARMS),
                               TimeTable[res2[hydrogen].get_full_id()][N],
                               atom_coords(res1, at_acceptor, TimeTable, N, PARMS),
                               h_bond_l, h_bond_angle, PARMS)
                Hbonds.append([at_acceptor, hydrogen, at_donor, tmp])
    return Hbonds


# Function to check whether there is a h_bood
def if_hbond(donor_coord, hydrogen_coord,
             acceptor_coord, h_bond_l, h_bond_angle, PARMS):
    try:
        kat = f_angl([donor_coord, hydrogen_coord, acceptor_coord])
        if PARMS["h_bond_atom"] == "donor":
            dlugosc = dist(donor_coord, acceptor_coord)
        else:
            dlugosc = dist(hydrogen_coord, acceptor_coord)
        return (kat, dlugosc)
    except Exception:
        return (-1, -1)


# Function to compte the distance
def dist(coord1, coord2):
    return ((coord1[0] - coord2[0])**2 +
            (coord1[1] - coord2[1])**2 +
            (coord1[2] - coord2[2])**2)**0.5


# Function to compute the angle
def f_angl(list_of_arrays):
    a = list_of_arrays[0] - list_of_arrays[1]
    b = list_of_arrays[2] - list_of_arrays[1]
    cos = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))

    #added because of python inaccuracy in the floating point operations
    if cos > 1:
        cos = 1
    else:
        if cos < -1:
            cos = -1
    try:
        ang_rad = math.degrees(np.arccos(cos))
    except RuntimeWarning:
        print "Cosinus, wrong ", cos
        ang_rad = math.degrees(ang_rad)
    return ang_rad


# Checking if the atom exists
def if_exists(res, at_name):
    try:
        res[at_name].get_coord()
        return True
    except Exception:
        return False


# Creating list of the hydrogensnames names
def hydrogens_names_for_donor_name(res, at_donor):
    inki = ["1", "2", "3"]
    hs = []
    k = at_donor[2:]
    # pusty
    h = "H" + at_donor[1] + k
    if if_exists(res, h) and h not in hs:
        hs.append(h)
    for i in inki:
        # beginning
        h = i + "H" + at_donor[1] + k
        if if_exists(res, h) and h not in hs:
            hs.append(h)
        # ending
        h = "H" + at_donor[1] + i + k
        if if_exists(res, h) and h not in hs:
            hs.append(h)
    return hs


# function gives list of atoms  - acceptors or donors
def donor_acceptor_list(resid, da, table, modify):
    if not modify:
        n = names(resid.get_resname())
	print "name ", n, "resid ", resid
	
#	for el in table:
#		print el[0]
        
	t = table[[i[0] for i in table].index(n)]
        tmp = ""
        if da == "a":
            for i in a:
                tmp += " " + t[table[0].index(i)]
        if da == "d":
            for i in d:
                tmp += " " + t[table[0].index(i)]
        # PODMIENIC
        out = []
        for i in tmp.split():
            if i not in out:
                out.append(i)
        return out

    elif modify:
        atoms_base = base_atoms(resid)
        hydrogens = [i for i in atoms_base if "H" in i.get_name()]
        nums = [j.get_name()[1] for j in hydrogens]
        out = []
        for n in [i.get_name() for i in atoms_base]:
            if "H" not in n and len(n) > 1:
                if da == "d" and list(n)[1] in nums:
                    out.append(n)
                elif da == "a" and list(n)[1] not in nums:
                    out.append(n)
        return out


def number_first(s):
    num = [i for i in s if i.isdigit()]
    letters = [i for i in s if not i.isdigit()]
    return "".join(num)+"".join(letters)


# changing name for residue
def names(resname):
    resname = resname.replace(" ", "")
    if resname in ["G", "U", "A", "C"]:
        return resname
    elif resname in ["GUA", "URA", "CYT", "ADE",
                     "G5", "G3", "A5", "A3", "C3",
                     "C5", "U5", "U3", "T5", "T3",
                     ""]:
        return resname[0]
    elif resname in ["DG5", "DG3", "DA5", "DA3", "DC3",
                     "DC5", "DU5", "DU3", "DT5", "DT3"]:
        return resname[:2]
    elif (resname.startswith("R") or
          resname.startswith("L")):
        return resname[1]
    else:
        return resname


# margin in %
def chcek_if_close_to_value(x, ref, margin):
    if (ref - ref*(1-margin) < x) and (ref*(1 + margin) > x):
        return True
    else:
        return False


# save to csv file
def list_to_csv_file(file_name, list_out):
    f = open(file_name + ".csv", "w")
    for i in range(0, len(list_out)):
        f.write(str(i) + " ; " + str(list_out[i]) + "\n")
    f.close()
    print "Saved csv file " + file_name + ".csv in directory: "


def nicely_print(l, LINE_LENGTH):
    out = ""
    word = " ".join(l)
    if len(l) > LINE_LENGTH:
        number_of_lines = len(l)/LINE_LENGTH
        for i in range(0, number_of_lines):
            out += str(i*LINE_LENGTH + 1) +\
                " "*(8 - len(str(i*LINE_LENGTH + 1))) +\
                word[i*LINE_LENGTH:(i + 1)*LINE_LENGTH] +\
                "   " + str((i + 1)*LINE_LENGTH) + "\n"
        # Adding ending
        last_one = number_of_lines*LINE_LENGTH+1
        out += str(last_one) + " "*(8-len(str(last_one))) +\
            word[last_one:] + " " + str(len(word)+1)+"\n"
    else:
        out += "1    " + word + "  " + str(len(l))
    return out


# @timer
def create_reps(l, color_number):
    trip = ""
    if l:
        for i in l:
            if i and "vmd" not in i:
                try:
                    tmp = i.split("]")[1]
                except Exception:
                    tmp = i
                tt = tmp.split("chain")[1:]
                if len(tt) == 1:
                    tt = "chain " + tt[0]
                else:
                    tt = "or ".join(["chain " + j for j in tt])
                trip += "mol selection {" + tt +\
                    "}\nmol color ColorID " +\
                    str(color_number) +\
                    "\nmol rep Licorice \nmol addrep top\n"
        return trip
    return " "


# @timer
def vmd_run(vmd, filename, PARMS):
    f = open(filename, "w")
    text = """mol load pdb """ + os.path.realpath(PARMS["file_name"]) + """
mol delrep 0 top
mol rep Lines
mol color colorID 2
mol addrep top \n""" + create_reps(vmd[0], 4) + create_reps(vmd[1], 5) +\
        create_reps(vmd[2], 6) + create_reps(vmd[3], 7)
    f.write(text)
    f.close()


def stacking_presentation(stacking_nucs, stacking_energies):
    out = " qn1 n2     Coul     VDW    sum\n"
    per_nucleotide = {}
    for i, n in enumerate(stacking_nucs):
        m = RS.resid_description(n[0]) + " " + RS.resid_description(n[1])
        add_safely_to_dicti(per_nucleotide,
                            RS.resid_description(n[0]),
                            stacking_energies[i])
        add_safely_to_dicti(per_nucleotide,
                            RS.resid_description(n[1]),
                            stacking_energies[i])
        if not IF_MODIFIED[n[0].get_id()] and not IF_MODIFIED[n[1].get_id()]:
            out += m + " "*(10 - len(m)) + str(stacking_energies[i]) + "\n"
        else:
            out += m + " "*(10 - len(m)) + " MODIFIED NUCLEOTIDES" + "\n"
    return (out, per_nucleotide)


def add_safely_to_dicti(dicti, k, v):
    if k in dicti.keys():
        tmp = dicti[k]
        dicti[k] = add_tuples(tmp, v)
    else:
        dicti[k] = v


def add_tuples(t1, t2):
    n = []
    for i in range(min(len(t1), len(t2))):
        n.append(round(t1[i]+t2[i], 2))
    return tuple(n)


# @timer
# file_name, cutoff, cutoff_stacking, h_bond_l, h_bond_angle, vmd
def single_frame_analysis(PARMS):
    nucleotides = get_nucleic_from_pdb(PARMS)
    if len(nucleotides) <= 1:
        PARMS["OUT_FILE"].write("Are there any nucleotides in the structure?\n")
        sys.exit()
    charges = read_in_charges(nucleotides, PARMS)
    TimeTable = TimeTable_from_nucleotides(nucleotides)
    megalist = measure_for_all(nucleotides, charges, PARMS, TimeTable, 0)
    struct_analysis = RS.analysis(megalist, nucleotides, True, IF_MODIFIED)
    filename = PARMS["out_name"] + "_vmd_run.tcl"

    out_file = (struct_analysis["fileText"] + "\n\nStacking pairs:" +
                stacking_presentation(megalist["stacking_nucs"],
                                      megalist["stacking_energies"])[0] +
                "\n\nStacking pi pairs\n" +
                stacking_pi_description(megalist["stacking_pi"]))

    vmd_run(struct_analysis["vmd"], filename, PARMS)
    if PARMS["vmd"]:
        print "Trying to run VMD ... "
        os.system("vmd -e "+filename)
        print " "

    LISTS = {}
    for i in ["arches", "dot_bracket", "trajectory", "translated_nucleotides",
              "stacking_pi", "pairs_in_time", "stacking"]:
        LISTS[i] = {}

    traj_len = 1
    frames = {1: (megalist, struct_analysis, out_file)}
    out_file_text = ""
    for key in frames.keys():
        frame = frames[key]
        out_file_text += frame[2] + "\n"
        for i in ["arches", "dot_bracket", "trajectory"]:
            LISTS[i][key] = frame[1][i]
        for i in ["translated_nucleotides", "stacking_pi"]:
            LISTS[i][key] = frame[0][i]
        LISTS["stacking"][key] = (frame[0]["stacking_nucs"],
                                  frame[0]["stacking_energies"])
        LISTS["pairs_in_time"][key] = frame[0]
    create_outputs(LISTS, traj_len, nucleotides,
                   struct_analysis["dot_bracket"], PARMS)
    PARMS["OUT_FILE"].write(out_file_text)
    return LISTS


def stacking_pi_description(l):
    out = ""
    for i in l:
        out += i + " " + " ".join([str(m) for m in l[i]]) + "\n"
    return out


# @timer
def single_frame_in_traj_analysis(nucleotides, charges,
                                  PARMS, TimeTable, N):
    megalist = measure_for_all(nucleotides, charges, PARMS, TimeTable, N)
    struct_analysis = [None]*10
    try:
        struct_analysis = RS.analysis(megalist, nucleotides,
                                      False, IF_MODIFIED)
        out_file = struct_analysis["fileText"] +\
            "\n\nStacking pairs:" +\
            stacking_presentation(megalist["stacking_nucs"],
                                  megalist["stacking_energies"])[0] +\
            "\n\nStacking pi pairs\n" +\
            stacking_pi_description(megalist["stacking_pi"])
    except Exception:
        f = open("error"+str(N), "w")
        f.write("problem in single_frame_in_traj")
        f.close()
        out_file = "Problem with frame"
    return (megalist, struct_analysis, out_file)


# @timer
def for_a_sub_traj(universe, nucleotides, index, CAT, NAME, charges,
                   PARMS, TimeTable):
    d = {}
    number_of_pickles = 0
    num_of_frames = float(len(index))
    for i, N in enumerate(index):
        print "frame number ", N, " "*(10-len(str(a))),\
            100*round((i+1)/num_of_frames, 3), "%\n"

        d[N] = single_frame_in_traj_analysis(nucleotides, charges,
                                             PARMS, TimeTable, N)
        if float(sys.getsizeof(d)/1024**3) > PARMS["max_memory_GB"]:
            pickle_file_name = (CAT + "/" + NAME + "_out_dictionary_pickle_" +
                                str_index(index) +
                                "_" + str(number_of_pickles))
            print "additional Pickle dumping " + pickle_file_name
            pickle.dump(d, open(pickle_file_name, "a"))
            number_of_pickles += 1
            del d  # delete dictionary
            d = {}
    if d:
        pickle_file_name = (CAT + "/" + NAME + "_out_dictionary_pickle_" +
                            str_index(index) + "_" +
                            str(number_of_pickles))
        print "Pickle dumping " + pickle_file_name
        pickle.dump(d, open(pickle_file_name, "a"))
    print "Done ", str_index(index)


def str_index(index):
    try:
        return str(index[0]) + "_"+str(index[-1])
    except Exception:
        return str(index)


def test(index):
    out = []
    for i in range(index[0], index[1]):
        out.append(i)
    return out


def filips_format(dot_brackets, fname):
    f = open(fname, "w")
    ll = [dot_brackets[i] for i in dot_brackets]
    f.write("\n".join(ll))
    f.close()


def divide_trajectory(PARMS):
    should = range(PARMS["first_frame"], PARMS["last_frame"],
                   PARMS["stride"])
    num_of_frames = len(should)
    # last frame
    trajs = []
    if PARMS["threads"] > 1:
        if PARMS["threads"] < num_of_frames:
            large_step = int(float(num_of_frames)/(PARMS["threads"]))
        else:
            PARMS["threads"] = num_of_frames
            large_step = 1
        for j in range(0, PARMS["threads"]-1):
            trajs.append(should[j*large_step:(j+1)*large_step])
        trajs.append(should[(j+1)*large_step:])
    else:
        trajs = [should]
    return (trajs, num_of_frames, should)


def get_a_catalogue(PARMS):
    NAME = PARMS["out_name"].split("/")[-1]
    ccat = "/".join(PARMS["out_name"].split("/")[:-1])
    cat = ccat+"/"+NAME+"_out_dictionaries_MINT"
    if not PARMS["only_analysis"]:
        os.system("rm -r " + cat)
        print "creating catalogue " + cat
        os.system("mkdir "+cat)
    return (cat, NAME)


def exact_description(at, rnum, segid, resname, atname):
    if segid.strip():
        return "nucleic and segid " + segid + " and resid " + rnum + \
            " and resname " + resname + " and name \"" + atname + "\""
    else:
        return "nucleic and resid " + rnum + " and resname " + resname +\
            " and name \'" + atname + "\'"


def TimeTable_from_nucleotides(nucleotides):
    TimeTable = {}
    for nuc in nucleotides:
        for at in nuc:
            full_id = at.get_full_id()
            TimeTable[full_id] = [at.get_coord()]
    return TimeTable


def ReadInTrajectory(nucleotides, universe, PARMS):
    TimeTable = {}
    for nuc in nucleotides:
        segid = nuc.get_segid().strip()
        resname = nuc.get_resname().strip()
        rnum = str(nuc.get_id()[1])
        for at in nuc:
            atname = at.get_id()
            full_id = at.get_full_id()
            desc = exact_description(at, rnum, segid, resname, atname)
            atoms = universe.select_atoms(desc)
            if len(atoms) == 0:
                desc = "resid "+rnum
                tmp = [n for n in universe.select_atoms(desc)
                       if (n.name == atname) and (n.resname == resname)]
                atoms = MDAnalysis.core.AtomGroup.AtomGroup(tmp)

            if len(atoms) == 1:
                try:
                    timeseries = universe.trajectory.timeseries(atoms)
                    TimeTable[full_id] = timeseries[0]
                except Exception:
                    TimeTable[full_id] = ManualTime_series(atoms[0].number,
                                                           universe)

            else:
                print "trajectory PROBLEM with ", at, nuc
                print desc
                print "atoms", [a for a in atoms]
    return TimeTable


def ManualTime_series(atomNumber, universe):
    return [np.array(ts[atomNumber].tolist())
            for ts in universe.trajectory]


# @timer
def trajectories(PARMS):
    CAT, NAME = get_a_catalogue(PARMS)
    nucleotides = get_nucleic_from_pdb(PARMS)
    if len(nucleotides) <= 1:
        PARMS["OUT_FILE"].write("Are there any nucleotides in the structure?\n")
        sys.exit()
    charges = read_in_charges(nucleotides, PARMS)
    universe = MDAnalysis.Universe(PARMS["file_name"], PARMS["file_dcd"])
    TimeTable = ReadInTrajectory(nucleotides, universe, PARMS)

    if PARMS["last_frame"] == -1:
        PARMS["last_frame"] = len(universe.trajectory)

    trajs, num_of_frames, should = divide_trajectory(PARMS)
    print "\n\nTrajectory has ", len(universe.trajectory),
    print " frames, running for ", num_of_frames
    if not PARMS["only_analysis"]:
        ths = []
        ths = [multiprocessing.Process(target=for_a_sub_traj,
                                       args=(universe, nucleotides,
                                             trajs[i], CAT, NAME, charges,
                                             PARMS, TimeTable))
               for i in range(PARMS["threads"])]
        for p in ths:
            p.start()

        for p in ths:
            p.join()
    frames = {}
    for i in trajs:
        print "Reading pickles .. ", str_index(i),
        ind = 0
        while True:
            try:
                d1 = pickle.load(open(CAT + "/" + NAME +
                                      "_out_dictionary_pickle_" +
                                      str_index(i) + "_" + str(ind), "r"))
                ind += 1
                print len(d1.keys()),
                frames.update(d1)
                print len(frames.keys())
            except Exception:
                break

    frames.keys().sort()
    if len(frames.keys()) != num_of_frames:
        print "Some Frames are missing!!"
        sys.exit()

    LISTS = {}
    for i in ["arches", "dot_bracket", "trajectory", "translated_nucleotides",
              "stacking_pi", "pairs_in_time", "stacking"]:
        LISTS[i] = {}
    out_file_text = ""
    traj_len = float(len(frames.keys()))
    for key in frames.keys():
        frame = frames[key]
        out_file_text += "\nframe number " + str(key) + "\n" + frame[2] + "\n"
        for i in ["arches", "dot_bracket", "trajectory"]:
#            print frame[1][i]
            LISTS[i][key] = frame[1][i]
        for i in ["translated_nucleotides", "stacking_pi"]:
            LISTS[i][key] = frame[0][i]
        LISTS["stacking"][key] = (frame[0]["stacking_nucs"],
                                  frame[0]["stacking_energies"])
        LISTS["pairs_in_time"][key] = frame[0]

    XML.RNAStructMl(PARMS["out_name"]+"_RNAStructMl.xml",
                    LISTS["dot_bracket"],
                    RS.get_sequence(nucleotides, IF_MODIFIED),
                    PARMS["file_name"])
    filips_format(LISTS["dot_bracket"], PARMS["out_name"] + "_dot_bracket.txt")
    avg_dot_bracket = average_secondary_structure(LISTS["arches"])
    PARMS["OUT_FILE"].write("\nAverage secondary structure :\n" +
                            RS.get_sequence(nucleotides, IF_MODIFIED) +
                            "\n" + avg_dot_bracket + "\n" + out_file_text)
    create_outputs(LISTS, traj_len, nucleotides, avg_dot_bracket, PARMS)


def create_outputs(LISTS, traj_len, nucleotides, avg_dot_bracket, PARMS):
    tmp = motifs_in_time(LISTS["trajectory"])
    per_nuc, out_per_nuc = per_nucleotide_char(LISTS["pairs_in_time"],
                                               traj_len, nucleotides, PARMS)

    OUTPUT_LISTS = {"pairing":  multi_to_dicti(LISTS["translated_nucleotides"]),
                    "nucleotides_characteristics": out_per_nuc,
                    "pairs": histogram_of_pairs(LISTS["pairs_in_time"],
                                                traj_len, nucleotides, PARMS),
                    "helices": helices_in_time_to_csv(tmp["helices"],
                                                      traj_len, PARMS),
                    "pseudoknots": pseudo_in_time_to_csv(tmp["pseudo"],
                                                    traj_len, PARMS),
                    "triplets": triplex_in_time_to_csv(tmp["triplexes"],
                                                      traj_len, PARMS),
                    "motifs": motifs_in_time_to_csv(tmp["motifs"],
                                                    traj_len, PARMS),
                    "pi_interactions": stacking_pi_in_time_to_csv(LISTS["stacking_pi"],
                                                                  traj_len, PARMS),
                    "stacking": stacking_in_time_to_csv(LISTS["stacking"],
                                                        traj_len, PARMS)}

    if PARMS["Mode"] == "Traj":
        maximal_clusters = MC.cluster(PARMS["margin"], traj_len,
                                      OUTPUT_LISTS["motifs"],
                                      PARMS["time_cutoff"])
        OUTPUT_LISTS["motifs_clusters"] = maximal_clusters[2]
        OUTPUT_LISTS["average_motifs"] = return_maximal_clusters(maximal_clusters, per_nuc)

    cX.list_to_xls(OUTPUT_LISTS, PARMS["out_name"], PARMS)
    if PARMS["create_csvs"]:
        cX.create_csvs(OUTPUT_LISTS, PARMS)
    parms_varna = {"nucleotides_eval": per_nuc,
                   "sequence": RS.get_sequence(nucleotides, IF_MODIFIED),
                   "dot_bracket": avg_dot_bracket,
                   "output": PARMS["out_name"] + "_varna.html",
                   "home": PARMS["MINT_home"],
                   "mode": PARMS["Mode"],
                   "output_pictures": PARMS["out_name"] + "_"}
    Varna.run(parms_varna)
    replace_beta_occup_with_vals(per_nuc, nucleotides, PARMS)


def average_secondary_structure(arches):
    all_pairs = []
    for key in arches:
        all_pairs.extend(pairs_from_arches(arches[key]))
    c = Counter(all_pairs)
    dot_bracket = ["."]*len(arches[key])
    for i in list(set(all_pairs)):
        if c[i] > 0.5*len(arches):
            dot_bracket[min(i)] = "("
            dot_bracket[max(i)] = ")"
    return "".join(dot_bracket)


def pairs_from_arches(single_list):
    pairs = []
    for i, n in enumerate(single_list):
        if n and set([i, n]) not in pairs:
            pairs.append(set([i, n]))
    pairs = [tuple(i) for i in pairs]
    return pairs


def return_maximal_clusters(maximal_cluster, per_nuc):
    licznik = 0
    out = []
    l = maximal_cluster[0]
    maximal_cluster = maximal_cluster[1]
    for i in maximal_cluster:
        ttmp = [licznik]
        ttmp.extend(l[i])
        out.append(ttmp)
        wc = ""
        non_wc = ""
        for m in l[i][1].split("-")[:-1]:
            wc += str(round(per_nuc[m][1], 1))+"  "
            non_wc += str(round(per_nuc[m][2], 1)) + "  "
        out.append(["WC", wc])
        out.append(["non-WC", non_wc])
        licznik += 1
    return out


def get_nucleic_from_its_description(nucleotides, description):
    descriptions = [RS.resid_description(i) for i in nucleotides]
    try:
        return nucleotides[descriptions.index(description)]
    except Exception:
        return False


def avg(ll):
    number_of_frames = len(ll)
    return round(float(sum(ll))/float(number_of_frames), 2)


def sd(ll, aa):
    ff = float(len(ll))
    return round(math.sqrt(sum([(i-aa)**2 for i in ll])/ff), 2)


def per_nucleotide_char(pairs_in_time, traj_len,
                        nucleotides, PARMS):
    nucls = {}
    for i in nucleotides:
        nucls[RS.resid_description(i)] = {"WC": [0]*int(traj_len),
                                          "non-WC": [0]*int(traj_len),
                                          "total Hbonds": [0]*int(traj_len),
                                          "Coulomb": [0]*int(traj_len),
                                          "VDW": [0]*int(traj_len),
                                          "Ssum": [0]*int(traj_len)}
    for i, ts in enumerate(pairs_in_time.keys()):
        # H - bonding
        l_pair = pairs_in_time[ts]["WCWC_nucs"] +\
            pairs_in_time[ts]["non_spec_nucs"]
        l_type = pairs_in_time[ts]["WCWC_clas"] +\
            pairs_in_time[ts]["non_spec_clas"]
        for n in range(len(l_pair)):
            tmp = [RS.resid_description(nucleotides[l_pair[n][0]]),
                   RS.resid_description(nucleotides[l_pair[n][1]])]
            typ = l_type[n].split("/")

            if len(typ) > 2 and typ[0] == "WC" and typ[1] == "WC":
                for t in tmp:
                    nucls[t]["WC"][i] = float(typ[2])
            elif len(typ) > 2:
                for t in tmp:
                    nucls[t]["non-WC"][i] = float(typ[2])

        for n in nucls.keys():
            nucls[n]["total Hbonds"][i] = nucls[n]["WC"][i] +\
                nucls[n]["non-WC"][i]
        # Stacking
        sta = stacking_presentation(pairs_in_time[ts]["stacking_nucs"],
                                    pairs_in_time[ts]["stacking_energies"])[1]
        for n in sta.keys():
            for ind, cvs in enumerate(["Coulomb", "VDW", "Ssum"]):
                nucls[n][cvs][i] = sta[n][ind]

        for n in nucls.keys():
            nucls[n]["total Hbonds"][i] = nucls[n]["WC"][i] +\
                nucls[n]["non-WC"][i]
        # Stacking
        sta = stacking_presentation(pairs_in_time[ts]["stacking_nucs"],
                                    pairs_in_time[ts]["stacking_energies"])[1]
        for n in sta.keys():
            for ind, cvs in enumerate(["Coulomb", "VDW", "Ssum"]):
                nucls[n][cvs][i] = sta[n][ind]

    nucls_ev = {}
    for nuc in nucls:
        tmp = []
        ssds = []
        for kk in ["WC", "non-WC", "total Hbonds", "Coulomb", "VDW"]:
            w = nucls[nuc][kk]
            if PARMS["Mode"] == "Traj":
                aa = avg(w)
                ssd = sd(w, aa)
                if kk in ["Coulomb", "VDW"]:
                    ssds.append(ssd)
                tmp.extend([aa, ssd])
            else:
                if len(nucls[nuc][kk]) <= 1:
                    tmp.append(nucls[nuc][kk][0])
                else:
                    print "PROBLEM!!"
                    sys.exit()
        # special treatament of sum:
        if PARMS["Mode"] == "Traj":
            aa = avg(nucls[nuc]["Ssum"])
            ssd = round(math.sqrt(sum([ss**2 for ss in ssds])), 2)
            tmp.extend([aa, ssd])
        else:
            tmp.append(nucls[nuc]["Ssum"][0])
        nucls_ev[nuc] = tmp
    out = []
    for nuc in nucls_ev:
        tmp = [int(nuc.split(":")[1]), nuc]
        tmp.extend(nucls_ev[nuc])
        tmp.append(RS.for_vmd([nuc]))
        out.append(tmp)
    out = sorted(out, key=itemgetter(0))
    return (nucls_ev, out)


# [WCWC_h_bonds, WCWC_clas, WCWC_conf, WCWC_nucs, non_spec_h_bonds,
# non_spec_clas, non_spec_spec, non_spec_nucs,
# translate_nucleotides_characteristics(nucleotides_characteristics)]
# @timer
def histogram_of_pairs(pairs_in_time, traj_len, nucleotides, PARMS):
    pairs = {}
    for ts in pairs_in_time.keys():
        l_pair = pairs_in_time[ts]["WCWC_nucs"] +\
            pairs_in_time[ts]["non_spec_nucs"]
        l_type = pairs_in_time[ts]["WCWC_clas"] +\
            pairs_in_time[ts]["non_spec_clas"]
        l_conf = pairs_in_time[ts]["WCWC_conf"] +\
            pairs_in_time[ts]["non_spec_spec"]
        for i in range(len(l_pair)):
            pair = "/".join([RS.resid_description(nucleotides[l_pair[i][0]]),
                             RS.resid_description(nucleotides[l_pair[i][1]])])
            new_key = pair + "," + "/".join(l_type[i].split("/")[:2]) + "," +\
                l_conf[i]

            if new_key not in pairs.keys():
                pairs[new_key] = [str(ts)]
            else:
                pairs[new_key].append(str(ts))

    out = []
    for key in pairs:
        tmp = key.split(",")[0].split("/")
        nuc = int(tmp[0].split(":")[1])
        if PARMS["Mode"] == "Traj":
            ttmp = [nuc]
            ttmp.extend(key.split(","))
            ttmp.extend([RS.for_vmd(tmp),
                        float(len(pairs[key]))/traj_len,
                        shorten_frmes_list(pairs[key])])
            out.append(ttmp)
        else:
            ttmp = [nuc]
            ttmp.extend(key.split(","))
            ttmp.append(RS.for_vmd(tmp))
            out.append(ttmp)
    if PARMS["Mode"] == "Traj":
        out = sorted(out, key=itemgetter(-2), reverse=True)
    else:
        out = sorted(out, key=itemgetter(0))
    return out


# @timer
def motifs_in_time(out):
    dicti = {"helices": {},
             "pseudo": {},
             "triplexes": {},
             "motifs": {}}

    for ts in out.keys():
        for n, k in enumerate(["helices", "pseudo",
                               "triplexes", "motifs"]):
            for motif in out[ts][n].split("\n")[2:]:
                tmp = motif.split("]")
                if len(tmp) > 1:
                    update_dictionary(dicti[k], tmp[1], ts)
    return dicti


# @timer
def triplex_in_time_to_csv(dicti, traj_len, PARMS):
    out = []
    for keys in dicti:
        resid = [i.replace(" ", "")[1:] for i in keys.split("-")
                 if i.replace(" ", "")[1:] != ""]
        if PARMS["Mode"] == "Traj":
            out.append([keys, RS.for_vmd(resid),
                        float(len(dicti[keys]))/traj_len,
                        shorten_frmes_list(dicti[keys])])
        else:
            out.append([keys, RS.for_vmd(resid)])
    out = sorted(out, key=itemgetter(-2), reverse=True)
    return out


# @timer
def pseudo_in_time_to_csv(dicti, traj_len, PARMS):
    out = []
    for keys in dicti:
        resid = [i[1:] for i in keys.split("-")]
        if PARMS["Mode"] == "Traj":
            out.append([keys, RS.for_vmd(resid),
                        float(len(dicti[keys]))/traj_len,
                        shorten_frmes_list(dicti[keys])])
        else:
            out.append([keys, RS.for_vmd(resid)])
    out = sorted(out, key=itemgetter(-2), reverse=True)
    return out


# @timer
def helices_in_time_to_csv(dicti, traj_len, PARMS):
    out = []
    for keys in dicti:
        beg = [i.split("-") for i in keys.split("->")]
        tmp = []
        for a in [0, 1]:
            tmp.extend(beg[a])
        tmp = RS.helices_to_vmd(tmp)
        if PARMS["Mode"] == "Traj":
            out.append([keys, tmp, float(len(dicti[keys]))/traj_len,
                        shorten_frmes_list(dicti[keys])])
        else:
            out.append([keys, tmp])
    out = sorted(out, key=itemgetter(-2), reverse=True)
    return out


# @timer
def motifs_in_time_to_csv(dicti, traj_len, PARMS):
    l = []
    for keys in dicti:
        tmp = keys.split()
        resids = [i for i in tmp[1].split("-") if i]
        if PARMS["Mode"] == "Traj":
            tmp.extend([RS.for_vmd(resids),
                        float(len(dicti[keys]))/traj_len,
                        shorten_frmes_list(dicti[keys])])
            l.append(tmp)
        else:
            tmp.append(RS.for_vmd(resids))
            l.append(tmp)
    l = sorted(l, key=itemgetter(-2), reverse=True)
    return l


def stacking_in_time_to_csv(stacking, traj_len, PARMS):
    dicti = {}
    for frame_number in stacking.keys():
        frame = stacking[frame_number]
        for pair_number, pair in enumerate(frame[0]):
            key = RS.resid_description(pair[0]) + "_" +\
                RS.resid_description(pair[1])
            if key in dicti.keys():
                for i in range(3):
                    dicti[key][i].append(frame[1][pair_number][i])
                dicti[key][3].append(frame_number)  # frame
            else:
                dicti[key] = [[frame[1][pair_number][0]],
                              [frame[1][pair_number][1]],
                              [frame[1][pair_number][2]],
                              [frame_number]]
    out = []
    for key in dicti:
        pair = dicti[key]
        tmp = []
        sds = []
        for i in range(2):
            if PARMS["Mode"] == "Traj":
                aa = avg(pair[i])
                ssd = sd(pair[i], aa)
                sds.append(ssd)
                tmp.extend([aa, ssd])
            else:
                tmp.extend(pair[i])
        # sum
        if PARMS["Mode"] == "Traj":
            aa = avg(pair[2])
            ssd = round(math.sqrt(sum([i**2 for i in sds])), 2)
            tmp.extend([aa, ssd])
        else:
            tmp.extend(pair[2])

        vmd = RS.for_vmd(key.split("_"))
        nuc = int(key.split("_")[0].split(":")[1])
        ttmp = [nuc, key.replace("_", "/")]
        ttmp.extend(tmp)

        if PARMS["Mode"] == "Traj":
            perc_frames = round(float(len(pair[3]))/traj_len, 2)
            ttmp.extend([vmd, perc_frames, shorten_frmes_list(pair[-1])])
        else:
            ttmp.extend([vmd])
        out.append(ttmp)
    if PARMS["Mode"] == "Traj":
        out = sorted(out, key=itemgetter(-2), reverse=True)
    else:
        out = sorted(out, key=itemgetter(0))
    return out


def stacking_pi_in_time_to_csv(stacking_pi, number_of_frames, PARMS):
    dicti = {}
    for frame_number in stacking_pi.keys():
        frame = stacking_pi[frame_number]
        for key in frame:
            if key in dicti.keys():
                dicti[key][0].append(frame[key][0])  # distances
                for i in range(1, 4):
                    dicti[key][i].append(frame[key][1][i-1])
                dicti[key][4].append(frame_number)   # frames
            else:
                dicti[key] = [[frame[key][0]],
                              [frame[key][1][0]],
                              [frame[key][1][1]],
                              [frame[key][1][2]],
                              [frame_number]]
    out = []
    for key in dicti:
        pair = dicti[key]
        tmp = []
        sds = []
        for i in range(3):
            if PARMS["Mode"] == "Traj":
                aa = avg(pair[i])
                ssd = sd(pair[i], aa)
                tmp.extend([aa, ssd])
                if i == [1, 2]:
                    sds.append(ssd)
            else:
                tmp.extend(pair[i])
        # sum
        if PARMS["Mode"] == "Traj":
            aa = avg(pair[3])
            ssd = round(math.sqrt(sum([i**2 for i in sds])), 2)
            tmp.extend([aa, ssd])
        else:
            tmp.extend(pair[3])

        vmd = RS.for_vmd(key.split("_")[:-1])
        nuc = int(key.split("_")[0].split(":")[1])

        if PARMS["Mode"] == "Traj":
            perc_frames = str(round(float(len(pair[4]))/number_of_frames, 2))
            ttmp = [nuc, key.replace("_", "/")]
            ttmp.extend(tmp)
            ttmp.extend([vmd, perc_frames, shorten_frmes_list(pair[4])])
            out.append(ttmp)
        else:
            ttmp = [nuc, key.replace("_", "/")]
            ttmp.extend(tmp)
            ttmp.append(vmd)
            out.append(ttmp)
    out = sorted(out, key=itemgetter(-2), reverse=True)
    return out


def update_dictionary(dicti, new_key, ts):
    if new_key not in dicti.keys():
        dicti[new_key] = [str(ts)]
    else:
        dicti[new_key].append(str(ts))


# @timer
def multi_to_dicti(nucleotide_characteristics):
    stats = {}
    k = nucleotide_characteristics.keys()[0]
    for nuc in nucleotide_characteristics[k]:
        tmp = []
        for key in nucleotide_characteristics.keys():
            i = nucleotide_characteristics[key]
            a = [str(j).replace(" ", "") for j in i[nuc]]
            tmp.append("-".join(a))
            stats[nuc] = [x + "-" + str(tmp.count(x)) for x in set(tmp)]

    out = []
    for nuc in stats:
        tmp = [int(nuc.split(":")[1]), nuc]
        tmp.extend(stats[nuc])
        out.append(tmp)
    out = sorted(out, key=itemgetter(0))
    return out


def name_of_atom(name):
    if len(name) == 3:
        return "  "+name+" "
    elif len(name) == 2:
        return "  "+name+"  "
    elif len(name) == 1:
        return "   "+name+"  "
    elif len(name) == 4:
        return " "+name+" "


def replace_beta_occup_with_vals(d, nucleotides, PARMS):
    # Use: beta_csv.py file.in file.pdb new.pdb
    # occup - WC-WC; beat = others
    if PARMS["Mode"] == "Traj":
        indeces = ["Hbonds-WC", "Hbonds-WC sd", "Hbonds-non-WC", "Hbonds-non-WC sd",
                   "Hbonds-total", "Hbonds-total sd", "Stacking-Coulomb", "Stacking-Coulomb sd",
                   "Stacking-VDW", "Stacking-VDW sd", "Stacking-total", "Stacking-total sd"]
    else:
        indeces = ["Hbonds-WC", "Hbonds-non-WC",
                   "Hbonds-total","Stacking-Coulomb",
                   "Stacking-VDW", "Stacking-total"]

    for ind, name in enumerate(["Hbonds-WC", "Hbonds-non-WC",
                   "Hbonds-total","Stacking-Coulomb",
                   "Stacking-VDW", "Stacking-total"]):
        f_pdb_out = open(PARMS["out_name"] + "_" + name + ".pdb", "w")
        snum = 1
        res_num = 1
        for nuc in nucleotides:
            nn = nuc.get_resname()
            res_name = nn + " "*(3 - len(nn))
            res_nnum = " "*(6 - len(str(round(res_num, 3)))) + str(res_num)
            ch = str(nuc.get_parent().get_id())
            beta = str(d[RS.resid_description(nuc)][indeces.index(name)])
            for atom in nuc:
                cc = "".join([" "*(8 - len(str(i)))+str(i)
                              for i in atom.get_coord()])
                nnum = " "*(5-len(str(snum))) + str(snum)
                at = name_of_atom(atom.get_name())
                l = "ATOM  " + nnum + at + res_name + " " + ch +\
                    res_nnum + "    " + cc + "  0.00" +\
                    " "*(6-len(beta))+beta+"\n"
                snum += 1
                f_pdb_out.write(l)
            res_num += 1
        print "New pdb file " +\
            PARMS["out_name"] + "_" + name + ".pdb" +\
            " created colored by " + name


def shorten_frmes_list(l):
    i = 0
    out = ""
    l.sort()
    while i < len(l):
        out += " " + str(l[i])
        if i + 1 < len(l) and int(l[i]) + 1 == int(l[i + 1]):

            while i + 1 < len(l) and int(l[i]) + 1 == int(l[i + 1]):
                i += 1
            out += "->" + str(l[i])
        i += 1
    return out


def print_parameters(PARMS):
    running = ""
    excluded = ["table_nuc_read", "dicti_modified_nucs", "OUT_FILE",
                "out_name_old", "out_file_name_general"]

    if PARMS["Mode"] == "Single":
        excluded.extend(["file_dcd", "first_frame", "last_frame",
                         "threads", "time_cutoff"])
    ll = PARMS.keys()
    ll.sort()
    running += "\n" + "MINT_home = " + PARMS["MINT_home"]
    i = ll.index("MINT_home")
    del ll[i]
    for i in ll:
        if i not in excluded:
            running += "\n" + i + "=" +\
                str(PARMS[i]).replace(PARMS["MINT_home"], "")
    return running


def run(PARMS):
    rs = [j for j in NUCLEOTIDES if j.startswith("R")]
    for j in rs:
        NUCLEOTIDES.append(j + "5")
        NUCLEOTIDES.append(j + "3")

    if not PARMS["only_analysis"]:
        # Trajectories
        if PARMS["Mode"] == "Traj":
            print_parameters(PARMS)
            PARMS["OUT_FILE"].write("Running with parameters: " +
                                    print_parameters(PARMS) + "\n")
            trajectories(PARMS)
            print "Description written to the file " +\
                PARMS["out_file_name_general"]

        # Single Frame
        if PARMS["Mode"] == "Single":
            PARMS["OUT_FILE"].write("Running with parameters: " +
                                    print_parameters(PARMS) + "\n")
            print "Running with parameters: \n" + print_parameters(PARMS)
            single_frame_analysis(PARMS)
            print "Output written to the " + PARMS["out_file_name_general"]

    else:
        print_parameters(PARMS)
        PARMS["OUT_FILE"].write("Running with parameters - ONLY ANALYSIS " +
                                print_parameters(PARMS) + "\n")
        trajectories(PARMS)


def reversed_h_bond(angle, distance):
    angle_cosinus = math.cos(math.radians(angle))
    x_square = distance**2+1-2*distance*angle_cosinus
    return math.sqrt(x_square)


def read_list_of_modified_nucs(fname):
    f = open(fname, "r")
    t = f.readlines()
    f.close()
    d = {}
    for i in t:
        tmp = i.split()
        d[tmp[0]] = tmp[1].strip().upper()
    return d


def inside_read_in_parms():
    PARMS = ParamReader.read_in_parameters_from_file(sys.argv[1])
    PARMS["MINT_home"] = os.path.dirname(os.path.realpath(__file__))
    PARMS["table_nucleotides"] = PARMS["MINT_home"] + "/" + \
        PARMS["table_nucleotides"]
    PARMS["table_charges"] = PARMS["MINT_home"] + "/" + \
        PARMS["table_charges"]
    PARMS["dicti_modified_nucs"] = read_list_of_modified_nucs(PARMS["MINT_home"] +
                                                              "/" +
                                                              PARMS["list_of_modified_nucs"])
    # reading in table of nucs
    PARMS["table_nuc_read"] = read_from_csv(PARMS["table_nucleotides"])
    return PARMS


def main():
    PARMS = inside_read_in_parms()
    # updadeating list of nuc
    tmp = [i[0] for i in PARMS["table_nuc_read"]][1:]
    NUCLEOTIDES.extend(tmp)

    # known modified nucleotides
    UNKNOWN_NUC.extend(PARMS["dicti_modified_nucs"].keys())

    if PARMS["Mode"] in ["Single", "Traj"]:
        PARMS["out_file_name_general"] = PARMS["out_name"] +\
            "_MINT_description.txt"
        PARMS["OUT_FILE"] = open(PARMS["out_file_name_general"], "w")
        run(PARMS)
        PARMS["OUT_FILE"].close()

    elif PARMS["Mode"] == "Download":
        for pdb_id in PARMS["pdb_list"]:
            if "#" not in pdb_id and len(pdb_id.split()) > 0:
                t = Download.run(pdb_id.split()[0], PARMS, "./  ")
                PARMS["Mode"] = "Single"
                PARMS["out_name"] = "./" + pdb_id.split()[0] + "/" + \
                    pdb_id.split()[0]+"_mint_analysis"
                PARMS["file_name"] = t
                PARMS["chains_names"] = []
                PARMS["OUT_FILE"] = open(PARMS["out_name"] + "_description",
                                         "w")
                PARMS["out_file_name_general"] = PARMS["out_name"] +\
                    "_description"
                run(PARMS)
                PARMS["OUT_FILE"].close()

if __name__ == "__main__":
    main()
    print "\nDone!"
