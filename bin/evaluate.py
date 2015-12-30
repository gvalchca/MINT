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

#!/usr/bin/python


def read_in_file(file_name):
    f = open(file_name)
    nuc_list = f.readlines()
    f.close()
    nucs = [i.replace("\n", "").split(";") for i in nuc_list]
    return nucs


def number_per_nuc(nuc, global_number_of_frames):
    sum = 0.0
    for i in nuc[2:]:
        tmp = i.split(":")
        number_of_frames = float(tmp[1])
        number_of_hydrogen_bonds = float(tmp[0].split("-")[1])
        sum += (float(number_of_frames) /
                float(global_number_of_frames) *
                number_of_hydrogen_bonds)
    return sum


def lista_number_per_nuc(l, global_number_of_frames):
    out = {}
    for nuc in l:
        out[nuc[0]] = number_per_nuc(nuc, global_number_of_frames)
    return out


def output_to_csv(t, file_name):
    f = open(file_name, "w")
    for key in t:
        f.write(key+";"+key[1:] + ";"+str(t[key])+"\n")
    f.close()
    print "Written to file " + file_name


def interscript(l, file_name, global_number_of_frames):
    t = lista_number_per_nuc(l, global_number_of_frames)
    output_to_csv(t, file_name.replace(".csv",
                                       "_nucleotide_window_evaluation.csv"))


def main():
    t = read_in_file("../dwa_tys/100ns_nucleotides_characteristics.csv")
    m = lista_number_per_nuc(t, 20000)
    output_to_csv(m, "ik")
