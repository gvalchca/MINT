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
import collections
import time
import math


def timer(f):
    def g(*args, **kwargs):
        start = time.clock()
        result = f(*args, **kwargs)
        open("../time.txt", 'a').write(f.__name__ + ": " +
                                       str(time.clock() - start) + "\n")
        return result
    return g


def print_one_dimension_list(l, begining):
    counter = 0
    arches = ""
    width = len(str(len(l)))
    ruler = ""
    seq = ""
    for i in l:
        seq += " " + str(counter + begining) +\
            " "*(width - len(str(counter + begining)))
        ruler += " " + str(counter) + " "*(width - len(str(counter)))
        counter += 1
        if isinstance(i, int):
            arches += " "+str(i)+" "*(width-len(str(i)))
        # if N
        else:
            arches += " " + "_" * width
    print seq
    print ruler
    print arches


# @timer
def create_list(WCWC_nucs, WCWC_clas, nucleotides):
    arches = [None]*len(nucleotides)
    WCWC_triplex = []
    for i in range(0, len(WCWC_nucs)):
        tmp = WCWC_clas[i].split("/")
        tmp_1 = WCWC_nucs[i]
        beg = tmp_1[0]
        end = tmp_1[1]
        # classical secondary structure - Watson Crick
        if (tmp[0] == tmp[1] == "WC") and (int(tmp[2]) >= 1):
            if end not in arches and beg not in arches:
                arches[end] = beg
                arches[beg] = end
            elif end in arches:
                WCWC_triplex.append([beg, end, arches[end]])
            elif beg in arches:
                WCWC_triplex.append([end, beg, arches[beg]])
    return (arches, WCWC_triplex)


# Helper
def get_all_indexes(l, element):
    return [item for item in range(len(l)) if l[item] == element]


#@timer
def triplexes(WCWC_nucs):
    t1 = []
    t2 = []
    for i in range(len(WCWC_nucs)):
        t1.append(WCWC_nucs[i][0])
        t2.append(WCWC_nucs[i][1])
    r = t1+t2
    cnt = collections.Counter(r)
    more_than_one_names = [k for k in cnt if cnt[k] > 1]
    out = []
    for i in more_than_one_names:
        ts = [t1, t2]
        for t in ts:
            for j in get_all_indexes(t, i):
                if t == t1:
                    other = t2
                else:
                    other = t1
                partner = other[j]
                out_tmp = [i, partner]
                 #track partner:
                if partner in more_than_one_names:
                    more_than_one_names.remove(partner)
                    m = [x for x in get_all_indexes(t, partner) if x != j]
                    if m:
                        out_tmp.append(other[m[0]])
                out.append(out_tmp)
    more_out = []
    i = 0
    while i < len(out)-1:
        if (set(out[i]) & set(out[i+1])):
            if len(out[i]) <= len(out[i+1]):
                out[i].reverse()
                tmp = out[i] + out[i+1]
            else:
                out[i+1].reverse()
                tmp = out[i+1] + out[i]
            more_out.append(unique_list(tmp))  # unique list
        else:
            more_out.extend([out[i], out[i+1]])
        i += 2
    return more_out


def unique_list(l):
    out = []
    for i in l:
        if i not in out:
            out.append(i)
    return out


def resid_description(res):
    if res:
        ch = str(res.get_parent().get_id())
        return ch + "|" + res.get_resname().strip()+":"+str(res.get_id()[1])
    else:
        return ""


def names(resname):
    resname = resname.replace(" ", "")
    if resname in ["G", "U", "A", "C"]:
        return resname
    elif resname in ["GUA", "URA", "CYT", "ADE"]:
        return resname[0]
    elif resname.startswith("R"):
        return resname[1]
    elif "U" in resname:
        return "U"
    elif "A" in resname:
        return "A"
    elif "G" in resname:
        return "G"
    elif "C" in resname:
        return "C"
    elif "DT" in resname:
        return "T"
    return resname[0]


def motif_numerical_description(d):
    s = d.split("P")
    t = [str(len(i)) for i in s[1:]]
    return "-".join(t)


def motif_description(m, l, nucleotides, if_print):
    out = "\nMotifs:\n"
    for i in range(len(m)):
        t = ""
        for j in m[i][1]:
            t = t + resid_description(nucleotides[j]) + "-"
        out += str(i+1) + "]  " + m[i][0] + "  " + t + "\n"
    return (out)


def pseudoknots_description(p, l, nucleotides, if_print):
    out = "\nPseudoknots:\n"
    for i in range(0, len(p)):
        for j in range(0, len(p[i])):
            p1 = resid_description(nucleotides[p[i][j][0]])
            p2 = resid_description(nucleotides[p[i][j][1]])
            out += p1 + "-" + p2 + " "
    return (out)


def helices_description(h, l, nucleotides, if_print):
    out = "\nHelices: \n"
    TMP = []
    for i in range(0, len(h)):
        l_beg = h[i][0]
        l_end = h[i][len(h[i])-1]
        if l_beg < l[l_beg]:
            tmp = h[i]
        else:
            tmp = [i for i in range(l[l_end], l[l_beg])]
        if tmp and not list((collections.Counter(tmp) &
                             collections.Counter(TMP)).elements()):
            beg = resid_description(nucleotides[l_beg])
            end = resid_description(nucleotides[l_end])
            beg_sec = resid_description(nucleotides[l[l_end]])
            end_sec = resid_description(nucleotides[l[l_beg]])
            if len(h[i]) > 1:
                TMP.extend(tmp)
                out += str(i + 1) + "] " + beg + "-" + end + " -> " +\
                    beg_sec + "-" + end_sec + "\n"
            else:
                TMP.extend(tmp)
                out += str(i + 1) + "] " + beg + " -> " + end_sec + "\n"
    return (out)


def number_from_resid(n):
    return n.split(":")[1]


def triplex_despription(t, nucleotides, if_print):
    out = "\nTriplexes: \n"
    for i in range(len(t)):
        triplex_tmp = ""
        for j in t[i]:
            triplex_tmp += resid_description(nucleotides[j]) + "-"
        out += str(i+1)+"]  " + triplex_tmp + "\n"
    return (out)


# @timer
def pseudo_knots_finder(l):
    ps = []
    ps_2 = []
    ps_1 = []
    for i in range(len(l)):
        if l[i] is None:
            ps_1.append(i)
            for j in range(i+1, len(l)):
                if l[j] is None:
                    if ((l[i] < l[j]) and (j < l[i])):
                        ps_2.append(j)
            if ps_1 and len(ps_1) < len(ps_2):
                ps.append(erase_ps(ps_1, l))
            elif ps_2:
                ps.append(erase_ps(ps_2, l))
            ps_1 = []
            ps_2 = []
    return ps


def erase_ps(l_tmp, l):
    ps = []
    for i in l_tmp:
        tmp = l[i]
        ps.append([i, tmp])
        l[i] = None
        l[tmp] = None
    return ps


#@timer
def structural_analizer(l):
    i = 0
    h = []
    m = []
    #  przeskocz niesparowany poczatek
    while l[i] is None and i < len(l)-1:
        i += 1
    while i < len(l)-1:
        helix_tmp = []
        while isinstance(l[i], int) and i < len(l)-1:
            helix_tmp.append(i)
            if l[i+1] is not None and l[i+1]+1 == l[i] and l[i] > i:
                i += 1
            else:
                h.append(helix_tmp)
                break
        #motifs finder
        if l[i] > i:
            motif = [i, l[i]]
            i = l[i]-1
            last_pair = i
            motif_class = "P"
            while i != motif[0]:
                if l[i] is None:
                    motif_class += "U"
                    motif.append(i)
                    i -= 1
                else:
                    motif_class += "P"
                    motif += [i, l[i]]
                    i = l[i]
                    last_pair = i
                    i -= 1
            motif.append(i)
            i = last_pair-1
            m.append([motif_numerical_description(motif_class), motif])
            i += 1
        else:
            i += 1
    return (h, m)


def arches_to_dot_bracket(arches):
    out = ["."] * len(arches)
    for i in range(len(arches)):
        if arches[i] is not None:
            if i < arches[i]:
                out[i] = "("
            elif i > arches[i]:
                out[i] = ")"
    return "".join(out)


def print_megalist_nicely(megalist, typ, nucleotides):
    out = ""
    if typ == "WCWC":
        WCWC_nucs = megalist["WCWC_nucs"]
        WCWC_h_bonds = megalist["WCWC_h_bonds"]
        WCWC_clas = megalist["WCWC_clas"]
        WCWC_conf = megalist["WCWC_conf"]
    else:
        WCWC_nucs = megalist["non_spec_nucs"]
        WCWC_h_bonds = megalist["non_spec_h_bonds"]
        WCWC_clas = megalist["non_spec_clas"]
        WCWC_conf = megalist["non_spec_spec"]
    for i in range(len(WCWC_nucs)):
        out += str(i) + "] " +\
            resid_description(nucleotides[WCWC_nucs[i][0]]) +\
            "/" + resid_description(nucleotides[WCWC_nucs[i][1]]) + "  " +\
            print_hbonds_nicely(WCWC_h_bonds[i]) + WCWC_clas[i] + "  " +\
            WCWC_conf[i] + "\n"
    return out


def reversed_h_bond(angle, distance):
        angle_cosinus = math.cos(math.radians(angle))
        x_square = distance**2+1-2*distance*angle_cosinus
        return math.sqrt(x_square)


def print_hbonds_nicely(h_bonds_list):
    tmp = []
    for i in h_bonds_list:
        x = reversed_h_bond(i[3][0], i[3][1])
        tmp.append(i[0] + "-" + i[1] + "-" + i[2] + " angle: " +
                   str(round(i[3][0], 2)) +
                   " distance: " + str(round(x, 2)) + " ")
    return " | ".join(tmp)


def get_sequence(nucleotides, IF_MODIFIED):
    sequence = ""
    for residue in nucleotides:
        if not IF_MODIFIED[residue.get_id()]:
            sequence += names(residue.get_resname())
        else:
            sequence += IF_MODIFIED[residue.get_id()].lower()
    return sequence


# @timer
def analysis(megalist, nucleotides, if_print, IF_MODIFIED):
    f = ""
    tmp_arches = create_list(megalist["WCWC_nucs"],
                             megalist["WCWC_clas"],
                             nucleotides)
    arches = tmp_arches[0]
    ps = pseudo_knots_finder(arches)
    tmp = structural_analizer(arches)
    t = triplexes(megalist["WCWC_nucs"])
    trajectory = []
    # (h,m)
    if tmp[0]:
        heli_desc = helices_description(tmp[0], arches, nucleotides, if_print)
        f += heli_desc[0] + heli_desc[1]
        trajectory.append(heli_desc[0])
    else:
        trajectory.append("")
    if ps:
        pseudo_desc = pseudoknots_description(ps, arches, nucleotides,
                                              if_print)
        f += (pseudo_desc[0] + pseudo_desc[1])
        trajectory.append(pseudo_desc[0])
    else:
        trajectory.append("")
    if t:
        triplex_desc = triplex_despription(t, nucleotides, if_print)
        f += triplex_desc[0] + triplex_desc[1]
        trajectory.append(triplex_desc[0])
    else:
        trajectory.append("")
    if tmp[1]:
        motif_desc = motif_description(tmp[1], arches, nucleotides, if_print)
        trajectory.append(motif_desc[0])
        f += motif_desc[0] + motif_desc[1]
    else:
        trajectory.append("")

    seq = get_sequence(nucleotides, IF_MODIFIED)
    f += "\nWC Base pairs: \n" +\
        print_megalist_nicely(megalist, "WCWC", nucleotides) +\
        "\nOther interactions: \n" +\
        print_megalist_nicely(megalist, "not-spec", nucleotides)
    f += "\nDot-Bracket\n" + seq+"\n" +\
        arches_to_dot_bracket(arches) +\
        "\n Modified nucleotides denoted by lower case letters\n"

    return [arches, ps, t, tmp[0], tmp[1], 
            trajectory, arches_to_dot_bracket(arches), f]
