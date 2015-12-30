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
#*   This program is distributed in the hope that it whileill be useful,    *
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
    motif_list = f.readlines()[1:]
    f.close()
    mots = [i.replace("\n", "").split(",") for i in motif_list]
    return mots


def create_table(l):
    t = []
    for i in range(len(l)):
        r = []
        a = l[i][1].split("-")[:-2]
        for j in range(len(l)):
            b = l[j][1].split("-")[:-2]
            r.append(intersect(a, b))
        t.append(r)
    return t


def cluster_from_table(t, l, margin):
    friends = {}
    for i in range(len(l)):
        friends[i] = [j for j in range(len(l)) if t[i][j] >=
                      margin * min(len(l[i][1].split("-")[:-2]),
                      len(l[j][1].split("-")[:-2]))]
    clusters = {}

    while friends:
        keys = friends.keys()
        keys.sort(key=lambda x: len(friends[x]), reverse=True)
        used = friends.pop(keys[0])
        clusters[keys[0]] = used
        for i in friends.keys():
            if i in used:
                friends.pop(i)
            else:
                tmp = remove_from_list(friends[i], used)
                if tmp:
                    friends[i] = tmp
    return clusters


def print_clusters_nicely(clusters, l):
    out = []
    keys = clusters.keys()
    clusters_maximas = []
    for c in range(len(keys)):
        cls = "cluster "+str(c)
        maximal = clusters[keys[c]][0]
        for i in clusters[keys[c]]:
            ttmp = [cls, i]
            ttmp.extend(l[i])
            out.append(ttmp)
            if len(l[i][1]) > len(l[maximal][1]):
                maximal = i
        clusters_maximas.append(maximal)
    return (out, clusters_maximas)


def remove_from_list(l, used):
    tmp = [i for i in l if i not in used]
    return tmp


def print_table(t):
    out = " ," + ",".join(int_list_to_string(range(len(t))))+"\n"
    for i in range(len(t)):
        out += str(i) + "," + ",".join(int_list_to_string(t[i]))+"\n"
    return out


def int_list_to_string(l):
    m = [str(i) for i in l]
    return m


def list_to_dictionary(mots):
    d = {}
    for i in mots:
        d[i[1]] = [[i[0]]+i[2:]]
    return d


def print_ditcti(d):
    for k in d:
        print k, d[k]


def intersect(a, b):
    return len(list(set(a) & set(b)))


def filter(l, percentage, overall):
    m = [i for i in l if float(i[3]) > percentage]
    return m


def cluster(margin, overall, l, percentage):
    l = filter(l, percentage, overall)
    t = create_table(l)
    clusters = cluster_from_table(t, l, margin)
    out = print_clusters_nicely(clusters, l)
    return [l, out[1], out[0]]
