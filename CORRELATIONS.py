#!/usr/bin/python
from operator import itemgetter
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from collections import Counter
from matplotlib import rc
rc('text', usetex=True)
font = {'family': 'serif',
        'serif': ['Computer Modern'],
        'weight': 'bold',
        'size': 30}
rc('mathtext', fontset='stixsans')
rc('font', **font)


def read_in_pairs_csv(filename):
    f = open(filename, "r")
    t = f.readlines()[1:]
    f.close()
    wcwc_pairs = []
    #  computing number of frames:
    tmp = t[2].split(",")
    i = 3
    while float(tmp[5]) < 0.5:
        tmp = t[i].split(",")
        i += 1
    number_of_frames = int(float(len(list_from_shorted(tmp[6])))/float(tmp[5]))
    number_of_frames = 3000
    for i in t:
        tmp = i.split(",")
        if float(tmp[5]) >= 0.3 and tmp[2] != "WC/WC":
            a = make_int(tmp[0])
            wcwc_pairs.append([a, tmp[1], list_from_shorted(tmp[6])])
    wcwc_pairs.sort(key=itemgetter(0))
    return (wcwc_pairs, number_of_frames)


def make_int(i):
    try:
        a = int(i)
    except Exception:
        a = int(i.split(":")[1])
    return a


def list_of_nucleotides_instead_of_pairs(wcwc_pairs, ranges):
    t = {}
    nucleotides = {}
    for pair in wcwc_pairs:
        for nuc in pair[1].split("/"):
            num = make_int(nuc)
            try:
                t[num].extend(pair[2])
            except Exception:
                t[num] = pair[2]
            if num not in nucleotides.keys():
                nucleotides[num] = nuc

    # only from ranges
    if not ranges:
        keys = t.keys()
        ranges = [(min(keys), max(keys))]

    new_t = {}
    for r in ranges:
        for a in range(r[0], r[1]+1):
            try:
                new_t[a] = t[a]
            except Exception:
                new_t[a] = []

    li = [[nuc, set(new_t[nuc])] for nuc in new_t]
    li.sort(key=itemgetter(0))
    return (li, nucleotides)


def list_from_shorted(l):
    out = []
    for i in l.split():
        if "->" in i:
            tmp = i.split("->")
            out.extend(range(int(tmp[0]), int(tmp[1])+1))
        else:
            out.extend([int(i)])
    return out


def heatmap(data, indexes, colormap, fileout):
    column_labels = indexes
    row_labels = indexes

    parms = [0.05, 0.15, 0.85, 0.90]
    s = float(len(indexes))
    x = int(20 * s/20 + 1)
    y = int(8 * s/20 + 1)
    fig, ax = plt.subplots(figsize=(x, y))
    colormap = matplotlib.colors.ListedColormap(["DodgerBlue", "white", "red"],
                                                name='from_list', N=None)
    fig.subplots_adjust(left=parms[0], bottom=parms[1],
                        right=parms[2], top=parms[3])
    fig.patch.set_facecolor('white')
    fig.patch.set_edgecolor('grey')

    heatmap = ax.pcolor(data, cmap=colormap, vmin=-1, vmax=1)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    # legend
    cbar = plt.colorbar(heatmap)
    cbar.ax.get_yaxis().set_ticks([])
    for j, lab in enumerate(['$-1$', '$0$', '$1$']):
        cbar.ax.text(0.5, (2*j+1) / 6.0, lab, ha='center',
                     va='center', fontsize=18, weight='bold')
    # cbar.ax.get_yaxis().set_ticks(['-1','0','1'])
    cbar.ax.set_ylabel(r'$\phi$', rotation=270, labelpad=40, fontsize=50)

    # cutting stuff
    plt.gca().set_xlim((0, len(indexes)))
    plt.gca().set_ylim((0, len(indexes)))

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    #  plt.show()
    print "Writting to file: ", fileout.replace(".csv", ".png")
    plt.savefig(fileout.replace(".csv",".png"))


def correlation_for_two(p1, p2, no_of_frames, cutoff, f):
    t1 = set(p1)
    t2 = set(p2)
    n11 = float(len(t1 & t2))  # when both are present
    n01 = float(len(t2.difference(t1)))  # when is pair in p2 and not in p1
    n10 = float(len(t1.difference(t2)))  # when is pair in p1 and not in pa
    # not present in the same time
    n00 = float(no_of_frames - len(t2.union(t1)))

    n1a = n11+n01
    n1b = n11+n10
    n0a = n00+n10
    n0b = n00+n01

    s = max(n1a*n1b*n0a*n0b, 1)

    # check only the changing ines
    if int(n00) == no_of_frames or int(n11) == no_of_frames:
        return 1.0
    if s > 0:
        m = ((n11*n00 - n10*n01)/math.sqrt(s))
        if m >= cutoff:
            return 1
        if m <= (-1)*cutoff:
            return -1
        else:
            return 0.0
    else:
        print n11, " ", n01, " ", n10, " ", n00
        print n1a, " ", n1b, " ", n0a, " ", n0b
        print "++++"
        return 0.0


def compute_correlation_nucleotides(nucleotides_list, no_of_frames, cutoff):
    outlist = []
    for i, p1 in enumerate(nucleotides_list):
        tmp = [0.0]*i
        for j, p2 in enumerate(nucleotides_list[i:]):
            tmp.append(correlation_for_two(p1[1], p2[1],
                                           no_of_frames,
                                           cutoff, 0))
        outlist.append(tmp)
    # outlist[-1][0] =-1
    return np.array(outlist)


def indexes_nucleotides(nucleotides_list, nucleotides):
    indexes = [i[0] for i in nucleotides_list]
    return indexes


def save_matrix_to_csv(corr, indexes, file_out):
    f = open(file_out, "w")
    for i in indexes:
        f.write(str(i) + ", ")
    f.write("\n")
    for ind, row in enumerate(corr):
        f.write(str(indexes[ind]) + ", ")
        for a in row:
            f.write(str(a) + ", ")
        f.write("\n")
    f.close()


def read_in_matrix(file):
    f = open(file, "r")
    t = f.readlines()
    f.close()
    indexes = t[0].split(",")
    out = []
    for i in t[1:]:
        tmp = i.split(",")[1:-1]
        out.append([float(j) for j in tmp])
    return (out, indexes)


def substract_two_matrixes(file1, file2, cutoff):
    tmp = read_in_matrix(file1)
    m1 = tmp[0]
    indexes1 = tmp[1][:-1]
    tmp = read_in_matrix(file2)
    m2 = tmp[0]
    indexes2 = tmp[1][:-1]

    # substracting

    if len(indexes1) != len(indexes2):
        print " Matrices differ in size !"
        sys.exit()

    substracted = []
    width = len(indexes1)

    for i, rows in enumerate(zip(m1, m2)):
        tmp = [0.0]*i
        tmp.append(round(-1.0+(i*(2/float(len(indexes1)-1))), 2))
        tmp.extend(substract_two_rows(rows, i, cutoff))
        substracted.append(tmp)

    heatmap(np.array(substracted), indexes1, plt.cm.PRGn)
    return (np.array(substracted), indexes1)


def substract_two_rows(rows, i, cutoff):
    row1 = rows[0]
    row2 = rows[1]
    out = []
    for i in range(i+1, len(row1)):
        #  out.append(round(row1[i]-row2[i] ,2))
        if ((row1[i] >= cutoff and row2[i] >= cutoff) or
           (row1[i] <= (-1)*cutoff and row2[i] <= (-1)*cutoff)):
            out.append(1)
        elif ((row1[i] <= cutoff and row1[i] >= (-1)*cutoff) or
              (row2[i] <= cutoff and row2[i] >= (-1)*cutoff)):
            out.append(0)
        elif ((row1[i] <= (-1)*cutoff and row2[i] >= cutoff) or
              (row2[i] <= (-1) * cutoff and row2[i] >= cutoff)):
            out.append(-1)
        else:
            out.append(0)
    return out


def run(file_in, ranges, file_out, cutoff):
    tmp = read_in_pairs_csv(file_in)
    wcwc_pairs = tmp[0]
    no_of_frames = tmp[1]
    tmp = list_of_nucleotides_instead_of_pairs(wcwc_pairs, ranges)
    nucleotides_list = tmp[0]
    nucleotides = tmp[1]
    print "no_of_frames ", no_of_frames
    corr = compute_correlation_nucleotides(nucleotides_list, no_of_frames,
                                           cutoff)
    indexes = indexes_nucleotides(nucleotides_list, nucleotides)
    save_matrix_to_csv(corr, indexes, file_out)
    heatmap(corr, indexes, plt.cm.seismic, file_out)


def clustering_analysis(file1):
    m1 = []
    f = open(file1, "r")
    lines = f.readlines()
    in1 = lines[0].split(",")[:-1]
    for row in lines[1:]:
        tmp = row.split(",")
        o = [tmp[0]]
        for a in tmp[1:-1]:
            o.append(float(a))
        m1.append(o)
    f.close()
    print len(m1)
    clusters = {}
    clusters[0] = [1]
    black = []
    while True:
        # co go jeszcze nie ma
        all_colored = []
        for i in clusters:
            all_colored.extend(clusters[i])
            gray = list(set(range(0,
                            len(m1[0])-1)).difference(set(all_colored)))
            if not gray:
                break
        num = gray[0]
        clusters[num] = []
        row = num+1
        while True:
            clusters[num].append(row-1)
            black.append(row-1)
            for i, val in enumerate(m1[row-1][1:]):
                if val:
                    clusters[num].append(i)
            gray = list(set(clusters[num]).difference(set(black)))
            if not gray:
                break
            row = gray[0]+1

    # translate clusters
    clust = {}
    for i in clusters:
        if clusters[i]:
            clust[i] = []
            for t in set(clusters[i]):
                clust[i].append(in1[t])
    for i in clust:
        if len(clust[i]) > 2:
            print i, " -> ", clust[i]
    all_coloured = []
    for i in clust:
        all_coloured.extend(clust[i])
    u = Counter(all_coloured)
    for i in u:
        if u[i] > 1:
            print i, " ", u[i]
    print len(all_coloured)
    print len(m1[0])


def next_nucleotide(m1, clusters):
    all_coloured = []
    for i in clusters:
        all_coloured.extend(clusters[i])
    n = set(range(len(m1[0])-1)).difference(set(all_coloured))
    rows = list(n)
    if rows:
        return rows[0]+1
    else:
        False


def main():
    try:
    # if True:
        cutoff = float(sys.argv[2])
        nucleotides = []
        if len(sys.argv) >= 4:
            for i in sys.argv[3].replace("[", "").replace("]", "").replace(" ", "").split(")"):
                if i:
                    tmp = [a for a in i.replace("(", "").split(",") if a]
                    nucleotides.append((int(tmp[0]), int(tmp[1])))

        file_in = sys.argv[1]
        file_out = file_in.replace("pairs.csv",
                                   "matrix_"+str(cutoff)+".csv")
        print """    running for file """ + file_in + """
            with cutoff """ + str(cutoff) + """
            for nucleotides """, nucleotides, """
            producing file """ + file_out

        run(file_in, nucleotides, file_out, cutoff)
    except Exception:
        print """
        specify the input file (pairs.csv),
        cutoff and if you whish the numbers of nucleotides to
        compute the correlation matrix example:
        python CORRELATIONS.py step6_pairs_in_time.csv 0.4 \"[(1210,1220), (985,995) , (1043,1047)]\"
        Do not forget \" around your nucleotide list!
        For more information please see MINT tutorial."""
        sys.exit()


main()
