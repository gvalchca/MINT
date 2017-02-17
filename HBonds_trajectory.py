# author: Anna Gorska
# !/usr/bin/python
import sys
import os
import pickle
import multiprocessing
import MDAnalysis
import numpy as np
import MINT as MINT
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/bin")


def for_a_sub_traj(nucleotides, charges, PARMS, TimeTable, ran, que, name):
    out = {}
    ppkl = name + "_" + str(min(ran)) + "_" + str(max(ran)) + ".pkl"
    if os.path.isfile(ppkl) and PARMS["only_analysis"]:
        que.append(ppkl)
        print "Not running for ", ppkl
    else:
        print "    running for ", ppkl
        for N in ran:
            dd = MINT.measure_for_all(nucleotides, charges, PARMS,
                                      TimeTable, N)
            out[N] = sum_of_hbonds(dd)
            out[N].extend(sum_of_stacking(dd))
        pickle.dump(out, open(ppkl, "wb"))
        que.append(ppkl)


def compute_rmsd(TimeTable, N, nucleotides):
    rrmsd = 0
    print "Frame ", N
    tmp = []
    for n in nucleotides:
        res = 0
        for a in n:
            if "H" not in a.get_name():
                mm = (TimeTable[a.get_full_id()][N] -
                      TimeTable[a.get_full_id()][0])
                rrmsd += np.linalg.norm(mm)
                res += np.linalg.norm(mm)
        tmp.append(round(res, 2))
    # print sum(tmp), tmp
    # print round(rrmsd, 2)
    return rrmsd


def run_for_a_trajectory(PARMS, filename, nucleotides, charges):
    universe = MDAnalysis.Universe(PARMS["file_name"], filename)
    TimeTable = MINT.ReadInTrajectory(nucleotides, universe, PARMS)
    PARMS["last_frame"] = len(universe.trajectory)
    manager = multiprocessing.Manager()
    trajs, num_of_frames, should = MINT.divide_trajectory(PARMS)
    print "Trajectory ", filename, "has ", len(universe.trajectory),
    name = filename.replace(".dcd", "")
    print " frames, running for ", num_of_frames
    que = manager.list()
    ths = [multiprocessing.Process(target=for_a_sub_traj,
                                   args=(nucleotides, charges,
                                         PARMS, TimeTable, trajs[i],
                                         que, name))
           for i in range(PARMS["threads"])]
    for p in ths:
        p.start()
    for p in ths:
        p.join()
    return que


def sum_of_hbonds(measure_for_all_dictionary):
    out = [0, 0]
    for i, kk in enumerate(["WCWC_h_bonds", "non_spec_h_bonds"]):
        for j in measure_for_all_dictionary[kk]:
            out[i] += len(j)
    ss = sum(out)
    out.append(ss)
    return out


def sum_of_stacking(measure_for_all_dictionary):
    out = [0, 0, 0]
    for j in range(3):
        out[j] = sum([nuc[j]
                      for nuc in
                      measure_for_all_dictionary["stacking_energies"]])
    return out


def put_together(pickles, fname, kk, working_dir):
    f = open(fname, "w")
    for k in kk:
        print "reading ", k
        nn = working_dir + k.replace(".dcd", "_")
        ll = [i.replace(nn, "") for i in pickles[k]]
        ll = sorted(ll, key=lambda el: int(el.split("_")[0]))
        for pkl in ll:
            oout = pickle.load(open(nn+pkl, "rb"))
            print "    reading ", nn+pkl
            oout_keys = oout.keys()
            oout_keys.sort()
            for frame in oout_keys:
                tmp = [str(i) for i in oout[frame]]
                f.write(str(frame)+","+",".join(tmp) + "\n")
    f.close()


def read_in_mint_pickles(ddir, fname):
    out = {}
    ffiles = [i for i in os.listdir(ddir) if "out_dictionary_pickle_" in i]
    for n in ffiles:
        LIST = pickle.load(open(ddir+"/"+n, "r"))
        kk = LIST.keys()
        print "Running for ", n,
        for N in kk:
            tmp = [N]
            tmp.extend(sum_of_hbonds(LIST[N][0]))
            tmp.extend(sum_of_stacking(LIST[N][0]))
            out[N] = tmp
        print " .. . done with", min(kk), "->", max(kk)
    f = open(fname, "w")
    kk = out.keys()
    kk.sort()
    print "Writing to file"
    for k in kk:
        f.write(",".join([str(i) for i in out[k]]) + "\n")
    f.close()
    print "Written to the ", fname


def run():
    PARMS = MINT.inside_read_in_parms()
    PARMS["OUT_FILE"] = open(PARMS["working_dir"]+"/"+PARMS["out_name"]+"_hbonds_log.txt", "w")
    nucleotides = MINT.get_nucleic_from_pdb(PARMS)
    charges = MINT.read_in_charges(nucleotides, PARMS)
    pickles = {}
    nuc_nums = []

    if PARMS["nucleotides"] != "":
        for i in PARMS["nucleotides"].split(";"):
            if i:
                tmp = [a for a in i.replace("(", "").replace(")", "").split('-')
                       if a]
                nuc_nums.extend(range(int(tmp[0]), int(tmp[1]) + 1))
        nuc = []
        for i in nucleotides:
            if i.get_id()[1] in nuc_nums:
                nuc.append(i)
        nucleotides = nuc
    if "out_dictionaries_MINT" in PARMS["files_dcd"][0]:
        read_in_mint_pickles(PARMS["working_dir"]+"/"+PARMS["files_dcd"][0],
                             PARMS["working_dir"]+"/"+PARMS["out_name"]+".csv")
    else:
        for filename in PARMS["files_dcd"]:
            filename_and_dir = PARMS["working_dir"] + filename
            pickles[filename] = run_for_a_trajectory(PARMS, filename_and_dir,
                                                     nucleotides, charges)
        put_together(pickles, PARMS["working_dir"] + PARMS["out_name"]+".csv",
                     PARMS["files_dcd"], PARMS["working_dir"])
        print "Written to file", PARMS["working_dir"] + PARMS["out_name"] + ".csv"


run()
