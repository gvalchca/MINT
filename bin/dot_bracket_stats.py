#!/usr/bin/python
import sys

f = open(sys.argv[1], "r")
seq = f.read()
f.close()

lines = seq.split("\n")
strs = []

for line in lines:
    if "<structure>" in line:
        strs.append(line.replace(">", "").replace("<", "").
                    replace("structure", "").replace("/", "").
                    replace("\n", "").replace(" ", ""))

d = {}
for line in strs:
    try:
        d[line] += 1
    except Exception:
        d[line] = 1

seq = seq.split("<freeSequence>")[1].split("</freeSequence>")[0]
l = len(seq)
frames = float(len(strs))

for j in sorted(d, key=d.get):
    ik = (float(d[j]) / frames) * 100.0
    if ik > 0.15:
        print j, "->", str(ik)
f = open("/workspace2/all-joanna-2/cytog_300ns/cytog-300ns_pairs_in_time.csv",
         "r")
avg = ["."] * l

for i in f.readlines()[1:]:
    tmp = i.split(",")
    if tmp[2] == "WC/WC" and float(tmp[5]) > 0.5:
        print tmp[:6]
        ik = tmp[4].replace("resid", " ").split()
        avg[int(ik[0])-1] = "("
        avg[int(ik[1])-1] = ")"
print seq
print "".join(avg)
f.close()
