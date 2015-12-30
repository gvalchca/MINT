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
# import os
from xlwt import easyxf, Workbook

HEADERS = {"motifs_clusters": ["cluster number", "motif number",
                               "motif topology", "cluster description", "VMD",
                               "% of frames", "frame numbers"],
           "average_motifs": ["cluster number", "topology",
                              "nucleotides", "VMD",
                              "% of frames", "frame numbers"],
           "motifs": ["motif topology", "nucleotides", "VMD",
                      "% of frames", "frame numbers"],
           "helices": ["nucleotides", "VMD", "% of frames",
                       "frame numbers"],
           "pseudoknots": ["nucleotides", "VMD", "% of frames",
                           "frame numbers"],
           "triplets": ["nucleotides", "VMD", "% of frames",
                        "frame numbers"],
           "stacking": ["1st nucleotide number", "nucleotides", "Coulomb avg",
                        "Coulomb sd", "VDW avg", "VDW sd",
                        "Stacking-total avg", "Stacking-total sd", "VMD",
                        "% of frames", "frame numbers"],
           "pi_interactions": ["1st nucleotide number", "nucleotides",
                               "distance avg", "distance sd", "Coulomb avg",
                               "Coulomb sd", "VDW avg", "Stacking-total sd",
                               "Stacking-total avg", "Sum sd", "VMD",
                               "% of frames", "frame numbers"],
           "pairing": ["nucleotide number", "nucleotide", "interactions"],
           "pairs": ["1st nucleotide number", "nucleotides",
                     "pair type", "pair configuration", "VMD",
                     "% of frames", "frame numbers"],
           "nucleotides_characteristics": ["nucleotide number", "nucleotide",
                                           "Hbonds-WC avg",
                                           "Hbonds-WC sd",
                                           "Hbonds-non-WC avg",
                                           "Hbonds-non-WC sd",
                                           "Hbonds-total avg",
                                           "Hbonds-total sd",
                                           "Stacking-Coulomb avg",
                                           "Stacking-Coulomb sd",
                                           "Stacking-VDW avg",
                                           "Stacking-VDW sd",
                                           "Stacking-total avg",
                                           "Stacking-total sd", "VMD"]}


HEADERS_sinlge = {"motifs": ["motif topology", "nucleotides", "VMD"],
                  "helices": ["nucleotides", "VMD"],
                  "stacking": ["1st nuclotide", "nucleotides", "Coulomb",
                               "VDW", "Stacking-total", "VMD"],
                  "pi_interactions": ["1st nucleotide number", "nucleotides",
                                      "distance", "Coulomb", "VDW",
                                      "Stacking-total", "VMD"],
                  "pairing": ["nucleotide number", "nucleotide",
                              "interactions"],
                  "pairs": ["1st nucleotide number", "nucleotides",
                            "pair type", "pair configuration", "VMD"],
                  "nucleotides_characteristics": ["nucleotide number",
                                                  "nucleotide", "Hbonds-WC",
                                                  "Hbonds-non-WC",
                                                  "Hbonds-total",
                                                  "stacking-Coulomb",
                                                  "stacking-VDW",
                                                  "stacking-total",
                                                  "VMD"],
                  "pseudoknots": ["nucleotides", "VMD"],
                  "triplets": ["nucleotides", "VMD"]}

fnames = {"helices": "_helices",
          "pseudoknots": "_pseudoknots",
          "triplets": "_triplets",
          "motifs": "_motifs",
          "pi_interactions": "_pi_interactions",
          "stacking": "_stacking",
          "pairs": "_pairs",
          "pairing": "_pairing",
          "motifs_clusters": "_motifs_clusters",
          "average_motifs": "_average_motifs",
          "nucleotides_characteristics": "_nucleotides_characteristics"}

LEGEND = {"motifs_clusters": "Secondary structure motifs clustered " +
                             "accordingly to given parameters.",
          "average_motifs": "List of average motifs, derived from the " +
                             "cluster list. Every average motif is  " +
                             "described by the vectors of the average " +
                             "numbers of hydrogen bonds per nucleotide " +
                             "forming WC and non-WC pairs.",
          "motifs": "List of all identified secondary strucuture motifs.",
          "helices": "List of all identified helices.",
          "triplets": "List of all identified triplets.",
          "pseudoknots": "List of all identified pseudoknots.",
          "stacking": "List of all nucleotide pairs recognized as stacked. " +
                      "All energies are expressed in kcal/mol.",
          "pi_interactions": "Analogous to the stacking, list of all " +
                             "ion-pi interactions. Distance is " +
                             "represented in Angstroms, energy in kcal/mol.",
          "pairing": "List of all of the hydrogen bonding partners for " +
                      "every nucleotide.The partners are enlisted " +
                      "vertically and coded: number of hydrogen bonds - " +
                      "nucleotide (name:number) - number of frames.",
          "pairs": "List of all hydrogen bonded nucleobase pairs that" +
          " appeared during the trajectory.",
          "nucleotides_characteristics": "List of all nucleotides described by" +
          "the number of hydrogen bonds in WC and " +
          "non-WC pairs and by the stacking energy broken " +
          "down into two energetic therms: " +
          "VDW and Coulomb.",
          "| notice |": "Only non-empty sheets are created."}


style_percent = easyxf('align: vert centre, horiz center;',
                       num_format_str='0.00%')
style_float = easyxf('align: vert centre, horiz center;',
                     num_format_str='##0.00')
style_index = easyxf('align: vert centre, horiz center;',
                     num_format_str='0')
style_str = easyxf('align: vert centre, horiz center;')
style_header = easyxf('font: bold on; align: vert centre, '
                      'horiz center;')


def create_csvs(LISTS, PARMS):
    if PARMS["Mode"] != "Traj":
        headers = HEADERS_sinlge
        for k in fnames.keys():
            fnames[k] = fnames[k].replace("_in_time", "")
    else:
        headers = HEADERS
    KK = LISTS.keys()
    KK.sort()

    for K in KK:
        if LISTS[K]:
            f = open(PARMS["out_name"]+fnames[K]+".csv", "w")
            f.write(",".join(headers[K]) + "\n")
            for ll in LISTS[K]:
                tmp = [str(i) for i in ll]
                f.write(",".join(tmp) + "\n")
            f.close()
            print "Written to file "+PARMS["out_name"]+fnames[K]+".csv"


def list_to_xls(LISTS, OUTNAME, PARMS):
   # try:
    list_to_xls_(LISTS, OUTNAME, PARMS)
   # except Exception:
    #    print "Your trajectory is too big to fit into the xls file",
     #   PARMS["create_csvs"] = True
      #  print ", I will create csvs."


def list_to_xls_(LISTS, OUTNAME, PARMS):
    # Create Excel workbook
    wb = Workbook(encoding='utf-8')
    KK = LISTS.keys()
    KK.sort()
    single_Legend_delete = []
    if PARMS["Mode"] != "Traj":
        headers = HEADERS_sinlge
        single_Legend_delete = ["motifs_clusters",
                                "average_motifs"]
        for k in fnames.keys():
            fnames[k] = fnames[k].replace("_in_time", "")
        for klucz in single_Legend_delete:
            if klucz in LEGEND.keys():
                del LEGEND[klucz]
    else:
        headers = HEADERS

    sheet = wb.add_sheet("LEGEND")
    kk = LEGEND.keys()
    kk.sort()
    sheet.write(0, 0, "sheet name", style_header)
    sheet.write(0, 1, "content", style_header)
    for i, k in enumerate(kk):
        i += 1
        sheet.write(i, 0, k, style_str)
        sheet.write(i, 1, LEGEND[k], style_str)

    for K in KK:
        if LISTS[K]:
            sheet = wb.add_sheet(K)
            for j, header in enumerate(headers[K]):
                if j >= 256:
                        print "too many columns ..," +\
                            " printing only 256 first ones ",
                        PARMS["create_csvs"] = True
                        continue
                else:
                    sheet.write(0, j, header, style_header)

            for i, ll in enumerate(LISTS[K]):
                i += 1
                if 'WC' in ll or 'non-WC' in ll:
                    tmp = [" "]
                    tmp.extend(ll)
                    ll = tmp
                for j, element in enumerate(ll):
                    if len(str(element)) > 32000:
                        sheet.write(i, j, "Too many frames to enumerate " +
                                    "in the xls file type." +
                                    "Please see the csv file!",
                                    style_header)
                        continue
                    if j >= 256:
                        print "too many columns .., " +\
                            "printing only 256 first ones ",
                        PARMS["create_csvs"] = True
                        continue
                    if K == "pairing" and j > 2:
                        column = headers[K][2]
                        try:
                            sheet.write(0, j, column, style_header)
                        except Exception:
                            pass
                    else:
                        column = headers[K][j]
                    if column == "% of frames":
                        sheet.write(i, j, element, style_percent)
                    elif isinstance(element, float):
                        sheet.write(i, j, element, style_float)
                    elif isinstance(element, int):
                        sheet.write(i, j, element, style_index)
                    else:
                        sheet.write(i, j, element, style_str)
    print "Saving to file: " + OUTNAME + "_MINT.xls"
    wb.save(OUTNAME + "_MINT.xls")
