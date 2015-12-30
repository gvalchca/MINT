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
#****************************************************************************

#!/usr/bin/python
from ftplib import FTP
import os
import platform


normal_db = "pub/pdb/data/structures/all/pdb/"
obsolete = "pub/pdb/data/structures/obsolete/pdb/"


#downloading the pdb file from the ftp of pdb database
def download(name, workingDirectory):
    ftp = FTP('ftp.wwpdb.org')
    ftp.login()
    filename = name + ".ent.gz"
    directory = 'pub/pdb/data/structures/all/pdb/'
    ftp.cwd(directory)
    if name not in os.listdir("./"):
        os.system("mkdir "+workingDirectory+"/"+name)
    ff = workingDirectory+"/" + name + "/" + filename
    f = open(ff, "wb")
    ftp.retrbinary('RETR ' + "pdb" + filename.lower(), f.write)
    f.close()
    ftp.quit()
    return ff


#function uncompressing the file
def uncompress(filename):
    os.system("gunzip " + filename)
    return filename.replace(".gz", "")


def protonate(filename, MINT_home):
    print "Protonating .. "
    ss = platform.system()
    if ss == "Linux":
        os.system(MINT_home+"/bin/reduce -build " + filename + " > " +
                  filename.replace(".ent", "_H.pdb"))
    elif ss == "Darwin":
        os.system(MINT_home+"/bin/reduce_mac -build " + filename + " > " +
                  filename.replace(".ent", "_H.pdb"))

    return filename.replace(".ent", "_H.pdb")


def run(pdb, parms, workingDirectory):
    print "\nDownloading "+pdb
    t = protonate(uncompress(download(pdb, workingDirectory)),
                  parms["MINT_home"])
    print t
    return t
