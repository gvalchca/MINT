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


def RNAStructMl(file_out, l, seq, strucuture_id):
    f = open(file_out, "w")
    f.write("""<?xml version="1.0" encoding="utf-8"?>
            <rnastructML xmlns="http://hobit.sourceforge.net/xsds/20060201/rnastructML">
            <rnastructure id=  " """ + strucuture_id +
            """ "><sequence seqID="Sequence_1"> """ +
            """ <freeSequence>""" + seq + """</freeSequence>  </sequence> \n""")
    for i in l:
        f.write(" <structure>" + l[i] + "</structure>\n")
    f.write("""</rnastructure> </rnastructML>""")
    f.close()
    print "Written dotBracket to " + file_out
