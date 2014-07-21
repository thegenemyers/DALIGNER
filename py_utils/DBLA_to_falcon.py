from DAPI import *
from ctypes import *
import LAPI
import sys


rcmap = dict(zip("ACGTacgtNn-","TGCATGCANN-"))
def rc(seq):
    return "".join([rcmap[c] for c in seq[::-1]])

ovl_data = LAPI.get_ovl_data(sys.argv[1])

db = HITS_DB()
open_DB(sys.argv[2], db)
trim_DB(db)
aln = LAPI.Alignment()
aln.aseq = LAPI.new_read_buffer(db)
aln.bseq = LAPI.new_read_buffer(db)

count = 0
for aread in ovl_data:
    LAPI.load_read(db, aread, aln.aseq, 2)
    aseq = cast( aln.aseq, c_char_p)
    aseq = aseq.value
    print "%08d" % aread, aseq
    for aln_data in ovl_data[aread]:
        aread, bread, acc, abpos, aepos, alen, comp, bbpos, bepos, blen = aln_data

        LAPI.load_read(db, bread, aln.bseq, 2)

        bseq = cast(aln.bseq, c_char_p)
        bseq = bseq.value
        bseq = bseq[bbpos:bepos]
        #load_read(db, ovl.bread, aln.bseq, 2)
        if comp == 1:
            bseq = rc(bseq)
        print bread, bseq
    print "+ +"
    count += 1
print "- -"

close_DB(db)
