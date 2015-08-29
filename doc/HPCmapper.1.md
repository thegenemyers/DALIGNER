% HPCmapper(1) 1.0
%
% August 2015

# NAME

HPCmapper - generate a script to map reads

# SYNOPSIS

**HPCmapper** [**-vb**] [**-k***int(20)*] [**-w***int(6)*] [**-h***int(50)*]
	[**-t***int*] [**-M***int*] [**-e***double(.85)*]
	[**-l***int(1000)*] [**-s***int(100)*] [**-H***int*]
    [**-m***track*]+ [**-dal***int(4)*] [**-deg***int(25)*]
	*ref:db|dam* *reads:db|dam* [*first:int*[-*last:int*]]

# DESCRIPTION

**HPCmapper** writes a UNIX shell script to the standard output that
consists of a sequence of commands that effectively "maps" every read in
the DB *reads* against a reference set of sequences in the DB *ref*,
recording all the found local alignments in the sequence of files
*ref.reads.1.las*, *ref.reads.2.las*, and so on, where *ref.reads.k.las*
contains the alignments between all of *ref* and the k'th block of
*reads*.  The parameters are exactly the same as for **HPCdaligner**(1)
save that the **-k**, **-h**, and **-e** defaults are set
appropriately for mapping, and the **-A** and **-I** options
make no sense as *ref* and *reads* are expected to be distinct
data sets.

If the integers *first* and *last* are missing, then the
script produced is for every block in the database *reads*.
If *first* is present then **HPCmapper** produces an script
that compares blocks *first* through *last* (*last* = *first*
if not present) against DAM *ref*.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAmerge**(1)
**LAshow**(1)
**LAcat**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCdaligner**(1)
