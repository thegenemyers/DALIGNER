% LAmerge(1) 1.0
%
% August 2015

# NAME

LAmerge - merge .las files into a single sorted file

# SYNOPSIS

**LAmerge** [**-v**] *merge:las* *parts:las* ...

# DESCRIPTION

Merge the .las files *parts* into a singled sorted file *merge*, where it is
assumed that the input *parts* files are sorted. Due to operating system
limits, the number of *parts* files must be <= 252. With the **-v** option
set, the program reports the number of records read and written.

Used correctly, **LAmerge** and **LAsort**(1) together allow one to perform
an "external" sort that produces a collection of sorted files containing in
aggregate all the local alignments found by the **daligner**(1), such that
their concatenation is sorted in order of (a,b,o,ab). In particular, this
means that all the alignments for a given a-read will be found consecutively
in one of the files. So computations that need to look at all the alignments
for a given read can operate in simple sequential scans of these sorted files.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAshow**(1)
**LAcat**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
