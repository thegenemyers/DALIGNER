% LAsplit(1) 1.0
%
% August 2015

# NAME

LAsplit - divide an .las alignment file

# SYNOPSIS

**LAsplit** *target:las* {*parts:int* | *path:db|dam*} < *source.las*

# DESCRIPTION

If the second argument is an integer n, then divide the alignment file
*source*, piped in through the standard input, as evenly as possible into n
alignment files with the name *target.i.las* for i in [1,n], subject to the
restriction that all alignment records for a given a-read are in the same file.

If the second argument refers to a database *path.db* that has been
partitioned, then divide the input alignment file into block .las files where
all records whose a-read is in *path.i.db* are in *align.i.las*.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAmerge**(1)
**LAshow**(1)
**LAcat**(1)
**LAcheck**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
