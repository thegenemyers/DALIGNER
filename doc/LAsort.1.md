% LAsort(1) 1.0
%
% August 2015

# NAME

LAsort - sort .las alignment files

# SYNOPSIS

**LAsort** [**-v**] *align:las* ...

# DESCRIPTION

Sort each .las alignment file specified on the command line. For each file
it reads in all the overlaps in the file and sorts them in lexicographical
order of (a,b,o,ab) assuming each alignment is recorded as
a[ab,ae] x b^o[bb,be]. It then writes them all to a file named *align.S.las*
(assuming that the input file was *align.las*). With the **-v** option set,
the program reports the number of records read and written.

## SEE ALSO

**daligner**(1)
**LAmerge**(1)
**LAshow**(1)
**LAcat**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
