% LAcat(1) 1.0
%
% August 2015

# NAME

LAcat - concatenate .las files

# SYNOPSIS

**LAcat** *source:las* > *target.las*

# DESCRIPTION

Given argument *source*, find all files *source*.1.las, *source*.2.las, ...
*source*.n.las where *source*.i.las exists for every i in [1,n]. Then
concatenate these files in order into a single .las file and pipe the result
to the standard output.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAmerge**(1)
**LAshow**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
