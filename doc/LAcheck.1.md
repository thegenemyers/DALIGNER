% LAcheck(1) 1.0
%
% August 2015

# NAME

LAcheck - verify structural integrity of .las files

# SYNOPSIS

**LAcheck** [**-vS**] *src1:db|dam* [*src2:db|dam*] *align:las* ...

# DESCRIPTION

LAcheck checks each .las file for structural integrity, where the a- and
b-sequences come from *src1* or from *src1* and *src2*, respectively.
That is, it makes sure each file makes sense as a plausible .las file,
e.g. values are not out of bound, the number of records is correct, the number
of trace points for a record is correct, and so on. The exit status is 0 if
every file is deemed good, and 1 if at least one of the files looks corrupted.

# OPTIONS

**-S**
:   Also check that the alignments are in sorted order

**-v**
:   Print a line for each .las file saying either the file is OK or reporting
	the first error. If the **-v** option is not set, then the program
	runs silently.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAmerge**(1)
**LAshow**(1)
**LAcat**(1)
**LAsplit**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
