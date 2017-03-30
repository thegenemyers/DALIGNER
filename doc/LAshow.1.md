% LAshow(1) 1.0
%
% August 2015

# NAME

LAshow - display local alignments from .las files

# SYNOPSIS

**LAshow** [**-caroUF**] [**-i***int(4)*] [**-w***int(100)*]
[**-b***int(10)*] *src1:db|dam* [*src2:db|dam*] *align:las*
[*reads:FILE* | *reads:range* ... ]

# DESCRIPTION

**LAshow** produces a printed listing of the local alignments contained in the
specified .las file, where the a- and b-reads come from *src1* or from *src1*
and *src2*, respectively. If a file or list of read ranges is given, then only
the overlaps for which the a-read is in the set specified by the file or list
are displayed. See **DBshow**(1) for an explanation of how the file and list
of read ranges are interpreted.

# OPTIONS

**-F**
:   Reverse the roles of the a- and b- reads in the display

**-U**
:   Use uppercase for DNA sequence instead of the default lowercase

**-i***indentation*
:   Set the indent for the cartoon and/or alignment displays if they are
	requested. (default: 4)

**-b***num_symbols*
:   Set the number of symbols on either side of the aligned segments in an
	alignment display. (default: 10)

**-w***w*
:   This parameter is used for the display modes specified
	by **-a** and **-r**. (default: 100)

**-o**
:   Only display alignments that are proper overlaps-- that is, where a
	sequence end occurs at each end of the alignment

## DISPLAY MODES

**-c**
:   Cartoon rendering of the alignment

**-a**, **-r**
:   Display an alignment of the local alignment

The **-a** option puts exactly *w* columns per segment of
the display, whereas the **-r** option puts exactly *w* a-read symbols in each
segment of the display. The **-r** display mode is useful when one wants to
visually compare two alignments involving the same a-read. If a combination of
the **-c**, **-a**, and **-r** flags is set, then the cartoon comes first,
then the **-a** alignment, and lastly the **-r** alignment.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAmerge**(1)
**LAcat**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
