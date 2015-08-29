% HPCdaligner(1) 1.0
%
% August 2015

# NAME

HPCdaligner - generate a script to run **daligner**(1)

# SYNOPSIS

**HPCdaligner** [**-vbAI**] [**-k***int(14)*] [**-w***int(6)*]
	[**-h***int(35)*] [**-t***int*] [**-M***int*]
	[**-e***double(.70)*] [**-l***int(1000)*] [**-s***int(100)*] [**-H***int*]
	[**-m***track*]+ [**-dal***int(4)*] [**-deg***int(25)*]
	*path:db|dam* [*first:int*[-*last:int*]]

# DESCRIPTION

**HPCdaligner** writes a UNIX shell script to the standard output that consists
of a sequence of commands that effectively run **daligner**(1) on all pairs of
blocks of a split database and then externally sorts and merges them using
**LAsort**(1) and **LAmerge**(1) into a collection of alignment files with
names *path.#.las* where # ranges from 1 to the number of blocks the database
is split into. These sorted files if concatenated by say **LAcat**(1)
would contain all the alignments in sorted order (of a-read, then b-read, and so on).
Moreover, all overlaps for a given a-read are guaranteed to not be split across
files, so one can run artifact analyzers or error correction on each sorted
file in parallel.

The database must have been previously split by **DBsplit**(1) and all the
parameters, except **-v**, **-dal**, and **-deg**, are passed through to the
calls to **daligner**(1). The defaults for these parameters are as for
**daligner**(1). The **-v** flag, for verbose-mode, is also passed to all
calls to **LAsort**(1) and **LAmerge**(1). **-dal** and **-deg** options are
described later.

For a database divided into N sub-blocks, the calls to **daligner**(1) will
produce in total 2TN^2 .las files assuming daligner runs with T threads.
These will then be sorted and merged into N^2 sorted .las files, one for each
block pair. These are then merged in ceil(log_deg N) phases where the number of
files decreases geometrically in **-deg** until there is 1 file per row of the
N x N block matrix. So at the end one has N sorted .las files that when
concatenated would give a single large sorted overlap file.

The **-dal** option (default 4) gives the desired number of block comparisons
per call to **daligner**(1). Some must contain *dal*-1 comparisons, and the
first *dal*-2 block comparisons even less, but the **HPCdaligner** "planner"
does the best it can to give an average load of dal block comparisons per
command. The **-deg** option (default 25) gives the maximum number of files
that will be merged in a single **LAmerge**(1) command. The planner makes the
most even k-ary tree of merges, where the number of levels is ceil(log_deg N).

If the integers *first* and *last* are missing, then the script produced
is for every block in the database. If *first* is present, then
**HPCdaligner** produces an incremental script that compares blocks *first*
through *last* (*last* = *first* if not present) against each other and all
previous blocks 1 through *first*-1, and then incrementally updates the .las
files for blocks 1 through *first*-1, and creates the .las files for
blocks *first* through *last*.

Each UNIX command line output by the **HPCdaligner** can be a batch job
(we use the && operator to combine several commands into one line to make this
so). Dependencies between jobs can be maintained simply by first running all
the **daligner**(1) jobs, then all the initial sort jobs, and then all the
jobs in each phase of the external merge sort. Each of these phases is
separated by an informative comment line for your scripting convenience.

# SEE ALSO

**daligner**(1)
**LAsort**(1)
**LAmerge**(1)
**LAshow**(1)
**LAcat**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCmapper**(1)
