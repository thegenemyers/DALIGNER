% DALIGNER(1) 1.0
%
% August 2015

# NAME

daligner - long read aligner

# SYNOPSIS

**daligner**
[**-vbAI**]
[**-k***int(14)*] [**-w***int(6)*] [**-h***int(35)*]
[**-t***int*] [**-M***int*]
[**-e***double(.70)*] [**-l***int(1000)*] [**-s***int(100)*] [**-H***int*]
[**-m***track*]+ *subject:db|dam* *target:db|dam* ...

# DESCRIPTION
Compare sequences in the trimmed *subject* block against those in the list of
*target* blocks searching for local alignments involving at least **-l** base
pairs (default 1000) or more, that have an average correlation rate of **-e**
(default 70%). The local alignments found will be output in a sparse encoding
where a trace point on the alignment is recorded every **-s** base pairs of
the a-read (default 100bp). Reads are compared in both orientations and local
alignments meeting the criteria are output to one of several created files
described below. The **-v** option turns on a verbose reporting mode that gives
statistics on each major step of the computation.

The options **-k**, **-h**, and **-w** control the initial filtration search
for possible matches between reads. Specifically, our search code looks for a
pair of diagonal bands of width 2^w (default 2^6 = 64) that contain a
collection of exact matching k-mers (default 14) between the two reads,
such that the total number of bases covered by the k-mer hits is h
(default 35). k cannot be larger than 32 in the current implementation.
If the **-b** option is set, then the **daligner** assumes the data has a
strong compositional bias (e.g. >65% AT rich), and at the cost of a bit more
time, dynamically adjusts k-mer sizes depending on compositional bias, so that
the mers used have an effective specificity of 4^k.

If there are one or more interval tracks specified with the **-m** option, then
the reads of the DB or DB's to which the mask applies are soft masked with
the union of the intervals of all the interval tracks that apply, that is any
k-mers that contain any bases in any of the masked intervals are ignored for
the purposes of seeding a match. An interval track is a track, such as the
"dust" track created by DBdust, that encodes a set of intervals over either
the untrimmed or trimmed DB.

Invariably, some k-mers are significantly over-represented (e.g. homopolymer
runs). These k-mers create an excessive number of matching k-mer pairs and
left unaddressed would cause daligner to overflow the available physical
memory.  One way to deal with this is to explicitly set the **-t** parameter
which suppresses the use of any k-mer that occurs more than *t* times in either
the subject or target block.  However, a better way to handle the situation is
to let the program automatically select a value of *t* that meets a given
memory usage limit specified (in Gb) by the **-M** parameter. By default
**daligner** will use the amount of physical memory as the choice for **-M**.
If you want to use less, say only 8Gb on a 24Gb HPC cluster node because you
want to run 3 **daligner** jobs on the node, then specify **-M***8*.
Specifying **-M***0* basically indicates that you do not want **daligner** to
self adjust k-mer suppression to fit within a given amount of memory.  

For each subject, target pair of blocks, say X and Y, the program reports
alignments where the a-read is in X and the b-read is in Y, and vice versa.
However, if the **-A** option is set ("A" for "asymmetric") then just overlaps
where the a-read is in X and the b-read is in Y are reported, and if X = Y,
then it further reports only those overlaps where the a-read index is less than the b-read index.  In either case, if the **-I** option is set
("I" for "identity") then when X = Y, overlaps between different
portions of the same read will also be found and reported.  

Each found alignment is recorded as -- a[ab,ae] x bo[bb,be] -- where a and b
are the indices (in the trimmed DB) of the reads that overlap, o indicates
whether the b-read is from the same or opposite strand, and [ab,ae] and
[bb,be] are the intervals of a and bo, respectively, that align. The program
places these alignment records in files whose name is of the form
X.Y.[C|N]#.las where C indicates that the b-reads are complemented and N
indicates they are not (both comparisons are performed) and # is
the thread that detected and wrote out the collection of alignments contained
in the file. That is the file X.Y.O#.las contains the alignments produced by
thread # for which the a-read is from X and the b-read is from Y and in
orientation O. The command
**daligner -A** *X* *Y* produces 2\*NTHREAD thread files X.Y.?.las and
**daligner** *X* *Y* produces 4\*NTHREAD files X.Y.?.las and Y.X.?.las
(unless *X*=*Y* in which case only NTHREAD files, X.X.?.las, are produced).

By default, **daligner** compares all overlaps between reads in the database
that are greater than the minimum cutoff set when the DB or DBs were split,
typically 1 or 2 Kbp. However, the HGAP assembly pipeline only wants to
correct large reads, say 8Kbp or over, and so needs only the overlaps where
the a-read is one of the large reads. By setting the **-H** parameter to say N,
one alters **daligner** so that it only reports overlaps where the a-read is
over N base-pairs long.

While the default parameter settings are good for raw Pacbio data, **daligner**
can be used for efficiently finding alignments in corrected reads or other
less noisy reads. For example, for mapping applications against .dams, we run

**daligner** **-k**20 **-h**60 **-e**.85

and on corrected reads, we typically run

**daligner** **-k**25 **-w**5 **-h**60 **-e**.95 **-s**500

and at these settings it is very fast.

# SEE ALSO

**LAsort**(1)
**LAmerge**(1)
**LAshow**(1)
**LAcat**(1)
**LAsplit**(1)
**LAcheck**(1)
**HPCdaligner**(1)
**HPCmapper**(1)
