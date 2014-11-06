/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Local alignment module.  Routines for finding local alignments given a seed position,
 *    representing such an l.a. with its interval and a set of pass-thru points, so that
 *    a detailed alignment can be efficiently computed on demand.
 *
 *  All routines work on a numeric representation of DNA sequences, i.e. 0 for A, 1 for C,
 *    2 for G, and 3 for T.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#include "DB.h"

#define TRACE_XOVR 125   //  If the trace spacing is not more than this value, then can
                         //    and do compress traces pts to 8-bit unsigned ints

#ifndef _A_MODULE

#define _A_MODULE

/*** PATH ABSTRACTION:

     Coordinates are *between* characters where 0 is the tick just before the first char,
     1 is the tick between the first and second character, and so on.  Our data structure
     is called a Path refering to its conceptualization in an edit graph.

     A local alignment is specified by the point '(abpos,bbpos)' at which its path in
     the underlying edit graph starts, and the point '(aepos,bepos)' at which it ends.
     In otherwords A[abpos+1..aepos] is aligned to B[bbpos+1..bepos] (assuming X[1] is
     the *first* character of X).

     There are 'diffs' differences in an optimal local alignment between the beginning and
     end points of the alignment (if computed by Compute_Trace), or nearly so (if computed
     by Local_Alignment).  

     Optionally, a Path can have additional information about the exact nature of the
     aligned substrings if the field 'trace' is not NULL.  Trace points to either an
     array of integers (if computed by a Compute_Trace routine), or an array of unsigned
     short integers (if computed by Local_Alignment).

     If computed by Local_Alignment 'trace' points at a list of 'tlen' (= n+1) short
     values:

            b_0, b_1, ... b_n-1, b_n

     to be interpreted as follows.  The alignment from (abpos,bbpos) to (aepos,bepos)
     passes through the n trace points for i in [1,n]:

            (a_i,b_i) where a_i = floor(abpos/TS)*TS + i*TS
                        and b_i = bbpos + (b_0 + b_1 + b_i-1)

     where also let a_0,b_0 = abpos,bbpos and a_(n+1),b_(n+1) = aepos,bepos.  That is, the
     interior (i.e. i != 0 and i != n+1) trace points pass through every TS'th position of
     the aread where TS is the "trace spacing" employed when finding the alignment (see
     New_Align_Spec).  Typically TS is 100.  These trace points allow the Compute_Trace
     routines to efficiently compute the exact alignment between the two reads by efficiently
     computing exact alignments between consecutive pairs of trace points.

     If computed by a Compute_Trace routine, 'trace' points at a list of 'tlen' integers
     < i1, i2, ... in > that encodes an exact alignment as follows.  A negative number j
     indicates that a dash should be placed before A[-j] and a positive number k indicates
     that a dash should be placed before B[k], where A and B are the two sequences of the
     overlap.  The indels occur in the trace in the order in which they occur along the
     alignment.  For a good example of how to "decode" a trace into an alignment, see the
     code for the routine Print_Alignment.

***/

typedef struct
  { void     *trace;
    READIDX   tlen;
    READIDX   diffs;
    READIDX   abpos, bbpos;
    READIDX   aepos, bepos;
  } Path;


/*** ALIGNMENT ABSTRACTION:

     An alignment is modeled by an Alignment record, which in addition to a *pointer* to a
     'path', gives pointers to the A and B sequences, their lengths, and indicates whether
     the B-sequence needs to be complemented ('comp' non-zero if so).  The 'trace' pointer
     of the 'path' subrecord can be either NULL, a list of pass-through points, or an exact
     trace depending on what routines have been called on the record.

     One can (1) compute a trace, with Compute_Trace, either from scratch if 'path.trace' = NULL,
     or using the sequence of pass-through points in trace, (2) print an ASCII representation
     of an alignment, or (3) reverse the roles of A and B, and (4) complement a sequence
     (which is a reversible process).

     If the alignment record shows the B sequence as complemented, *** THEN IT IS THE
     RESPONSIBILITY OF THE CALLER *** to make sure that bseq points at a complement of
     the sequence before calling Compute_Trace or Print_Alignment.  Complement_Seq complements
     the sequence a.  The operation does the complementation/reversal in place.  Calling it a
     second time on a given fragment restores it to its original state.
***/

#define COMP(x)  ((x) & 0x1)

#define COMP_FLAG 0x1

typedef struct
  { Path *path;
    char   *aseq;         /* Pointer to A sequence                              */
    char   *bseq;         /* Pointer to B sequence                              */
    READIDX alen;         /* Length of A sequence                               */
    READIDX blen;         /* Length of B sequence                               */
    int     flags;        /* Pipeline status and complementation flags          */
  } Alignment;

void Complement_Seq(char *a);

  /* Many routines like Local_Alignment, Compute_Trace, and Print_Alignment need working
     storage that is more efficiently reused with each call, rather than being allocated anew
     with each call.  Each *thread* can create a Work_Data object with New_Work_Data and this
     object holds and retains the working storage for routines of this module between calls
     to the routines.  Free_Work_Data frees a Work_Data object and all working storage
     held by it.
  */

  typedef void Work_Data;

  Work_Data *New_Work_Data();

  void       Free_Work_Data(Work_Data *work);

  /* Local_Alignment seeks local alignments of a quality determined by a number of parameters.
     These are coded in an Align_Spec object that can be created with New_Align_Spec and
     freed with Free_Align_Spec when no longer needed.  There are 4 essential parameters:

     ave_corr:    the average correlation (1 - 2*error_rate) for the sought alignments.  For Pacbio
                    data we set this to .70 assuming an average of 15% error in each read.
     trace_space: the spacing interval for keeping trace points and segment differences (see
                    description of 'trace' for Paths above)
     freq[4]:     a 4-element vector where afreq[0] = frequency of A, f(A), freq[1] = f(C),
                    freq[2] = f(G), and freq[3] = f(T).  This vector is part of the header
                    of every HITS database (see db.h).

     If an alignment cannot reach the boundary of the d.p. matrix with this condition (i.e.
     overlap), then the last/first 30 columns of the alignment are guaranteed to be
     suffix/prefix positive at correlation ave_corr * g(freq) where g is an empirically
     measured function that increases from 1 as the entropy of freq decreases.

     You can get back the original parameters used to create an Align_Spec with the simple
     utility functions below.
  */

  typedef void Align_Spec;

  Align_Spec *New_Align_Spec(double ave_corr, int trace_space, float *freq);

  void        Free_Align_Spec(Align_Spec *spec);

  int    Trace_Spacing      (Align_Spec *spec);
  double Average_Correlation(Align_Spec *spec);
  float *Base_Frequencies   (Align_Spec *spec);

  /* Local_Alignment finds the longest significant local alignment between the sequences in
     'align' subject to the alignment criterion given by the Align_Spec 'spec' that passes
     through the point '(x,y)' within the underlying dynamic programming matrix.  The path
     record of 'align' has its 'trace' filled from the point of view of an overlap between
     the aread and the bread.  In addition a Path record from the point of view of the bread
     versus the aread is returned by the function, with this Path's 'trace' filled in
     appropriately.  The space for the returned path and the two 'trace's are in the
     working storage supplied by the Work_Data packet and this space is reused with each call,
     so if one wants to retain the bread-path and the two trace point sequences, then they
     must be copied to user-allocated storage before calling the routine again.
  */

  Path *Local_Alignment(Alignment *align, Align_Spec *spec, Work_Data *work, int x, int y);

  /* Given a legitimate Alignment object, Compute_Trace_X computes an exact trace for the alignment.
     If 'path.trace' is non-NULL, then it is assumed to be a sequence of pass-through points
     and diff levels computed by Local_Alignment.  In either case 'path.trace' is set
     to point at an integer array within the storage of the Work_Data packet encoding an
     exact optimal trace from the start to end points.  If the trace is needed beyond the
     next call to a routine that sets it, then it should be copied to an array allocated
     and managed by the caller.

     Compute_Trace_ALL does not require a sequence of pass-through points, as it computes the
     best alignment between (path->abpos,path->bbpos) and (path->aepos,path->bepos) in the
     edit graph between the sequences.  Compute_Trace_PTS computes a trace by computing the
     trace between successive pass through points.  It is much, much faster than Compute_Trace_ALL
     but at the tradeoff of not necessarily being optimal as pass-through points are not all
     perfect.  Compute_Trace_MID computes a trace by computing the trace between the mid-points
     of alignments between two adjacent pairs of pass through points.  It is generally twice as
     slow as Compute_Trace_PTS, but it produces nearer optimal alignments.
  */

  void Compute_Trace_ALL(Alignment *align, Work_Data *work);
  void Compute_Trace_PTS(Alignment *align, Work_Data *work, int trace_spacing);
  void Compute_Trace_MID(Alignment *align, Work_Data *work, int trace_spacing);

  /* Print_Acartoon prints an ASCII representation of the overlap relationhip between the
     two reads of 'align' to the given 'file' indented by 'indent' space.

     If the alignment trace is an exact trace, then one can ask Print_Alignment to print an
     ASCII representation of the alignment 'align' to the file 'file'.  Indent the display
     by "indent" spaces and put "width" columns per line in the display.  Show "border"
     characters of sequence on each side of the aligned region.  If upper is non-zero then
     display bases in upper case.  If coord is greater than 0, then the positions of the
     first character in A and B in the given rwo is displayed with a field width given by
     coord's value.

     Print_Reference is like Print_Alignment but rather than printing exaclty "width" columns
     per segment, it prints "block" characters of the A sequence in each segment.  This results
     in segments of different lengths, but is convenient when looking at two alignments involving
     A as segments are guaranteed to cover the same interval of A in a segment.
  */

  void Print_ACartoon(FILE *file, Alignment *align, int indent);

  void Print_Alignment(FILE *file, Alignment *align, Work_Data *work,
                       int indent, int width, int border, int upper, int coord);

  void Print_Reference(FILE *file, Alignment *align, Work_Data *work,
                       int indent, int block, int border, int upper, int coord);


/*** OVERLAP ABSTRACTION:

     Externally, between modules an Alignment is modeled by an "Overlap" record, which
     replaces the pointers to the two sequences with their ID's in the HITS data bases,
     and contains its path as a subrecord rather than as a pointer (indeed, typically the
     corresponding Alignment record points at the Overlap's path sub-record).  One can read
     and write binary records of an "Overlap" and produce an ASCI print-outs of them.
***/

typedef struct {
  Path    path;         /* Path: begin- and end-point of alignment + diffs    */
  int     aread;        /* Id # of A sequence                                 */
  int     bread;        /* Id # of B sequence                                 */
  READIDX alen;         /* Length of A sequence                               */
  READIDX blen;         /* Length of B sequence                               */
  int     flags;        /* Pipeline status and complementation flags          */
} Overlap;

  /* Read_Overlap reads the next Overlap record from stream 'input', not including the trace
     (if any), and without modifying 'ovl's trace pointer.  Read_Trace reads the ensuing trace
     into the memory pointed at by the trace field of 'ovl'.  It is assumed to be big enough to
     accommodate the trace where each value take 'tbytes' bytes.

     Write_Overlap write 'ovl' to stream 'output' followed by its trace vector (if any) that
     occupies 'tbytes' bytes per value.  

     Print_Overlap prints an ASCII version of the contents of 'ovl' to stream 'output'
     indented from the left margin by 'indent' spaces.

     Compress_TraceTo8 converts a trace fo 16-bit values to 8-bit values in place, and
     Decompress_TraceTo16 does the reverse conversion.

     Check_Trace_Points checks that the number of trace points is correct and that the sum
     of the b-read displacements equals the b-read alignment interval, assuming the trace
     spacing is 'tspace'.  It reports an error message if there is a problem and 'verbose'
     is non-zero.  The 'ovl' came from the file names 'fname'.

     Print_OCartoon is the same as Print_ACartoon except it takes an Overlap pointer.
  */

  int Read_Overlap(FILE *input, Overlap *ovl);
  int Read_Trace(FILE *innput, Overlap *ovl, int tbytes);

  void Write_Overlap(FILE *output, Overlap *ovl, int tbytes);
  void Print_Overlap(FILE *output, Overlap *ovl, int indent);

  void Compress_TraceTo8(Overlap *ovl);
  void Decompress_TraceTo16(Overlap *ovl);

  int  Check_Trace_Points(Overlap *ovl, int tspace, int verbose, char *fname);

  void Print_OCartoon(FILE *file, Overlap *ovl, int indent);

#endif // _A_MODULE
