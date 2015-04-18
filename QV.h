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
 *  Compressor/decompressor for .quiv files: customized Huffman codes for each stream based on
 *    the histogram of values occuring in a given file.  The two low complexity streams
 *    (deletionQV and substitutionQV) use a Huffman coding of the run length of the prevelant
 *    character.
 *
 *  Author:   Gene Myers
 *  Date:     Jan 18, 2014
 *  Modified: July 25, 2014
 *
 ********************************************************************************************/

#ifndef _QV_COMPRESSOR

#define _QV_COMPRESSOR

  //  A PacBio compression scheme

typedef struct
  { void    *delScheme;   //  Huffman scheme for deletion QVs
    void    *insScheme;   //  Huffman scheme for insertion QVs
    void    *mrgScheme;   //  Huffman scheme for merge QVs
    void    *subScheme;   //  Huffman scheme for substitution QVs
    void    *dRunScheme;  //  Huffman scheme for deletion run lengths (if delChar > 0)
    void    *sRunScheme;  //  Huffman scheme for substitution run lengths (if subChar > 0)
    int      delChar;     //  If > 0, run-encoded deletion value
    int      subChar;     //  If > 0, run-encoded substitution value
    int      flip;        //  Need to flip multi-byte integers
    char    *prefix;      //  Header line prefix
  } QVcoding;

  // Read the next nlines of input, and QVentry returns a pointer to the first line if needed.

int       Read_Lines(FILE *input, int nlines);
char     *QVentry();

  // Read the .quiva file on input and record frequency statistics.

void     QVcoding_Scan(FILE *input);

  // Given QVcoding_Scan has been called at least once, create an encoding scheme based on
  //   the accumulated statistics and return a pointer to it.  The returned encoding object
  //   is *statically* allocated within the routine.  If lossy is set then use a lossy scaling
  //   for the insertion and merge streams.

QVcoding *Create_QVcoding(int lossy);

  //  Read/write a coding scheme to input/output.  The encoding object returned by the reader
  //    is *statically* allocated within the routine.

QVcoding *Read_QVcoding(FILE *input);
void      Write_QVcoding(FILE *output, QVcoding *coding);

  //  Free all the auxiliary storage associated with coding (but not the object itself!)

void      Free_QVcoding(QVcoding *coding);

  //  Assuming the file pointer is positioned just beyond an entry header line, read the
  //    next set of 5 QV lines, compress them according to 'coding', and output.  If lossy
  //    is set then the scheme is a lossy one.

void      Compress_Next_QVentry(FILE *input, FILE *output, QVcoding *coding, int lossy);

  //  Assuming the input is positioned just beyond the compressed encoding of an entry header,
  //    read the set of compressed encodings for the ensuing 5 QV vectors, decompress them,
  //    and place their decompressed values into entry which is a 5 element array of character
  //    pointers.  The parameter rlen computed from the preceeding header line, critically
  //    provides the length of each of the 5 vectors.

void      Uncompress_Next_QVentry(FILE *input, char **entry, QVcoding *coding, int rlen);

#endif // _QV_COMPRESSOR
