from ctypes import *
from DAPI import *

_READIDX = c_uint16
_TRACE_XOVR = 125

class Path(Structure):
    _fields_ = [("trace", POINTER(c_uint16)),
                ("tlen", _READIDX),
                ("diffs", _READIDX),
                ("abpos", _READIDX),
                ("bbpos", _READIDX),
                ("aepos", _READIDX),
                ("bepos", _READIDX)]


class Alignment(Structure):
    _fields_ = [("path", POINTER(Path)),
                ("aseq", POINTER(c_char)),
                ("bseq", POINTER(c_char)),
                ("alen", c_int),
                ("blen", c_int),
                ("flag", c_int)]

class Overlap(Structure):
    _fields_ = [("path", Path),
                ("aread", c_int),
                ("bread", c_int),
                ("alen", _READIDX),
                ("blen", _READIDX),
                ("flags", c_int)]

libc = CDLL("libc.so.6")

fopen = libc.fopen
fclose = libc.fclose
fread = libc.fread
fread.argtypes = [c_void_p, c_size_t, c_size_t, c_void_p]

"""void *realloc(void *ptr, size_t size);"""
realloc = libc.realloc
realloc.argtypes = [c_void_p, c_size_t]
realloc.restype = c_void_p

"""void *malloc(size_t size);"""
malloc = libc.malloc
malloc.argtypes = [c_size_t]
malloc.restype = c_void_p

ptr_size = sizeof(c_void_p)
ovl_IO_size = sizeof(Overlap) - ptr_size

def _read_overlap(in_f, ovl):
    p = ovl.path.trace
    fread( cast( addressof(ovl) + ptr_size, c_void_p ), ovl_IO_size, 1, in_f )
    ovl.path.trace = p

def _read_trace(in_f, ovl, tbytes):
    fread( cast(ovl.path.trace, c_void_p), tbytes, ovl.path.tlen, in_f )

def get_ovl_data(fn):

    in_f = fopen(fn, "r")

    novl = c_int64()
    tspace = c_int()
    fread(addressof(novl) , sizeof(c_int64), 1, in_f)
    fread(addressof(tspace) , sizeof(c_int), 1, in_f)

    if tspace.value < _TRACE_XOVR:
        small  = 1
        tbytes = sizeof(c_uint8)
    else:
        small  = 0
        tbytes = sizeof(c_uint16)

    tmax = 1000
    trace = cast( malloc( sizeof(c_uint16) * tmax ), POINTER(c_uint16) )

    ovl = Overlap()
    ovl_data = {}

    for j in xrange(novl.value):
        _read_overlap(in_f, ovl)

        if ovl.path.tlen > tmax:
            tmax = 1.2*ovl.path.tlen + 100
            trace = cast( realloc( trace, sizeof(c_uint16) * tmax ),  POINTER(c_uint16) )

        ovl.path.trace = trace
        _read_trace(in_f, ovl, tbytes)

        if ovl.alen < 8000:
            continue
        if ovl.path.abpos > 50 and ovl.path.bbpos > 50:
            continue
        if ovl.alen - ovl.path.aepos > 50 and ovl.blen - ovl.path.bepos > 50:
            continue
        comp = ovl.flags & 0x1
        bbpos, bepos  = ovl.path.bbpos, ovl.path.bepos
        if comp == 1:
            bbpos, bepos = ovl.blen - bepos, ovl.blen - bbpos
        acc = 100 - (200.0 * ovl.path.diffs / ( ovl.path.aepos - ovl.path.abpos + ovl.path.aepos - ovl.path.abpos ))

        ovl_data.setdefault(ovl.aread,[])

        ovl_data[ovl.aread].append( (ovl.aread, ovl.bread, acc, ovl.path.abpos,  ovl.path.aepos, ovl.alen, comp, bbpos, bepos, ovl.blen) )
    fclose(in_f)

    return ovl_data
