from ctypes import *

_READIDX = c_uint16
_TRACE_XOVR = 125

class HITS_READ(Structure):
    _fields_ = [ ("origin", c_int),
                 ("begin", _READIDX),
                 ("end", _READIDX),
                 ("boff", c_int64),
                 ("coff", c_int64),
                 ("flags", c_int)]


class HITS_TRACK(Structure):
    pass

HITS_TRACK._fields_ = [ ("_track", POINTER(HITS_TRACK)),
                        ("name", c_char_p),
                        ("size", c_int),
                        ("anno", c_void_p),
                        ("data", c_void_p)]


class HITS_DB(Structure):
    _fields_ = [ ( "oreads", c_int ),
                 ( "breads", c_int ),
                 ( "cutoff", c_int ),
                 ( "all", c_int),
                 ( "freq", c_float * 4),
                 ( "maxlen", c_int),
                 ( "totlen", c_int64),
                 ( "nreads", c_int),
                 ( "trimmed", c_int),
                 ( "part", c_int),
                 ( "ofirst", c_int),
                 ( "bfirst", c_int),
                 ( "path", c_char_p),
                 ( "loaded", c_int),
                 ( "bases", c_void_p),
                 ( "reads", POINTER(HITS_READ)),
                 ( "tracks", POINTER(HITS_TRACK)) ]


DB = CDLL("./DB.so")

libc = CDLL("libc.so.6")

fopen = libc.fopen
fclose = libc.fclose
fread = libc.fread
fread.argtypes = [c_void_p, c_size_t, c_size_t, c_void_p]

open_DB = DB.Open_DB
open_DB.argtypes = [c_char_p, POINTER(HITS_DB)]
open_DB.restype = c_int

load_read  = DB.Load_Read
load_read.argtypes = [POINTER(HITS_DB), c_int, c_char_p, c_int]
load_read.restype = c_int

close_DB = DB.Close_DB
close_DB.argtypes = [POINTER(HITS_DB)]
close_DB.restype = c_int

trim_DB = DB.Trim_DB
trim_DB.argtypes = [POINTER(HITS_DB)]
trim_DB.restype = c_int

new_read_buffer = DB.New_Read_Buffer
new_read_buffer.argtypes = [ POINTER(HITS_DB) ]
new_read_buffer.restype = POINTER(c_char)
