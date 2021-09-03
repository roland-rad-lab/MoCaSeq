"""
***********************************************
- PROGRAM: gjainPyLib.py
***********************************************
"""
__version_info__ = ('1', '0', '1')
__version__      = '.'.join(__version_info__)

import string
import glob
import random
import matplotlib
import itertools
import types
import argparse
import os.path
import sys
import textwrap
import re
import tempfile
import os
import atexit
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess
import pprint as pp
import scipy.stats as stats
import seaborn as sns
import pip
from sinfo import sinfo

#######################################################
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
#######################################################

def import_or_install(package):
    '''
     Import a package, where package is of type str,
     and if it is unable to, calls pip and attempt
     to install it from there
     source: https://stackoverflow.com/questions/4527554/check-if-module-exists-if-not-install-it
    '''
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', '--user', package])

def get_file_info(file_name_with_path):
    ''' Get path, basename, ext, path+basename and basename+ext of a file '''

    # get the path of the file
    input_file_path = os.path.dirname(os.path.abspath(file_name_with_path))

    # get the filename without path but with extension
    input_file_name = os.path.basename(file_name_with_path)

    # get basename of the file (no path or extension)
    input_file_basename = os.path.splitext(input_file_name)[0]

    # get extension of the file
    input_file_extension = os.path.splitext(input_file_name)[1]
    # remove "."
    input_file_extension = re.sub(r"^\.","",input_file_extension)

    # get the path+filename
    path_filename = input_file_path + "/" + input_file_basename

    return (input_file_path, input_file_basename, input_file_extension, path_filename, input_file_name)

def create_dir(dirname):
    '''Creates a directory if not exists. Returns relevant status messages on error'''
    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise

def get_fileslist_in_dir(dirname, *exts):
    '''
    	1) Get the list of files inside a directory
    	2) If extension is passed then return files for that extension
    '''

    # Get all files in the directory
    filesList = [ f for f in os.listdir(dirname) if os.path.isfile(os.path.join(dirname, f))]

    # If extension is passed then parse files with the provided extension (.bed, .txt)
    specificFiles = list()
    for ext in exts:
        files_ext = [i for i in filesList if i.endswith(".{0}".format(ext))]
        specificFiles.extend(files_ext)

    if exts:
        return specificFiles
    else:
        return filesList

def print_initial_arguments(parser):
    ''' Prints all the initial arguments entered '''

    print("\n---------- Script ----------")
    print(sys.argv[0])

    print("\n---------- Session info ----------")
    print(sinfo())
    print("*"*50)

    print("\n------Input Arguments ------")
    opts = vars(parser.parse_args())
    maxl = max(len(text) for text in opts.keys())
    for k,v in opts.items():
        print("%-*s = %s" % (maxl, k, v))
    print("-"*29, "\n")

def print_help():
    ''' Print system help '''
    print >> sys.stderr, "\n ----------------- HELP ------------------\n", parser.print_help(), "\n"

def nanmean(data, **args):
    '''Get the means of a vector with nans'''
    return np.ma.filled(np.ma.masked_array(data,np.isnan(data)).mean(**args), fill_value=np.nan)

def useful_lines(fh):
    '''Skips all the comments (right now #) from the file and returns the generator object of lines without comments.
    It remembers the state where it left. So it will after the last executed line '''
    for line in (l.strip() for l in fh if not l.startswith('#')):
        yield line

def _get_unique_list(seq, idfun=None):
    ''' Originally proposed by Andrew Dalke '''
    seen = set()
    if idfun is None:
        for x in seq:
            if x not in seen:
                seen.add(x)
                yield x
    else:
        for x in seq:
            x = idfun(x)
            if x not in seen:
                seen.add(x)
                yield x

def get_unique_list(seq, idfun=None):
    '''
     - Get unique element in a list with the order of elements preserved
     - Source: https://stackoverflow.com/questions/4459703/how-to-make-lists-contain-only-distinct-element-in-python
    '''
    # Order preserving
    return list(_get_unique_list(seq, idfun))

def get_line_counts(filename):
    '''
     - Get the number of lines in a file efficiently
     - Source: https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
     - rawincounts()
    '''
    from itertools import (takewhile,repeat)
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )

def _get_line_counts_makegen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)
def get_line_counts_pygen(filename):
    '''
     - Get the number of lines in a file efficiently
     - Source: https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
     - rawgencount()
     - For huge files like 100+ GB in size Using a separate generator function, this runs a smidge faster:
    '''
    f = open(filename, 'rb')
    f_gen = _get_line_counts_makegen(f.raw.read)
    return sum( buf.count(b'\n') for buf in f_gen )

def get_line_counts_mapcount(filename):
    '''
     - Get the number of lines in a file efficiently
     - Source: https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
     - mapcount()
    '''
    import time
    import mmap
    import random
    from collections import defaultdict
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

def get_temp_file(suffix='', prefix='tmp', dir=None):
    '''
     - Create a temporary file that will be removed at exit
     - Returns a path to the file
    '''
    _f, path = tempfile.mkstemp(suffix, prefix, dir)
    os.close(_f)
    remove_at_exit(path)
    return path

# Register a file to be removed at exit
def remove_at_exit(path):
    atexit.register(os.remove, path)

# class for Logging
class Log(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        pass
