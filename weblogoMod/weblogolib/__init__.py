#!/usr/bin/env python

# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

# Replicates README.txt

"""
WebLogo (http://code.google.com/p/weblogo/) is a tool for creating sequence 
logos from biological sequence alignments.  It can be run on the command line,
as a standalone webserver, as a CGI webapp, or as a python library.

The main WebLogo webserver is located at http://weblogo.threeplusone.com

Please consult the manual for installation instructions and more information:
(Also located in the weblogolib/htdocs subdirectory.)

    http://weblogo.threeplusone.com/manual.html

For help on the command line interface run
    ./weblogo --help

To build a simple logo run
    ./weblogo  < cap.fa > logo0.eps
    
To run as a standalone webserver at localhost:8080 
    ./weblogo --serve

To create a logo in python code:
    >>> from weblogolib import *
    >>> fin = open('cap.fa')
    >>> seqs = read_seq_data(fin) 
    >>> data = LogoData.from_seqs(seqs)
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options)
    >>> fout = open('cap.eps', 'w') 
    >>> eps_formatter( data, format, fout)


-- Distribution and Modification --
This package is distributed under the new BSD Open Source License. 
Please see the LICENSE.txt file for details on copyright and licensing.
The WebLogo source code can be downloaded from 
http://code.google.com/p/weblogo/

WebLogo requires Python 2.5, 2.6 or 2.7, and the python
array package 'numpy' (http://www.scipy.org/Download)

Generating logos in PDF or bitmap graphics formats require that the ghostscript
program 'gs' be installed. Scalable Vector Graphics (SVG) format also requires 
the program 'pdf2svg'.

"""

import sys
import copy
import os
from datetime import datetime
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from math import sqrt

from  weblogoMod.corebio.data import rna_letters, dna_letters, amino_acid_letters


from string import Template
from subprocess import *
from weblogoMod.corebio.utils import resource_string, resource_filename

from math import log, sqrt, exp

# Avoid 'from numpy import *' since numpy has lots of names defined
from numpy import array, asarray, float64, ones, zeros, int32,all,any, shape
import numpy as na


from color import *
from colorscheme import *
from weblogoMod.corebio.seq import Alphabet, Seq, SeqList
from weblogoMod.corebio import seq_io
from weblogoMod.corebio.utils import isfloat, find_command, ArgumentError, stdrepr
from weblogoMod.corebio.moremath import *
from weblogoMod.corebio.data import amino_acid_composition
from weblogoMod.corebio.seq import unambiguous_rna_alphabet, unambiguous_dna_alphabet, unambiguous_protein_alphabet

import weblogoMod.corebio

from logomath import Dirichlet


# ------ META DATA ------

__all__ = [ 'LogoOptions', 
            'description', 
            '__version__', 
            'LogoFormat',
            'LogoData',
            'GhostscriptAPI',
            'std_color_schemes',
            'default_color_schemes',
            'classic',
            'std_units',
            'std_sizes',
            'std_alphabets',
            'std_percentCG',
            'pdf_formatter',
            'jpeg_formatter',
            'png_formatter',
            'png_print_formatter',
            'txt_formatter',
            'eps_formatter',
            'formatters',
            'default_formatter',
            'base_distribution',
            'equiprobable_distribution',
            'read_seq_data',
            'color',
            'colorscheme',
            'logomath',
            ]

description  = "Create sequence logos from biological sequence alignments." 

__version__ = weblogoMod.corebio.__version__

# These keywords are substituted by subversion.
# The date and revision will only tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
release_date ="$Date: 2012-07-02 19:28:12 -0700 (Mon, 02 Jul 2012) $".split()[1]
release_build = "$Revision: 145 $".split()[1]
release_description = "WebLogo %s (%s)" % (__version__,  release_date)



def cgi(htdocs_directory) :
    import weblogolib._cgi
    weblogolib._cgi.main(htdocs_directory)
        
class GhostscriptAPI(object) :
    """Interface to the command line program Ghostscript ('gs')"""
    
    formats = ('png', 'pdf', 'jpeg')
    
    def __init__(self, path=None) :
        try:
            command = find_command('gs', path=path)
        except EnvironmentError:
            try:
                command = find_command('gswin32c.exe', path=path)
            except EnvironmentError:
                raise EnvironmentError("Could not find Ghostscript on path."
                " There should be either a gs executable or a gswin32c.exe on your system's path")
        
        self.command = command        
    
    def version(self) :
        args = [self.command, '--version']
        try :
            p = Popen(args, stdout=PIPE)
            (out,err) = p.communicate() 
        except OSError :
            raise RuntimeError("Cannot communicate with ghostscript.")  
        return out.strip()
       
    def convert(self, format, fin, fout,  width,  height, resolution=300) :
        device_map = { 'png':'png16m',  'pdf':'pdfwrite', 'jpeg':'jpeg'}
       
        try :
            device = device_map[format]
        except KeyError:
            raise ValueError("Unsupported format.")
        
        args = [self.command, 
            "-sDEVICE=%s" % device, 
            "-dPDFSETTINGS=/screen", #Modification printer changed to screen
            #"-q",   # Quite: Do not dump messages to stdout.
            "-sstdout=%stderr", # Redirect messages and errors to stderr
            "-sOutputFile=-", # Stdout
            "-dUseCIEColor", #Modification
            "-dDEVICEWIDTHPOINTS=%s" % str(width),  
            "-dDEVICEHEIGHTPOINTS=%s" % str(height),  
            "-dSAFER",  # For added security
            "-dNOPAUSE",]
            
        if device != 'pdf' :
            args.append("-r%s" % str(resolution) ) 
            if resolution < 300 : # Antialias if resolution is Less than 300 DPI
                args.append("-dGraphicsAlphaBits=4")
                args.append("-dTextAlphaBits=4")
                args.append("-dAlignToPixels=0")
        
        args.append("-")  # Read from stdin. Must be last argument.
        
        error_msg = "Unrecoverable error : Ghostscript conversion failed " \
                    "(Invalid postscript?). %s" % " ".join(args) 

        source = fin.read()
        try :
            p = Popen(args, stdin=PIPE, stdout = PIPE, stderr= PIPE) 
            (out,err) = p.communicate(source) 
        except OSError :
            raise RuntimeError(error_msg)
    
        if p.returncode != 0 :
            print("COMMAND " + str(self.command))
            print("ARGS" +  str(args))
            error_msg += '\nReturn code: %i\n' % p.returncode 
            if err is not None : error_msg += err
            raise RuntimeError(error_msg)
        
        print >>fout, out
# end class Ghostscript


aa_composition = [ amino_acid_composition[_k] for _k in
                    unambiguous_protein_alphabet]



# ------  DATA ------

classic = ColorScheme([
    ColorGroup("G",  "orange" ),
    ColorGroup("TU", "red"),
    ColorGroup("C",  "blue"),
    ColorGroup("A",  "green")
    ] )
    
std_color_schemes = {"auto":            None,   # Depends on sequence type
                     "monochrome":      monochrome,
                     "base pairing":    base_pairing,
                     "classic":         classic,
                     "hydrophobicity" : hydrophobicity,  
                     "chemistry" :      chemistry,
                     "charge" :         charge,
                     }#

default_color_schemes = {
    unambiguous_protein_alphabet: hydrophobicity, 
    unambiguous_rna_alphabet: base_pairing,
    unambiguous_dna_alphabet: base_pairing 
}

std_units = {
    "bits"    : 1./log(2),
    "nats"    : 1.,
    "digits"  : 1./log(10),
    "kT"      : 1.,
    "kJ/mol"  : 8.314472 *298.15 /1000.,
    "kcal/mol": 1.987 *298.15  /1000.,
    "probability" : None,
}

# The base stack width is set equal to 9pt Courier. 
# (Courier has a width equal to 3/5 of the point size.)
# Check that can get 80 characters in journal page @small
# 40 characters in a journal column
std_sizes = {
    "small" : 5.4 , 
    "medium" : 5.4*2,
    "large"  : 5.4*3
}
            
std_alphabets = {
    'protein': unambiguous_protein_alphabet, 
    'rna': unambiguous_rna_alphabet,
    'dna': unambiguous_dna_alphabet}

std_percentCG = {
    'H. sapiens'    : 40.,
    'E. coli'       : 50.5,
    'S. cerevisiae' : 38.,
    'C. elegans'    : 36.,
    'D. melanogaster': 43.,
    'M. musculus'   :  42.,
    'T. thermophilus' : 69.4,
}

# Thermus thermophilus: Henne A, Bruggemann H, Raasch C, Wiezer A, Hartsch T,
# Liesegang H, Johann A, Lienard T, Gohl O, Martinez-Arias R, Jacobi C, 
# Starkuviene V, Schlenczeck S, Dencker S, Huber R, Klenk HP, Kramer W, 
# Merkl R, Gottschalk G, Fritz HJ: The genome sequence of the extreme 
# thermophile Thermus thermophilus.
# Nat Biotechnol 2004, 22:547-53
            

  
class LogoOptions(object) :
    """ A container for all logo formatting options. Not all of these
    are directly accessible through the CLI or web interfaces. 
    
    To display LogoOption defaults:
    >>> from weblogolib import *
    >>> LogoOptions()
    
    All physical lengths are measured in points. (72 points per inch, 28.3 points per cm)
      
    String attributes:
        o creator_text      -- Embedded as comment in figures.
        o logo_title             
        o logo_label
        o unit_name         -- See std_units for options. (Default 'bits') 
        o yaxis_label       -- Defaults to unit_name      
        o xaxis_label
        o fineprint         -- Defaults to WebLogo name and version
        
    Boolean attributes:
        o show_yaxis 
        o show_xaxis         
        o show_ends 
        o show_fineprint        
        o show_errorbars    -- Draw errorbars (default: False)
        o show_boxes        -- Draw boxes around stack characters (default: True)
        o debug             -- Draw extra graphics debugging information. 
        o rotate_numbers    -- Draw xaxis numbers with vertical orientation? 
        o scale_width       -- boolean, scale width of characters proportional to ungaps
        o pad_right         -- Make a single line logo the same width as multiline logos (default: False)                           
                                
    Other attributes:
        o stacks_per_line
        
        o yaxis_tic_interval 
        o yaxis_minor_tic_ratio 
        o yaxis_scale
        o xaxis_tic_interval 
        o number_interval

        o shrink_fraction       -- Proportional shrinkage of characters if show_boxes is true.

        o errorbar_fraction 
        o errorbar_width_fraction 
        o errorbar_gray 

        o resolution             -- Dots per inch (default: 96). Used for bitmapped output formats
        
        o default_color 
        o color_scheme 
        
        o stack_width           --
        o stack_aspect_ratio    -- Ratio of stack height to width (default: 5)

        o logo_margin           -- Default: 2 pts
        o stroke_width          -- Default: 0.5 pts
        o tic_length            -- Default: 5 pts
        o stack_margin          -- Default: 0.5 pts
        
        o small_fontsize        -- Small text font size in points
        o fontsize              -- Regular text font size in points
        o title_fontsize        -- Title text font size in points
        o number_fontsize       -- Font size for axis-numbers, in points.
        
        o text_font
        o logo_font
        o title_font
        
        o first_index
        o logo_start
        o logo_end

    """       

    def __init__(self, **kwargs) :
        """ Create a new LogoOptions instance.
        
        >>> L = LogoOptions(logo_title = "Some Title String")
        >>> L.show_yaxis = False
        >>> repr(L)
        """

        self.alphabet = None

        self.creator_text = release_description
        
        self.logo_title = ""
        self.logo_label = ""
        self.stacks_per_line = 40
        
        self.unit_name = "bits"
     
        self.show_yaxis = True
        # yaxis_lable default depends on other settings. See LogoFormat
        self.yaxis_label = None
        self.yaxis_tic_interval = 1.
        self.yaxis_minor_tic_ratio = 5
        self.yaxis_scale = None

        self.show_xaxis = True
        self.xaxis_label = ""
        self.xaxis_tic_interval =1
        self.rotate_numbers = False
        self.number_interval = 5            
        self.show_ends = False
        self.annotate = None
        
          
        self.show_fineprint = True
        self.fineprint =  "Based on WebLogo "+__version__ 
    
        self.show_boxes = False
        self.shrink_fraction = 0.5
  
        self.show_errorbars = True
        self.errorbar_fraction = 0.90
        self.errorbar_width_fraction = 0.25      
        self.errorbar_gray = 0.75
 
        self.resolution  = 96.     # Dots per inch
           
        self.default_color = Color.by_name("black")    
        self.color_scheme = None
        #self.show_color_key = False # NOT yet implemented
        
        self.debug = False
        
        self.logo_margin = 2
        self.stroke_width = 0.5
        self.tic_length = 5
        
        self.stack_width = std_sizes["large"]        
        self.stack_aspect_ratio = 5
        
        self.stack_margin = 0.5
        self.pad_right = False  

        self.small_fontsize = 6
        self.fontsize = 10
        self.title_fontsize = 12
        self.number_fontsize = 8

        self.text_font    = "ArialMT"
        self.logo_font    = "Arial-BoldMT"
        self.title_font   = "ArialMT"

        self.first_index = 1        
        self.logo_start = None      
        self.logo_end=None          
        self.scale_width = True

        self.reverse_stacks = True       # If true, draw stacks with largest letters on top.

        from weblogoMod.corebio.utils import update
        update(self, **kwargs)

    def __repr__(self) :
        from weblogoMod.corebio.util import stdrepr
        return stdrepr( self) 

    def __repr__(self) :
        attributes = vars(self).keys()
        attributes.sort()
        return stdrepr(self, attributes )
 
# End class LogoOptions


        

class LogoFormat(LogoOptions) :     
    """ Specifies the format of the logo. Requires LogoData and LogoOptions 
    objects.
    
    >>> data = LogoData.from_seqs(seqs )
    >>> options = LogoOptions()
    >>> options.title = "A Logo Title"
    >>> format = LogoFormat(data, options) 
    
    Raises an ArgumentError if arguments are invalid.
    """
    def __init__(self, data, options= None) :
        """ Create a new LogoFormat instance.
        
        """
        LogoOptions.__init__(self)
        
        if options is not None :
            self.__dict__.update(options.__dict__)
                 
        self.alphabet = data.alphabet
        self.seqlen = data.length
        
        
        # Derived parameters.
        self.show_title = False
        self.show_xaxis_label = False
        self.yaxis_minor_tic_interval = None
        self.lines_per_logo       = None
        self.char_width       = None        # Maximum character width. Stack width minus margins. 
        self.line_margin_left = None
        self.line_margin_right    = None
        self.line_margin_bottom = None
        self.line_margin_top = None
        self.title_height = None
        self.xaxis_label_height = None
        self.line_height = None
        self.line_width = None
        self.logo_height = None
        self.logo_width = None
        self.creation_date = None
        self.end_type = None

        self.stack_height = self.stack_width * self.stack_aspect_ratio
        
        # Attribute to test, test, error message
        arg_conditions = (
            ("stacks_per_line",     lambda x: x>0 ,     "Stacks per line must be positive."),
            ("stack_width",         lambda x: x>0.0,    "Stack width must be greater than zero."),
            ("stack_aspect_ratio" , lambda x: x>0,    "Stack aspect ratio must be greater than zero."),
            ("fontsize" ,           lambda x: x>0 ,     "Font sizes must be positive."),
            ("small_fontsize" ,     lambda x: x>0 ,     "Font sizes must be positive."),
            ("title_fontsize" ,     lambda x: x>0 ,     "Font sizes must be positive."),
            ("errorbar_fraction" ,  lambda x: x>=0.0 and x<=1.0, 
                "The visible fraction of the error bar must be between zero and one."),
            ("yaxis_tic_interval" , lambda x: x>=0.0 , "The yaxis tic interval cannot be negative."),                                    
            ("yaxis_minor_tic_interval" , lambda x: not (x and x<0.0) , "Distances cannot be negative."),
            ("xaxis_tic_interval" , lambda x: x>0.0 ,   "Tic interval must be greater than zero."),
            ("number_interval" ,    lambda x: x>0.0 ,      "Invalid interval between numbers."),
            ("shrink_fraction" ,    lambda x: x>=0.0 and x<=1.0 , "Invalid shrink fraction."),
            ("stack_margin" ,       lambda x: x>0.0 , "Invalid stack margin."),
            ("logo_margin" ,        lambda x: x>0.0 , "Invalid logo margin."),
            ("stroke_width",       lambda x: x>0.0 , "Invalid stroke width."),
            ("tic_length" ,         lambda x: x>0.0 , "Invalid tic length."),
        )
        
        # Run arguments tests. The second, attribute argument to the ArgumentError is 
        # used by the UI to provide user feedback.
        # FIXME: More validation        
        for test in arg_conditions :
            if not test[1]( getattr(self,test[0]) ) : raise ArgumentError(test[2], test[0])
   
         
        # Inclusive upper and lower bounds
        # FIXME: Validate here. Move from eps_formatter        
        if self.logo_start is  None: self.logo_start = self.first_index
        
        if self.logo_end is  None : 
            self.logo_end = self.seqlen + self.first_index -1 
        
        self.total_stacks = self.logo_end - self.logo_start +1

        if self.logo_start - self.first_index <0 :
            raise ArgumentError(
                "Logo range extends before start of available sequence.",
                'logo_range')
        
        if self.logo_end - self.first_index  >= self.seqlen  : 
            raise ArgumentError(
                "Logo range extends beyond end of available sequence.",
                'logo_range')
    
        if self.logo_title      : self.show_title = True
        if not self.fineprint   : self.show_fineprint = False
        if self.xaxis_label     : self.show_xaxis_label = True

        if self.yaxis_label is None : 
            self.yaxis_label = self.unit_name
        
        if self.yaxis_label : 
            self.show_yaxis_label = True
        else :
            self.show_yaxis_label = False
            self.show_ends = False
              
        if not self.yaxis_scale :
            conversion_factor = std_units[self.unit_name]
            if conversion_factor :
                self.yaxis_scale=log(len(self.alphabet))*conversion_factor
            else :
                self.yaxis_scale=1.0    # probability units

        if self.yaxis_scale<=0.0 : 
            raise ArgumentError("Invalid yaxis scale", 'yaxis_scale',)

        if self.yaxis_tic_interval >= self.yaxis_scale:
            self.yaxis_tic_interval /= 2. 

        self.yaxis_minor_tic_interval \
            = float(self.yaxis_tic_interval)/self.yaxis_minor_tic_ratio
                      
        if self.color_scheme is None :
            if self.alphabet in default_color_schemes :
                self.color_scheme = default_color_schemes[self.alphabet]
            else :
                self.color_scheme = monochrome

        self.lines_per_logo = 1+ ( (self.total_stacks-1) / self.stacks_per_line)
    
        if self.lines_per_logo==1 and not self.pad_right:
            self.stacks_per_line = min(self.stacks_per_line, self.total_stacks)

        self.char_width = self.stack_width - 2* self.stack_margin
    
    
        if self.show_yaxis :
            self.line_margin_left = self.fontsize * 3.0
        else :
            self.line_margin_left = 0

        if self.show_ends :
            self.line_margin_right = self.fontsize *1.5 
        else :
            self.line_margin_right = self.fontsize             

        if self.show_xaxis :
            if self.rotate_numbers :
                self.line_margin_bottom = self.number_fontsize *2.5
            else:
                self.line_margin_bottom = self.number_fontsize *1.5
        else :
            self.line_margin_bottom = 4

        self.line_margin_top = 4
    
        if self.show_title :
            self.title_height = self.title_fontsize 
        else :
            self.title_height = 0

        self.xaxis_label_height =0.
        if self.show_xaxis_label :
            self.xaxis_label_height += self.fontsize
        if self.show_fineprint :
            self.xaxis_label_height += self.small_fontsize

        self.line_height = (self.stack_height + self.line_margin_top +    
                            self.line_margin_bottom )
        self.line_width  = (self.stack_width*self.stacks_per_line + 
                            self.line_margin_left + self.line_margin_right )

        self.logo_height = int(2*self.logo_margin + self.title_height \
            + self.xaxis_label_height + self.line_height*self.lines_per_logo) 
        self.logo_width = int(2*self.logo_margin + self.line_width )


        self.creation_date = datetime.now().isoformat(' ') 

        end_type = '-'
        end_types = {
            unambiguous_protein_alphabet: 'p', 
            unambiguous_rna_alphabet: '-',
            unambiguous_dna_alphabet: 'd' 
        }
        if self.show_ends and self.alphabet in end_types:
            end_type = end_types[self.alphabet]
        self.end_type = end_type

        
        if self.annotate is None :
            self.annotate = []
            for i in range(self.seqlen):
                index = i + self.first_index
                if index % self.number_interval == 0 : 
                    self.annotate.append( "%d"%index)
                else :
                    self.annotate.append("")
                    
        if len(self.annotate)!=self.seqlen  :
            raise ArgumentError(
                  "Annotations must be same length as sequences.",
                  'annotate')
            

    # End __init__
# End class LogoFormat



# ------ Logo Formaters ------
# Each formatter is a function f(LogoData, LogoFormat, output file).
# that draws a representation of the logo into the given file.
# The main graphical formatter is eps_formatter. A mapping 'formatters'
# containing all available formatters is located after the formatter
# definitions. 

def pdf_formatter(data, format, fout) :
    """ Generate a logo in PDF format."""
    
    feps = StringIO()
    eps_formatter(data, format, feps)
    feps.seek(0)
    
    gs = GhostscriptAPI()    
    gs.convert('pdf', feps, fout, format.logo_width, format.logo_height)


def _bitmap_formatter(data, format, fout, device) :
    feps = StringIO()
    eps_formatter(data, format, feps)
    feps.seek(0)
    
    gs = GhostscriptAPI()    
    gs.convert(device, feps, fout, 
        format.logo_width, format.logo_height, format.resolution)


def jpeg_formatter(data, format, fout) : 
    """ Generate a logo in JPEG format."""
    _bitmap_formatter(data, format, fout, device="jpeg")

def svg_formatter(data, format, fout) : 
    """ Generate a logo in Scalable Vector Graphics (SVG) format.
    Requires the program 'pdf2svg' be installed.
    """

    fpdf = StringIO()
    pdf_formatter(data, format, fpdf)
    fpdf.seek(0)
    
    try:
   	    command = find_command('pdf2svg') 
    except EnvironmentError:
    	raise EnvironmentError("Scalable Vector Graphics (SVG) format requires the program 'pdf2svg'. "
    	                        "Cannot find 'pdf2svg' on search path.")

    import tempfile, os
    fpdfi, fname_pdf = tempfile.mkstemp(suffix=".pdf")
    fsvgi, fname_svg = tempfile.mkstemp(suffix=".svg")
    try:
        
        fpdf2 = open(fname_pdf, 'w')
        fpdf2.write(fpdf.getvalue() )
        fpdf2.seek(0)
  
        args = [command, fname_pdf, fname_svg]
        p = Popen(args)
        (out,err) = p.communicate() 

        fsvg = open(fname_svg)
        fout.write(fsvg.read())
    finally:
        os.remove(fname_svg)
        os.remove(fname_pdf)


def png_formatter(data, format, fout) : 
    """ Generate a logo in PNG format."""
    _bitmap_formatter(data, format, fout, device="png")


def png_print_formatter(data, format, fout) : 
    """ Generate a logo in PNG format with print quality (600 DPI) resolution."""
    format.resolution = 600
    _bitmap_formatter(data, format, fout, device="png")


def txt_formatter( logodata, format, fout) :
    """ Create a text representation of the logo data. 
    """
    print >>fout, str(logodata)

   

    
def eps_formatter( logodata, format, fout) :
    """ Generate a logo in Encapsulated Postscript (EPS)"""
    
    substitutions = {}
    from_format =[
        "creation_date",    "logo_width",           "logo_height",      
        "lines_per_logo",   "line_width",           "line_height",
        "line_margin_right","line_margin_left",     "line_margin_bottom",
        "line_margin_top",  "title_height",         "xaxis_label_height",
        "creator_text",     "logo_title",           "logo_margin",
        "stroke_width",     "tic_length",           
        "stacks_per_line",  "stack_margin",
        "yaxis_label",      "yaxis_tic_interval",   "yaxis_minor_tic_interval",
        "xaxis_label",      "xaxis_tic_interval",   "number_interval",
        "fineprint",        "shrink_fraction",      "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",    "small_fontsize",       "fontsize",
        "title_fontsize",   "number_fontsize",      "text_font",
        "logo_font",        "title_font",          
        "logo_label",       "yaxis_scale",          "end_type",
        "debug",            "show_title",           "show_xaxis",
        "show_xaxis_label", "show_yaxis",           "show_yaxis_label",
        "show_boxes",       "show_errorbars",       "show_fineprint",
        "rotate_numbers",   "show_ends",            "stack_height",
        "stack_width"
        ]
   
    for s in from_format :
        substitutions[s] = getattr(format,s)

    substitutions["shrink"] = str(format.show_boxes).lower()


    # --------- COLORS --------------
    def format_color(color):
        return  " ".join( ("[",str(color.red) , str(color.green), 
            str(color.blue), "]"))  

    substitutions["default_color"] = format_color(format.default_color)

    colors = []  
    for group in format.color_scheme.groups :
        cf = format_color(group.color)
        for s in group.symbols :
            colors.append( "  ("+s+") " + cf )
    substitutions["color_dict"] = "\n".join(colors)
        
    data = []
    
    # Unit conversion. 'None' for probability units
    conv_factor = std_units[format.unit_name]
    
    data.append("StartLine")

    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index 
        stack_index = seq_index - seq_from
        
        if stack_index!=0 and (stack_index % format.stacks_per_line) ==0 :
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")
        data.append("0 0 0 setrgbcolor\n(%s) StartStack" % format.annotate[seq_index] )
        # if format.annotate[seq_index][-1] == "*":
        #     data.append("0 0 1 setrgbcolor\n(%s) StartStack" % format.annotate[seq_index] )
        # else:
        #     data.append("0 0 0 setrgbcolor\n(%s) StartStack" % format.annotate[seq_index] )

        if conv_factor: 
            stack_height = logodata.entropy[seq_index] * std_units[format.unit_name]
        else :
            stack_height = 1.0 # Probability

        s = zip(logodata.counts[seq_index], logodata.alphabet)
        def mycmp( c1, c2 ) :
            # Sort by frequency. If equal frequency then reverse alphabetic
            if c1[0] == c2[0] : return cmp(c2[1], c1[1])
            return cmp(c1[0], c2[0])
        
        s.sort(mycmp)
        if not format.reverse_stacks: s.reverse()

        C = float(sum(logodata.counts[seq_index])) 
        if C > 0.0 :
            fraction_width = 1.0
            if format.scale_width :
                fraction_width = logodata.weight[seq_index] 
            # print >>sys.stderr, fraction_width
            for c in s:
                data.append(" %f %f (%s) ShowSymbol" % (fraction_width, c[0]*stack_height/C, c[1]) )

        # Draw error bar on top of logo. Replaced by DrawErrorbarFirst above.
        if logodata.entropy_interval is not None and conv_factor and C>0.0:

            low, high = logodata.entropy_interval[seq_index]
            center = logodata.entropy[seq_index]
            low *= conv_factor
            high *= conv_factor
            center *=conv_factor
            if high> format.yaxis_scale : high = format.yaxis_scale 

            down = (center - low) 
            up   = (high - center) 
            data.append(" %f %f DrawErrorbar" % (down, up) )
            
        data.append("EndStack")
        data.append("")
               
    data.append("EndLine")
    substitutions["logo_data"] = "\n".join(data)  


    # Create and output logo
    template = resource_string( __name__, 'template.eps', __file__)
    logo = Template(template).substitute(substitutions)
    print >>fout, logo
 

# map between output format names and logo  
formatters = {
    'eps': eps_formatter, 
    'pdf': pdf_formatter,
    'png': png_formatter,
    'png_print' : png_print_formatter,
    'jpeg'  : jpeg_formatter,
    'svg'   : svg_formatter,
    'logodata' : txt_formatter,          
    }     
    
default_formatter = eps_formatter





def parse_prior(composition, alphabet, weight=None) :
    """ Parse a description of the expected monomer distribution of a sequence.
    
    Valid compositions:
    
    - None or 'none' :      No composition sepecified 
    - 'auto' or 'automatic': Use the typical average distribution
                            for proteins and an equiprobable distribution for
                            everything else.    
    - 'equiprobable' :      All monomers have the same probability.
    - a percentage, e.g. '45%' or a fraction '0.45':
                            The fraction of CG bases for nucleotide alphabets
    - a species name, e.g. 'E. coli', 'H. sapiens' :
                            Use the average CG percentage for the specie's      
                            genome.
    - An explicit distribution,  e.g. {'A':10, 'C':40, 'G':40, 'T':10}
    """
    if composition is None: return None
    comp = composition.strip()
    
    if comp.lower() == 'none': return None
    
    
    if weight is None and alphabet is not None: 
        weight = sqrt(float(len(alphabet)))

    if weight<0 : raise ValueError("Weight cannot be negative.")
    
    
    if comp.lower() == 'equiprobable' :
        prior = weight * equiprobable_distribution(len(alphabet)) 
    elif comp.lower() == 'auto' or comp.lower() == 'automatic':
        if alphabet == unambiguous_protein_alphabet :
            prior =  weight * asarray(aa_composition, float64)
        else :
            prior = weight * equiprobable_distribution(len(alphabet)) 
    
    elif comp in std_percentCG :
        prior = weight * base_distribution(std_percentCG[comp])

    elif comp[-1] == '%' :
        prior = weight * base_distribution( float(comp[:-1]))

    elif isfloat(comp) :
        prior = weight * base_distribution( float(comp)*100. )

    elif composition[0] == '{' and composition[-1] == '}' : 
        explicit = composition[1: -1]
        explicit = explicit.replace(',',' ').replace("'", ' ').replace('"',' ').replace(':', ' ').split()
        
        if len(explicit) != len(alphabet)*2 :
            #print explicit
            raise ValueError("Explicit prior does not match length of alphabet")
        prior = - ones(len(alphabet), float64) 
        try :
            for r in range(len(explicit)/2) :
                letter = explicit[r*2]
                index = alphabet.ord(letter)
                value = float(explicit[r*2 +1])
                prior[index] = value
        except ValueError :
            raise ValueError("Cannot parse explicit composition")
    
        if any(prior==-1.) :
            raise ValueError("Explicit prior does not match alphabet") 
        prior/= sum(prior)
        prior *= weight
        
        
    else : 
        raise ValueError("Unknown or malformed composition: %s"%composition)
    
    if len(prior) != len(alphabet) :
        raise ValueError(
            "The sequence alphabet and composition are incompatible.")
    return prior

    
def base_distribution(percentCG) :
    A = (1. - (percentCG/100.))/2.
    C = (percentCG/100.)/2.
    G = (percentCG/100.)/2.
    T = (1. - (percentCG/100))/2.
    return asarray((A,C,G,T), float64)   

def equiprobable_distribution( length) :
    return ones( (length), float64) /length   
  

def read_seq_data(lines, 
                input_parser=seq_io.read, 
                alphabet=None, 
                ignore_lower_case=False, 
                max_file_size=0):
    """ Read sequence data from the input stream and return a seqs object. 
    
    The environment variable WEBLOGO_MAX_FILE_SIZE overides the max_file_size argument.
    Used to limit the load on the WebLogo webserver.
    """

    seqs = input_parser(lines)

    if seqs is None or len(seqs) ==0 :
        raise ValueError("Please provide a multiple sequence alignment")
    
    if ignore_lower_case :
        # Case is significant. Do not count lower case letters.
        for i,s in enumerate(seqs) :
            seqs[i] = s.mask()

    # Add alphabet to seqs.
    if alphabet :
        seqs.alphabet = alphabet 
    else :
        seqs.alphabet = Alphabet.which(seqs)
    return seqs


    
class LogoData(object) :
    """The data needed to generate a sequence logo.
       
    - alphabet 
    - length
    - counts  -- An array of character counts
    - entropy -- The relative entropy of each column
    - entropy_interval -- entropy confidence interval
     """
     
    def __init__(self, length=None, alphabet = None, counts =None, 
            entropy =None, entropy_interval = None, weight=None) :
        """Creates a new LogoData object"""
        self.length = length
        self.alphabet = alphabet
        self.counts = counts
        self.entropy = entropy
        self.entropy_interval = entropy_interval
        self.weight = weight
        

    @classmethod    
    def from_counts(cls, alphabet, counts, prior= None):
        """Build a LogoData object from counts."""
        # Counts is a Motif object?
        #counts = counts.array
        
        seq_length, A = counts.shape
        
        if prior is not None: prior = array(prior, float64)
        
        if prior is None or sum(prior)==0.0:
            R = log(A)
            ent = zeros(  seq_length, float64)
            entropy_interval = None    
            for i in range (0, seq_length) :
                C = sum(counts[i]) 
                #FIXME: fixup corebio.moremath.entropy()?
                if C == 0 :
                    ent[i] = 0.0
                else :
                    ent[i] = R - entropy(counts[i])
        else :
            ent = zeros(  seq_length, float64)
            entropy_interval = zeros( (seq_length,2) , float64)
        
            R = log(A)
            
            for i in range (0, seq_length) :
                alpha = array(counts[i] , float64)
                alpha += prior
                
                posterior = Dirichlet(alpha)
                ent[i] = posterior.mean_relative_entropy(prior/sum(prior)) 
                entropy_interval[i][0], entropy_interval[i][1] = \
                    posterior.interval_relative_entropy(prior/sum(prior), 0.95) 
 
        weight = array( na.sum(counts,axis=1) , float) 
        weight /= max(weight)
 
        return cls(seq_length, alphabet, counts, ent, entropy_interval, weight)


    @classmethod    
    def from_seqs(cls, seqs, prior= None):
        """Build a LogoData object from a SeqList, a list of sequences."""
        # --- VALIDATE DATA ---
        # check that at least one sequence of length at least 1 long
        if len(seqs)==0 or len(seqs[0]) ==0:
            raise ValueError("No sequence data found.")   
    
        # Check sequence lengths    
        seq_length = len(seqs[0])
        for i,s in enumerate(seqs) :
            #print i,s, len(s)
            #TODO: Redundant? Should be checked in SeqList?
            if seq_length != len(s) :
                raise ArgumentError(
    "Sequence number %d differs in length from the previous sequences" % (i+1) ,
                        'sequences')

        # FIXME: Check seqs.alphabet?
        
        counts = seqs.profile()
        return cls.from_counts(seqs.alphabet, counts, prior)


    def __str__(self) :
        out = StringIO()
        print >>out, '## LogoData'
        print >>out, '# First column is position number, counting from zero'
        print >>out, '# Subsequent columns are raw symbol counts'
        print >>out, '# Entropy is mean entropy measured in nats.' 
        print >>out, '# Low and High are the 95% confidence limits.'
        print >>out, '# Weight is the fraction of non-gap symbols in the column.'
        print >>out, '#\t'
        print >>out, '#\t',
        for a in self.alphabet :
            print >>out, a, '\t',
        print >>out, 'Entropy\tLow\tHigh\tWeight'
        
        for i in range(self.length) :
            print >>out, i+1, '\t',
            for c in self.counts[i] : print >>out, c, '\t',
            print >>out, "%6.4f" % self.entropy[i], '\t',
            if self.entropy_interval is not None:
                print >>out, "%6.4f" % self.entropy_interval[i][0], '\t', 
                print >>out, "%6.4f" % self.entropy_interval[i][1], '\t',
            else :
                print >>out, '\t','\t',
            if self.weight is not None :
                print >>out, "%6.4f" % self.weight[i],
            print >>out, ''
        print >>out, '# End LogoData'
        
        return out.getvalue()





