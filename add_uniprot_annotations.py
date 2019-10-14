"""add_protein_annotations.py
Adds richer annotations from UniProt Swiss-Prot DAT files (flat text format)
to proteins specified as a list of accession strings read in from the clipboard.
Annatations are added as table columns and written back to the clipboard in the
original order. This program is designed to function as a companion to Excel so
that annotations can be added along with any other processing needed to prepare
PAW results files for clients.

MIT License

Copyright (c) 2019 Phillip Wilmarth, OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

written by Kyra Patton, OHSU, summer 2016.
additional changes by Phil Wilmarth, PSR Core, OHSU, 2016.

Version 3-4:
    Debugged keywords processing, added some GUI borders.
    -PW 8/4/2016
Version 5:
    DAT file parsing results pickled.
    Started some keyword frequency analysis (need to save results)
    Added exclusions for less informative keywords
    -PW 8/5/2016
Version 6:
    Supports Python 2 or 3
    Added support for GO terms
    -KP 8/10/2016

    changed accessions from clipboard parsing
    -PW 12/10/2018

Renamed to "add_uniprot_annotations.py" for distribution on Github -PW 20191008
    Changed some behavior to facilitate annotating more lists in one session
    Works with 3-species DAT files (see "keywlist_download.py")
    
To-Do:    
    KW lines can wrap and then ECO codes don't get removed
       need to parse out all terms and put the lines together (10/6 PW)

Completed:
    DONE - Make KW and GO reports
    DONE - Implement CC -!- PATHWAY information and correlate with REACTOME info
    DONE - Add support for BLAST ortholog mapping
"""

"""Issue 6/22/2018:
Excel Mac 2016 seems to have a bunch of unicode junk at end of clipboard contents
The last read accession gets messed up. A blank cell at end of selection does not fix.
Have to add a physical dummy accession. I was not getting good formatted data back on the
clipboard for export.

ISSUE SEEMS TO BE RESOLVED FROM MICROSOFT IN NEWER EXCEL VERSIONS.
"""
import os
import sys
import gzip
import re
import time
try:
    import cPickle as pickle
except ImportError:
    import pickle

from tkinter import *
from tkinter import filedialog

import numpy as np
import pandas as pd


# module-wide function definitions
def get_file(default_location, ext_list=[('All files', '*.*')], title_string="Select a file"):
    """Dialog box to browse to a folder.  Returns full file path."""
    # set up GUI elements
    root = Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass

    # check location
    if not os.path.exists(default_location):
        default_location = os.getcwd()
        
    # create dialog box for file selection and return selection
    root.update()
    selection = filedialog.askopenfilename(parent=root, initialdir=default_location,
                                           filetypes=ext_list, title=title_string)
    return selection
    
def get_folder(default_location, title_string=""):
    """Dialog box to browse to a folder.  Returns full folder path."""    
    # set up GUI elements
    root = Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string if not passed
    if title_string == "":   
        title_string = 'Select parent folder containing desired folders'
    
    # create dialog box for folder selection
    root.update()   # helps make sure dialog box goes away after selection
    full_folder_name = filedialog.askdirectory(parent=root, initialdir=default_location,
                                               title=title_string, mustexist=True)    
    # return full folder name
    return full_folder_name   
    
    
# class definitions:            
class OneKeyWord:
    """Data container for one UniProt keyword definition."""
    
    def __init__(self):
        """Attributes and method name dictionary."""
        self.identifier = None  # keyword name
        self.accession = None   # keyword accession
        self.definition = None  # keyword description
        self.synonyms = []      # keyword synonyms
        self.GO_mapping = {}    # cross-reference to GO term
        self.hierarchy = None   # heirarchy
        self.web = None         # web link
        self.category = None    # one of ten keywordcategories

        # method switchyard
        self.methods = {'ID   ': self.parse_ID,
                        'IC   ': self.parse_ID,
                        'AC   ': self.parse_AC,
                        'DE   ': self.parse_DE,
                        'SY   ': self.parse_SY,
                        'GO   ': self.parse_GO,
                        'CA   ': self.parse_CA}
        return
                        
    def parse_ID(self, content):
        """Parses ID (or IC) lines."""
        self.identifier = content.rstrip('.')
        return
        
    def parse_AC(self, content):
        """Parses AC lines."""
        self.accession = content
        return
        
    def parse_DE(self, content):
        """Parses DE lines."""
        content = content.replace(os.linesep, ' ')
        self.definition = content
        return
        
    def parse_SY(self, content):
        """Parses SY lines."""
        content = content.replace(os.linesep, ' ')
        self.synonyms = [x.strip().rstrip('.') for x in content.split(';') if x.strip()]
        return
        
    def parse_GO(self, content):
        """Parses GO lines."""
        for line in content.split(os.linesep):
            parts = [x for x in line.split('GO:') if x]
            for part in parts:
                self.GO_mapping['GO:' + part.split(';')[0]] = part.split(';')[1].strip()
        return
        
    def parse_CA(self, content):
        """Parses AC lines."""
        self.category = content.rstrip('.')
        return

class KeyWords:
    """Container for parsed keywords from UniProt keyword documentation files."""
        
    def __init__(self):
        """Basic constructor."""
        self.keywords = {}      # store OneKeyWord objects in dictionary by keyword
        self.categories = []    # alpabetical list of keyword categories 
        return

    def parse_file(self, keyword_file):
        """Parses a keyword list file and creates OneKeyWord objects."""
        blocks = {}
        keyword = OneKeyWord() # object for first record
        for line in (open(keyword_file, 'r')):
            if line.startswith('ID   '): # reset blocks once we are into records
                blocks = {}
            if line.startswith("//"):   # record delimiter occurs after records
                self._parse_blocks(keyword, blocks)
                self.keywords[keyword.identifier] = keyword
                keyword = OneKeyWord() # reset object for next record
                continue # this way we don't get "//" lines in blocks
            key, value = line[:5], line[5:].rstrip() # separate code from content
            if key in blocks:
                blocks[key] += os.linesep + value # add platform line separator for multiple lines
            else:
                blocks[key] = value
        self.set_categories()
        return
                
    def _parse_blocks(self, keyword, blocks):
        """Calls parsing methods for each block in blocks dictionary."""
        for key, block in blocks.items():
            try:
                keyword.methods[key](block) 
            except KeyError:
                pass
        return                   

    def set_categories(self):
        """Gets the list of categories from the OneKeyWord objects."""
        self.categories = sorted(set([self.keywords[k].category for k in self.keywords]))
        
    def put_keywords_in_categories(self, keyword_list):
        """Given a list of keywords, returns a list of keywords grouped into 10 categories."""
        by_category_list = [[] for _ in self.categories]
        for keyword in keyword_list:
            index = self.categories.index(self.keywords[keyword].category)
            by_category_list[index].append(keyword)
        return ['; '.join(sorted(x)) for x in by_category_list]
        
    def summary_stats(self):
        """Computes frequencies of category terms."""
        import collections
        cat_freq = collections.Counter([self.keywords[k].category for k in self.keywords])
        print('\nCategory Frequencies')
        for cat in sorted(cat_freq):
            print('...category: %s, frequency: %i' % (cat, cat_freq[cat]))
        print()
        return 
        
class Annotations:
    """Object containing all of the different annotations for a protein record."""
    def __init__(self):
        # define the attributes
        self.identifier = None      # UniProt identifier string
        self.db = None              # "sp" or "tr"
        self.accession = None       # primary accession
        self.other_accessions = []  # other accessions
        self.fasta_accession = None # compund accession like in the FASTA files
        self.name = None            # primary protein name
        self.other_names = []       # alternative protein names
        self.flags = []             # flags from DE lines
        self.gene = None              # UniProt gene name
        self.other_genes = []       # other gene synonyms
        self.os = None              # species name
        self.ox = None              # taxonomy number
        self.mgi_acc = None         # mouse gene index cross-reference
        self.mgi_gene = None        # MGI gene name (may differ from UniProt)
        self.keywords = []          # keyword list
        self.go = GOTerms()         # GO term object
#        self.cc = None
        self.pathway = PathWays()   # container for CC PATHWAY and Reactome annotations

        # compile re pattern for evidence codes
        self.eco = re.compile(r' {?ECO:(.)*[},]') # ' {ECO:' followed by zero or more characters then '}'
        
        # non-informative keywords to exclude:
        self.excluded_keywords = ['Reference proteome', 'Complete proteome',
                                  'Direct protein sequencing']

        # methods switchyard dictionary
        self.parse_method = {'ID': self.get_identifier,
                             'AC': self.get_accessions,
                             'DE': self.get_names,
                             'GN': self.get_gene_name,
                             'OS': self.get_os,
                             'OX': self.get_ox,
                             'CC': self.get_cc,
                             'DR': self.get_databases,
                             'KW': self.get_keywords}
        return

    def parse_record(self, prot_rec):
        """Parses annotation fields from protein records.
        prot_rec: a list of strings, protein record"""        
        # index the record to speed parsing
        idx = self.make_index(prot_rec)

        # skip to the right part of the record and call its method
        for key in idx:
            try:
                self.parse_method[key](prot_rec[idx[key]:])
            except KeyError:
                pass    # skip annotation lines we are not parsing        
        self.fasta_accession = '|'.join([self.db, self.accession, self.identifier])        
        return        

    def make_index(self, prot_rec):
        """Indexes the major sections of the protein record.'
        prot_rec: a list of strings, protein record"""
        index = {}
        for i, line in enumerate(prot_rec):
            try:
                if line[0].isalpha() and line[1].isalpha() and line[2:5] == '   ':
                    key = line[0:2]
                    if key not in index:
                        index[key] = i
            except IndexError:
                pass
        return index
        
    def get_identifier(self, prot_rec):
        """Gets identifier and DB from ID line.
        prot_rec: a list of strings, protein record"""
        for line in prot_rec:
            if line.startswith('ID   '):
                self.identifier = line[5:].split()[0]
                if line[5:].split()[1].rstrip(';') == 'Reviewed':
                    self.db = 'sp'
                else:
                    self.db = 'tr'
                return
    
    def get_accessions(self, prot_rec):
        """Gets primary accession and list of secondary accessions.
        prot_rec: a list of strings, protein record"""
        other_accessions = []
        for line in prot_rec:
            if line.startswith('AC   '):
                other_accessions += [acc.strip() for acc in line[5:].split(';') if acc.strip()]
            else:
                break
        self.accession = other_accessions.pop(0)
        self.other_accessions = other_accessions
        return
        
    def get_names(self, prot_rec):
        """Gets protein names (primary and alternatives, full and short).
        prot_rec: a list of strings, protein record"""
        short_names = []
        for line in prot_rec:
            if not line.startswith('DE'):
                break
            line = self.eco.sub('', line[5:])     # lines start with two uppercase letters and 3 spaces
            if 'RecName:' in line:
                self.name = line.split('Full=')[1].rstrip(';')
            elif 'AltName:' in line:
                if short_names and not self.other_names:
                    self.name += ' (' + '; '.join(short_names) + ')'
                    short_names = []
                elif short_names:
                    self.other_names[-1] += ' (' + '; '.join(short_names) + ')'
                    short_names = []
                try:
                    self.other_names.append(line.split('Full=')[1].rstrip(';'))
                except IndexError:
                    self.other_names.append(line.split('AltName:')[1].strip().rstrip(';'))
            elif 'Flags:' in line:
                self.flags = [x.strip() for x in line.split(':')[1].split(';') if x.strip()]
            elif 'Short=' in line:
                short_names.append(line.split('Short=')[1].rstrip(';'))
    
        # may have some short names still to process
        if short_names and not self.other_names:
            self.name += ' (' + '; '.join(short_names) + ')'
            return
        elif short_names:
            self.other_names[-1] += ' (' + '; '.join(short_names) + ')'
            return
    
    def get_gene_name(self, prot_rec):
        """Return gene name from DAT file.
        prot_rec: a list of strings, protein record"""
        for line in prot_rec:
            if line.startswith('GN   '):
                line = self.eco.sub('', line)   # remove evidence codes
                self.gene = line[5:].split(';')[0].replace('Name=', '')
                synonyms = [x for x in line[5:].split(';') if 'Synonyms' in x]
                if synonyms:
                    other_genes = synonyms[0].replace('Synonyms=', '')
                    self.other_genes = [x.strip() for x in other_genes.split(',')]
                return
    
    def get_os(self, prot_rec):
        """Return the organism species
        prot_rec: a list of strings, protein record"""
        for line in prot_rec:
            if line.startswith('OS   '):
                self.os = line[5:].rstrip('.')
                return
    
    def get_ox(self, prot_rec):
        """Return the Taxonomy ID
        prot_rec: a list of strings, protein record"""
        for line in prot_rec:
            if line.startswith('OX   '):
                self.ox = line.split('=')[1].split()[0].rstrip(';')
                return    

    def get_databases(self, prot_rec):
        """Calls any database cross reference block parsings."""
        self.pathway.parse_reactome(prot_rec)    
        self.get_mgi(prot_rec)
        self.get_GO(prot_rec)
        return
                
    def get_mgi(self, prot_rec):
        """Returns the MGI accession number and gene name from the DAT file.
        prot_rec: a list of strings, protein record"""
        for line in prot_rec:
            if line.startswith('DR   MGI;'):
                self.mgi_acc = line[5:].split(';')[1].strip()
                self.mgi_gene = line[5:].split(';')[2].strip().rstrip('.')
                return
        return "N/A"
    
    def get_GO(self, prot_rec):
        """Populates a GOTerm object.
        prot_rec: a list of strings, protein record"""
        self.go.parse_GO_terms(prot_rec)
        return

    def get_cc(self, prot_rec):
        """Parses CC PATHWAY info."""
        self.pathway.parse_cc_pathway(prot_rec)
       
    def get_keywords(self, prot_rec):
        """Return list of keywords from DAT file.
        prot_rec: a list of strings, protein record"""
        keyword_lst = []
        for line in prot_rec:
            if line.startswith('KW   '):
                line = self.eco.sub('', line)   # remove evidence codes
                [keyword_lst.append(k.strip()) for k in line[5:].split(';') 
                if (k.strip() and k.strip().replace('.', '') not in self.excluded_keywords)]
            else:
                break
        if keyword_lst:
            keyword_lst[-1] = keyword_lst[-1].rstrip('.')
        self.keywords = keyword_lst
        return
        
    def _snoop(self):
        """Prints itself for diagnostics."""
        print('identifier:', self.identifier)
        print('db:', self.db)
        print('main acc:', self.accession)
        print('other acc:', self.other_accessions)
        print('fasta acc:', self.fasta_accession)
        print('gene name:', self.gene_name)
        print('other genes:', self.other_genes)
        print('species:', self.os)
        print('taxonomy:', self.ox)
        print('protein name:', self.name)
        print('aternative names:', self.other_names)
        print('mgi acc:', self.mgi_acc)
        print('mgi gene:', self.mgi_gene)
        print('key words:', self.keywords)
        return        
        
class PathWays:
    """Container for Reactome and CC PATHWAY annotations."""
    
    def __init__(self):
        """Basoc constructor."""
        self.react_acc = []     # list of Reactome database keys
        self.react_desc = []    # list of Reactome description strings
        self.react_string = ''  # formatted reactome info
        self.cc_string = ''     # string to collect CC description lines
        
        # compile re pattern for evidence codes
        self.eco = re.compile(r' {?ECO:(.)*[},]') # ' {ECO:' followed by zero or more characters then '}'
        return
        
    def parse_reactome(self, prot_rec):
        """Parses Reatome DR lines into paired acc and desc lists."""
        for line in prot_rec:
            line = self.eco.sub('', line)
            if line.startswith('DR   Reactome;'):
                line_split = [x.strip() for x in line.split(';')]
                self.react_acc.append(line_split[1])
                if line_split[2].endswith('.'):
                    line_split[2] = line_split[2][:-1]
                self.react_desc.append(line_split[2])
        self.react_string = '; '.join(['%s {%s}' % (desc, acc) for (desc, acc) in zip(self.react_desc, self.react_acc)])
        return
        
    def parse_cc_pathway(self, prot_rec):
        """Parses CC PATHWAY lines from protein record."""
        in_pathway = False
        for line in prot_rec:
            line = self.eco.sub('', line)
            if '-!- PATHWAY:' in line:
                in_pathway = True
                self.cc_string = line.split('-!- PATHWAY: ')[1] + ' '
            elif in_pathway and ('-!-' in line or '----------' in line):
                break
            elif in_pathway and '-!-' not in line:
                if line == 'CC      .':
                    continue
                try:
                    self.cc_string += line.split('CC   ')[1].lstrip() + ' '
                except IndexError:
                    break
        if self.cc_string.endswith(' '):
            self.cc_string = self.cc_string[:-1]
        return
 
class GOTerms:
    """Object containing GO terms out of DAT file."""
    def __init__(self):
        """prot_rec: a list of strings, protein record."""
        self._go_num = []               # GO accession
        self._go_type = []              # GO category
        self._go_desc = []              # GO term
        self.molecular_function = ''    # collects MF GO terms
        self.cellular_component = ''    # collects CC GO terms
        self.biological_process = ''    # collects BP GO terms
        return
        
    def parse_GO_terms(self, prot_rec):
        """Retrieve the GO Terms from the DAT file
        prot_rec: a list of strings, protein record"""
        for line in prot_rec:
            if line.startswith('DR   GO;'):
                terms = line.split('; ')[1:]
                self._go_num.append(terms[0][3:])
                self._go_type.append(terms[1][0])
                self._go_desc.append(terms[1][2:])
        self._sort_GO_terms()
        return
                
    def _sort_GO_terms(self):
        """Sort the GO Terms based on their category:
        Biological process, Molecular function, or
        Cellular component."""
        molecular_function = []
        cellular_component = []
        biological_process = []
        for i, item in enumerate(self._go_type):
            if item == 'F':
                molecular_function.append('%s {GO:%s}' % (self._go_desc[i], self._go_num[i]))
            if item == 'C':
                cellular_component.append('%s {GO:%s}' % (self._go_desc[i], self._go_num[i]))
            if item == 'P':
                biological_process.append('%s {GO:%s}' % (self._go_desc[i], self._go_num[i]))
        
        # combine terms into single strings
        self.molecular_function = '; '.join(molecular_function)
        self.cellular_component = '; '.join(cellular_component)
        self.biological_process = '; '.join(biological_process)
        return        
                
class AnnotationPickle:
      """Container for parsed annotation dictionary of DAT file.
      Includes DAT file name, and DAT file creation date to test if pickle seems viable.
      """
      def __init__(self, dat_file, dat_date, annotate_dict):
          """Basic constructor."""
          self.dat_file = dat_file
          self.dat_date = dat_date
          self.annotate_dict = annotate_dict
          return
         
          
# GUI classes
class StatusBar(Frame):
    """window status bar class with set and clear methods.
    Code from "http://effbot.org/tkinterbook/text.htm"
    """
    def __init__(self, master):
        Frame.__init__(self, master)
        self.label = Label(self, bd=1, padx=5, pady=1, relief=SUNKEN, anchor=W)
        self.label.pack(fill=X)
        return
               
    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()
        return
            
    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()
        return     
    # end StatusBar class  

class ProteinAnnotator:
    """Object creating the main GUI window."""
    def __init__(self):
        self.root = Tk()
        self.root.title('Protein Annotator')
        self.myFrame = Frame(self.root)
        self.myFrame.pack()
        self.radio_var = IntVar()
        self.kw_var = IntVar()
        self.pw_var = IntVar()
        self.go_var = IntVar()
        self.sf_var = IntVar()
        
        # set default values
        self.radio_var.set(1)
        self.kw_var.set(1)
        self.pw_var.set(1)
        self.go_var.set(1)
        self.sf_var.set(0)
        
        # create a button toolbar
        self.toolbar = Frame(self.myFrame)

        self.b1 = self.make_toolbar_button('Get accessions', self.get_accessions)
        self.b2 = self.make_toolbar_button('Parse DAT file', self.parse_dat_file, width=13)
        self.b3 = self.make_toolbar_button('Blast mapping', self.blast_mapping)
        self.b4 = self.make_toolbar_button('Add annotations', self.add_annotations)
        self.b5 = self.make_toolbar_button('Reset', self.clear_data, width=8)
        self.b6 = self.make_toolbar_button('Help', self.print_help, width=8)
        self.b7 = self.make_toolbar_button('Quit', self.quit_me, width=8)

        self.toolbar.pack(side=TOP, fill=X)
        
        #create an option bar
        self.option_bar = Frame(self.myFrame)
        self.rb_frame = Frame(self.option_bar, bd=2, relief=SUNKEN)  # radiobuttons frame

        self.rb1 = self.make_radiobutton('human', 1)
        self.rb2 = self.make_radiobutton('mouse', 2)
        self.rb3 = self.make_radiobutton('arabidopsis', 3)
        
        self.rb_frame.pack(side=LEFT, fill=X, padx=5, pady=5)
                
        self.cb_frame = Frame(self.option_bar, bd=2, relief=SUNKEN)  # checkboxes frame

        self.options = IntVar()        
        self.cb4 = self.make_checkbutton('Summary Files', self.sf_var)
        self.cb3 = self.make_checkbutton('Pathways', self.pw_var)
        self.cb2 = self.make_checkbutton('GO Terms', self.go_var)
        self.cb1 = self.make_checkbutton('Keywords', self.kw_var)



        
        self.cb_frame.pack(side=RIGHT, fill=X)
        self.option_bar.pack(side=TOP, fill=X, padx=5, pady=5)
        
        # add a multi-line text box with scrollbars
        self.textFrame = Frame(self.myFrame, bd=2, relief=SUNKEN)
        self.textFrame.grid_rowconfigure(0, weight=1)
        self.textFrame.grid_columnconfigure(0, weight=1)
        
        self.text = Text(self.textFrame, height=40, width=132,
                         wrap=NONE, padx=5, pady=5)
        self.xscroll = Scrollbar(self.textFrame, orient=HORIZONTAL,
                                 command=self.text.xview)
        self.xscroll.grid(row=1, column=0, sticky=W+E)
        self.yscroll = Scrollbar(self.textFrame, command=self.text.yview)
        self.yscroll.grid(row=0, column=1, sticky=N+S)
        
        self.text.configure(xscrollcommand=self.xscroll.set)
        self.text.configure(yscrollcommand=self.yscroll.set)
        self.text.grid(row=0, column=0, sticky=N+S+E+W)
        
        self.textFrame.pack()

        # add a status line
        self.status = StatusBar(self.myFrame)
        self.status.pack(side=BOTTOM, fill=X)
        self.status.set("%s", "Status line")
        self.print_help()

        # define the actual structures for the data
        self.accessions = []        # holds list of accessions
        self.acc_read = False       # flag for if accessions are loaded
        self.dat_file = None        # DAT file path and name
        self.dat_read = False       # flag for if DAT file parsed
        self.default = os.getcwd()  # can set a default location here
        self._annotate_dict = {}    # maps all possible accessions to annotations
        self.annotations = []       # list of matching annotations 
        self.blast_map = {}         # optional BLAST ortholog mapping
        self.blast_matches = {}     # some BLAST match information
        self.blast_read = False     # flag for if BLAST map was read in
        
        # enter main loop
        self.root.mainloop()

    def make_toolbar_button(self, _text, _command, width=12):
        """Makes toolbar buttons
        _text: string; button text
        _command: method; method to be executed when button is clicked
        width: int; desired width of button
        """
        b = Button(self.toolbar, text=_text, width=width, command=_command,
                   borderwidth=2, relief=RAISED)
        b.pack(side=LEFT, padx=5, pady=5)
        return(b)
    
    def make_radiobutton(self, _text, _value):
        """Makes a radio button for option bar'
        _text: string; button text
        _value: int; radiobutton value"""
        rb = Radiobutton(self.rb_frame, text = _text, variable = self.radio_var,
                         value = _value, width = 12, borderwidth = 2)
        rb.pack(side=LEFT, padx=5, pady=5)
        return(rb)
        
    def make_checkbutton(self, _text, var):
        """Makes a checkbutton for the option bar
        _text: string; checkbutton text
        var: int; variable value"""
        cb = Checkbutton(self.cb_frame, text = _text, width = 13, variable = var)
        cb.pack(side=RIGHT, padx=5, pady=5)
        return(cb)    

    def _parse_accessions(self):
        """Helper function to parse an accessions from clipboard."""
        acc = []
        for line in self.accessions:
            line = line.strip().split()[0]
            if not line or line.lower().startswith('accession') or line.lower() == 'acc':
                continue
            if line.endswith('_family'):
                line = line.replace('_family', '')
            acc.append(line)
        return acc
    
    # toolbar button functions

    def get_accessions(self):
        """Gets column of accessions from the clipboard."""
        self.accessions = self.root.clipboard_get()
        self.accessions = self.accessions.splitlines()
        if len(self.accessions) == 0:
            self.clear_screen()
            self.accessions = []
            self.acc_read = False
            self.text.insert("1.0", 'WARNING: Clipboard was empty!')
            self.status.set("%s", "%s accessions read from clipboard" % len(self.accessions))
            return

        # parse the accessions
        self.accessions = pd.DataFrame({'Accession': self._parse_accessions()})
        self.acc_read = True

        # echo the accessions to the screen
        self.echo_dataframe(self.accessions)
        self.status.set("%s", "%s accessions read from clipboard" % len(self.accessions))
        return

    def select_dat_file(self):
        """Get UniProt flat format text file (DAT file)."""        
        ext_list = [('GZip files', '*.gz'), ('DAT files', '*.dat')]
        message = 'Select a UniProt DAT file'
        self.dat_file = get_file(self.default, ext_list, message)   # self.root is the root window

    def parse_dat_file(self):
        """Gets UniProt DB file and makes accession maps."""
        # browse to DAT file
        self.select_dat_file()
        self.status.set("%s", "parsing DAT file or reloading")
        
        # look for pickled annotation dictionary and reload if it exists
        if os.path.exists(self.dat_file + '.pk'):
            read_pk = True
            pickled_anno = pickle.load(open(self.dat_file + '.pk', 'rb'))
            try:
                if (pickled_anno.dat_file != self.dat_file or 
                    pickled_anno.dat_date != os.path.getctime(self.dat_file)):
                    read_pk = False
            except AttributeError:
                read_pk = False
            if read_pk:
                self._annotate_dict = pickled_anno.annotate_dict
                count = int(len(self._annotate_dict) / 3)
        if not read_pk:
            count, self._annotate_dict = self._process_dat_records() # parse DAT file
            
            # save the parsed file results for next time
            pickled_anno = AnnotationPickle(self.dat_file, os.path.getctime(self.dat_file), 
                                            self._annotate_dict)
            pickle.dump(pickled_anno, open(self.dat_file + '.pk', 'wb'), protocol = 2)    # protocol 2 to make compatible with python 2

        print('DB count: %s, dict size: %s' % (count, len(self._annotate_dict)))
        print("DONE processing annotations")
        self.status.set("%s", 'DB count: %s, dict size: %s' % (count, len(self._annotate_dict)))
        self.dat_read = True
##        writeFile(self)   # optionally write annotations to separate files by category
##        self.acc_mapping()  # lookup the annotations for the accessions
##        self.status.set("%s", "%s protein annotation records parsed" % len(self.annotations))
        
    # DAT file processing
    def _process_dat_records(self):
        """Parses all records in a gzipped DAT file."""
        buff = []
        count = 0
        dat_dict = {}
        try:
            for line in gzip.open(self.dat_file, 'rt'):
                line = line.rstrip()
                if line == '//':
                    annotations = Annotations()
                    annotations.parse_record(buff)
                    dat_dict[annotations.identifier] = annotations
                    dat_dict[annotations.accession] = annotations
                    dat_dict[annotations.fasta_accession] = annotations
                    buff = []
                    count += 1
                else:
                    buff.append(line)
        except ValueError:
            for line in gzip.open(self.dat_file):
                line = line.rstrip()
                if line == '//':
                    annotations = Annotations()
                    annotations.parse_record(buff)
                    dat_dict[annotations.identifier] = annotations
                    dat_dict[annotations.accession] = annotations
                    dat_dict[annotations.fasta_accession] = annotations
                    buff = []
                    count += 1
                else:
                    buff.append(line)
        except OSError:
            for line in open(self.dat_file):
                line = line.rstrip()
                if line == '//':
                    annotations = Annotations()
                    annotations.parse_record(buff)
                    dat_dict[annotations.identifier] = annotations
                    dat_dict[annotations.accession] = annotations
                    dat_dict[annotations.fasta_accession] = annotations
                    buff = []
                    count += 1
                else:
                    buff.append(line)
            
        return count, dat_dict
               
    def acc_mapping(self):
        """acc: a list of accessions
        db_dict: DAT file object"""
        # save annotations in a list (should be matched to self.accessions)
        self.annotations = []
        fail_count = 0
        for acc in self.accessions.iloc[:, 0]:
            if acc in self.blast_map:
                acc = self.blast_map[acc]   # work with ortholog accession if it exists
            if acc in self._annotate_dict:
                self.annotations.append(self._annotate_dict[acc])
            else:
                acc_parts = acc.split('|')
                if len(acc_parts) == 3:     # this is assuming UniProt format
                    if acc_parts[2] in self._annotate_dict:
                        temp_acc = acc_parts[2]
                    else:
                        temp_acc = acc_parts[1]
                if temp_acc in self._annotate_dict:
                    self.annotations.append(self._annotate_dict[temp_acc])
                else:
                    print('failed lookup:', acc)
                    fail_count += 1
                    # need a blank annotation object
                    self.annotations.append(Annotations())
                     
        print('\nDone fetching annotations for %s accessions' % len(self.annotations))
        print('...%d lookups failed' % fail_count)
        self.status.set("%s", "%s protein accessions looked up" % len(self.accessions))
        return self.annotations    
       
    def blast_mapping(self):
        """Use Blast ortholog mapping to model organism's DAT file"""
        # check is accessions have been loaded
        if not self.acc_read:
            self.clear_screen()
            self.print_string('Please load some accessions from the clipboard!')
            return

        # browse to BLAST map TXT file
        ext_list = [('Text files', '*.txt')]
        message = 'Select a BLAST mapping file'
        blast_map_file = get_file(self.default, ext_list, message)
        if not blast_map_file: return   # cancel button response
        
        # read the mapping file
        blast = pd.read_table(blast_map_file, skiprows=5)
        blast = blast.dropna(thresh=4) # this should drop rows after the main table
        
        # drop some columns we do not need
        if 'match_status' in blast.columns:
            keep = ['query_acc', 'hit_acc', 'hit_desc', 'blast_scores', 'match_status']
        else:
            keep = ['query_acc', 'hit_acc', 'hit_desc', 'blast_scores', 'status']   # older BLAST map files

        blast_brief = blast[keep]

        # make the accession mapping dictionary (allow for parsed accessions) 
        for query, hit in zip(blast_brief['query_acc'], blast_brief['hit_acc']):
            self.blast_map[query] = hit
            if len(query.split('|')) == 3:     # UniProt format
                   self.blast_map[query.split('|')][1] = hit
                   self.blast_map[query.split('|')][2] = hit
            if (len(query.split('|')) == 1) and ('.' in query): # NCBI and Ensembl format
                self.blast_map[query.split('.')[0]] = hit
            
            
        # save some of the BLAST results in a dictionary keyed by query_acc
        blast_matches = {}
        for row_tuple in blast_brief.iterrows():
            row = row_tuple[1]
            acc = row['query_acc']
            self.blast_matches[row['query_acc']] = [str(x) for x in row]
            if len(acc.split('|')) == 3:     # UniProt format
                   self.blast_matches[acc.split('|')][1] = [str(x) for x in row]
                   self.blast_matches[acc.split('|')][2] = [str(x) for x in row]
            if (len(acc.split('|')) == 1) and ('.' in acc): # NCBI and Ensembl format
                self.blast_matches[acc.split('.')[0]] = [str(x) for x in row]
                
        # organize by accession list (if loaded from clipboard)
        if self.acc_read:
            keys = list(self.accessions['Accession'])
        else:
            keys = self.blast_matches.keys()

        # create a pandas dataframe for the BLAST data
        rows = []
        for key in keys:
            try:
                rows.append([key] + self.blast_matches[key])
            except KeyError:
                rows.append([key, 'NA', 'NA', 'NA', 'NA', 'NA'])

        # make a pandas dataframe
        self.blast_table = pd.DataFrame(rows, columns = ['Index', 'query_acc', 'hit_acc', 'hit_desc', 'blast_scores', 'match_status'])

        # write BLAST results to screen and cliboard
        self.echo_dataframe(self.blast_table)
        self.root.clipboard_clear()
        self.root.clipboard_append(self.blast_table.to_csv(sep='\t', line_terminator='\r', index=False))
        self.status.set("%s", "%s BLAST mappings read in (echoed to screen and clipboard)" % len(self.blast_map))
        self.blast_read = True
        return    
    
    def clear_screen(self):
        """Clears the window."""
        self.text.delete("1.0", END)
        self.status.set("%s", "Screen Cleared")
            
    def clear_data(self):
        """Clears the window contents and clipboard."""
        self.text.delete("1.0", END)
        self.root.clipboard_clear()
        self.root.clipboard_append("")
        self.dat_file = None
        self.dat_read = False
        self.accessions = []
        self.acc_read = False
        self.blast_read = False
        self.status.set("%s", "Data and screen cleared")
    
    def print_help(self):
        """Prints some help info to the window."""
        self.clear_screen()
        help_text = \
"""Program add_protein_annotations.py, version 1.0

Takes a list of protein accessions from a PAW results file (typically Excel) for human,
mouse, or arabidopsis and adds additional annotations from UniProt (Swiss-Prot) DAT files.

"Get accessions" => reads UniProt accessions from the clipboard and echos them to the screen.

"Parse DAT file" => parses a selected UniProt DAT file and extracts annotations of interest.

"Blast mapping" => uses ortholog mapping information from BLAST to bootstrap proteins
from non-model organisms to human, mouse, or arabidopsis orthologs.

"Add annotations" => maps the UniPort accessions from the clipboard to the
annotations from the UniProt DAT file. Writes results to screen and clipboard.

"Reset" => clears data and the screen text window.

"Help" => prints this information.

"Quit" => ends the application.

Written by Kyra Patton, OHSU, 2016.

(and Phil Wilmarth, OHSU, 2016-2019)"""
        self.text.insert("1.0", help_text)
        self.status.set("%s", "Help Text")

    def echo_dataframe(self, frame):
        """Echoes the data from the clipboard to the window."""
        self.clear_screen()
        self.text.insert(CURRENT, frame.to_string(index=False))
        self.text.insert(CURRENT, '\n')

    def print_string(self, string):
        self.text.insert(CURRENT, string)
        self.text.insert(CURRENT, '\n')        

    def add_annotations(self):
        """Prints annotations to the window."""
        # make sure we have accessions and have parsed a DAT file
        return_flag = False
        self.clear_screen()
        if not self.acc_read:
            self.print_string('Please load some accessions from the clipboard!')
            return_flag = True
        if not self.dat_read:
            self.print_string('Please parse a DAT file!')
            return_flag = True
        if return_flag:
            self.status.set("%s", "Annotation lookup failed")
            return

        # lookup the annotations for the accessions
        self.acc_mapping()
        print("%s protein annotation records parsed" % len(self.annotations))
        self.status.set("%s", "%s protein annotation records parsed" % len(self.annotations))

        # format the annotation table (with or without the BLAST mapping info)
        if self.blast_read:
            annot_table = pd.merge(self.blast_table, AnnotationTable(self).table, on='Index')
        else:
            annot_table = AnnotationTable(self).table

        # write annotations to screen and clipboard
        self.echo_dataframe(annot_table)
        self.root.clipboard_clear()
        self.root.clipboard_append(annot_table.to_csv(sep='\t', line_terminator='\r', index=False))
        self.status.set("%s", "Annotations shown above and written to clipboard")
        print('Annotations added for %s proteins' % len(self.accessions))
        
    def quit_me(self):
        """Quits the application."""
        self.status.set("%s", "Bye")
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.quit()
        return

class AnnotationTable:
    """Collects all of the annotations into dataframes."""
    def __init__(self, parent):
        """parent: calling object"""
        self.parent = parent                        # pointer to calling object
        self.dat_file = self.parent.dat_file        # Swiss-Prot DAT file
        self.accessions = self.parent.accessions    # accessions to be annotated
        self.annotations = self.parent.annotations  # annotations from DAT file
        self.reports_folder = None                  # folder to write rports to
        self.default = self.parent.default          # set a default for dialog boxes
        
        self.table = None           # final table for display and export to clipboard
        self.basic_table = None     # table for the general annotations
        self.kw_table = None        # table for the key words
        self.mgi_table = None       # table for MGI info if mouse
        self.pw_table = None        # pathways table
        self.go_table = None        # GO terms table
        
        self.make_main_table()
        if self.parent.radio_var.get() == 2:
            self.make_mgi_table()
        if self.parent.kw_var.get() == 1:
            self.make_kw_table(self.accessions)
        if self.parent.pw_var.get() == 1:
            self.make_pw_table()
        if self.parent.go_var.get() == 1:
            self.make_go_table()
        self.concatenate()
        return
    
    def make_main_table(self):
        """Creates a dataframe of all the non-optional annotations."""
        keys = ['Index', 'Primary Protein Name', 'Alternative Protein Names', 'Identifier', 'Accession',
                'Other Accessions', 'UniProt Gene Name', 'Other Gene Synonyms', 
                'Species Name', 'Taxonomy Number', 'Key Words']
        df_dict = {k: [] for k in keys}
        for acc in self.accessions.iloc[:, 0]:
            df_dict['Index'].append(acc)
        for anno in self.annotations:
            df_dict['Primary Protein Name'].append(anno.name)
            df_dict['Alternative Protein Names'].append(anno.other_names)
            df_dict['Accession'].append(anno.accession)
            df_dict['Identifier'].append(anno.identifier)
            df_dict['Other Accessions'].append(anno.other_accessions)
            if anno.gene:
                df_dict['UniProt Gene Name'].append('="%s"' % anno.gene)
            else:
                df_dict['UniProt Gene Name'].append('na')
            df_dict['Other Gene Synonyms'].append(anno.other_genes)
            df_dict['Species Name'].append(anno.os)
            df_dict['Taxonomy Number'].append(anno.ox)
            df_dict['Key Words'].append(anno.keywords)
            
            #convert annotations that are type list to semicolon separated strings
            for k in keys:
                if type(df_dict[k][-1]) == list:
                    df_dict[k][-1] = '; '.join(df_dict[k][-1]) # use join method

        for k in keys:
            self.basic_table = pd.DataFrame.from_dict(df_dict)
        self.basic_table = self.basic_table[keys]   # put the table columns in the order in keys
        self.basic_table['UniProt Link'] = self.basic_table['Accession'].apply(self.add_uniprot_hyperlinks)
        if self.parent.kw_var.get() == 1:
            self.kw_table = pd.DataFrame(df_dict['Key Words'], columns=['Key Words'])
        self.basic_table = self.basic_table.fillna('na')
        self.basic_table[self.basic_table == ''] = 'na'
        return
            
    def add_uniprot_hyperlinks(self, acc):
        """Adds hyperlinks to UniProt website."""
        if not acc:
            return None
        else:
            return ('=hyperlink("http://www.uniprot.org/uniprot/' + acc + '", "' + acc + '")') 
            
    def make_mgi_table(self):
        """Makes table of MGI annotations."""
        mgi_dict = {}
        mgi_dict['MGI Accession'] = []
        mgi_dict['MGI Gene Name'] = []
        mgi_dict['MGI Link'] = []
        for anno in self.annotations:
            mgi_dict['MGI Accession'].append(anno.mgi_acc)
            if anno.mgi_acc:
                mgi_dict['MGI Link'].append('=HYPERLINK("http://www.informatics.jax.org/marker/' + anno.mgi_acc + '", "' + anno.mgi_acc + '")')
            else:
                mgi_dict['MGI Link'].append('na')
            if anno.mgi_gene:
                mgi_dict['MGI Gene Name'].append('="%s"' % anno.mgi_gene)
            else:
                mgi_dict['MGI Gene Name'].append('na')
        self.mgi_table = pd.DataFrame.from_dict(mgi_dict)
        self.mgi_table = self.mgi_table.fillna('na')
        self.mgi_table[self.mgi_table == ''] = 'na'
        return
            
    def make_kw_table(self, accessions):
        """Makes table of keyword annotations."""
        self.kw = KeyWords()
        try:
            self.kw.parse_file(os.path.join(os.path.dirname(self.dat_file), 'keywlist.txt'))  # parse keyword definitions file
        except:
            """Need to browse to file if not found!"""
            print('\nWARNING: key word list definition file not found\n')
            return
            
        # analyze the keyword frequencies and associated proteins
        if self.parent.sf_var.get() == 1:
            self.analyze_keywords()
            
        # put the keywords into their 10 categories
        keywords_by_category = []
        for anno in self.annotations:
            keywords_by_category.append(self.kw.put_keywords_in_categories(anno.keywords))
        cat_table = pd.DataFrame(keywords_by_category, columns=['KW: '+cat for cat in self.kw.categories])
        self.kw_table = pd.concat([self.kw_table, cat_table], axis=1)        
        self.kw_table = self.kw_table.fillna('na')
        self.kw_table[self.kw_table == ''] = 'na'
        return
    
    def analyze_keywords(self):
        """Frequency anaylsis of key words and inverse mapping to proteins."""
        keyword_freq = {}
        for ident, keywords in zip(self.basic_table['Identifier'], self.basic_table['Key Words']):
            if ident == 'na':
                continue
            for keyword in keywords.split('; '):
                if keyword == 'na':
                    continue
                if keyword in keyword_freq:
                    keyword_freq[keyword].append(ident)
                else:
                    keyword_freq[keyword] = [ident]
        keyword_items = keyword_freq.items()
        keyword_rows = [(x, len(y), '; '.join(y)) for (x, y) in keyword_items]
        keyword_rows = sorted(keyword_rows, key=lambda x: x[1], reverse=True)
        
        # write frequency report to KW report file
        if not self.reports_folder:
            self.reports_folder = get_folder(self.default, 'Select a folder for reports')
            if not self.reports_folder:
                return
        report = open(os.path.join(self.reports_folder, 'keyword_report.txt'), 'w')
        print('KeyWord Report generated on:', time.ctime(), file=report)
        print('Total number of key words was:', len(keyword_rows), '\n', file=report)
        columns = ['Keyword', 'Category', 'Description', 'Synonyms', 'Frequency', 'Proteins']
        print('\t'.join(columns), file=report)
        for kw, freq, prots in keyword_rows:
            kw_obj = self.kw.keywords[kw]
            row = [kw, kw_obj.category, kw_obj.definition, '; '.join(kw_obj.synonyms), str(freq), prots]
            if freq > 1:
                print('\t'.join(row), file=report)
        report.close()
        return
    
    def make_pw_table(self):
        "Make a dataframe of the pathway annotations."
        pw_dict = {}
        pw_dict['CC Pathway'] = []
        pw_dict['Reactome Pathway'] = []
        for anno in self.annotations:
            pw_dict['CC Pathway'].append(anno.pathway.cc_string)
            pw_dict['Reactome Pathway'].append(anno.pathway.react_string)
        self.pw_table = pd.DataFrame.from_dict(pw_dict)
        self.pw_table = self.pw_table.fillna('na')
        self.pw_table[self.pw_table == ''] = 'na'
        
        # analyze and write report of Reactome pathways
        if self.parent.sf_var.get() == 1:
            self.analyze_pathways()
        return
        
    def analyze_pathways(self):
        """Frequency anaylsis of Reactome pathways and inverse mapping to proteins."""
        pathway_freq = {}
        for ident, pathways in zip(self.basic_table['Identifier'], self.pw_table['Reactome Pathway']):
            if ident == 'na':
                continue
            for pathway in pathways.split('; '):
                if pathway == 'na':
                    continue
                if pathway in pathway_freq:
                    pathway_freq[pathway].append(ident)
                else:
                    pathway_freq[pathway] = [ident]
        pathway_items = pathway_freq.items()
        pathway_rows = [(x, len(y), '; '.join(y)) for (x, y) in pathway_items]
        pathway_rows = sorted(pathway_rows, key=lambda x: x[1], reverse=True)
        
        # write frequency report to pathway report file
        if not self.reports_folder:
            self.reports_folder = get_folder(self.default, 'Select a folder for reports')
            if not self.reports_folder:
                return
        report = open(os.path.join(self.reports_folder, 'pathway_report.txt'), 'w')
        print('Pathway Report generated on:', time.ctime(), file=report)
        print('Total number of pathways was:', len(pathway_rows), '\n', file=report)
        columns = ['Identifier', 'Description', 'Link', 'Frequency', 'Proteins']
        print('\t'.join(columns), file=report)
        for pw, freq, prots in pathway_rows:
            desc, ident = pw.split('{')[0], pw.split('{')[1][:-1]
            link = '=hyperlink("http://www.reactome.org/content/detail/' + ident + '", "' + ident + '")'
            row = [ident, desc, link, str(freq), prots]
            if freq > 1:
                print('\t'.join(row), file=report)
        report.close()
        return
    
    def make_go_table(self):
        """Make a dataframe of the go terms."""
        go_dict = {}
        for key in ['GO: Biological Process', 'GO: Cellular Component', 'GO: Molecular Function']:
            go_dict[key] = []
        for anno in self.annotations:
            # not all protein records have GO terms, so need to handle missing data 
            try:
                go_dict['GO: Biological Process'].append(anno.go.biological_process)
                go_dict['GO: Cellular Component'].append(anno.go.cellular_component)
                go_dict['GO: Molecular Function'].append(anno.go.molecular_function)
            except AttributeError:
                go_dict['GO: Biological Process'].append('na')
                go_dict['GO: Cellular Component'].append('na')
                go_dict['GO: Molecular Function'].append('na')
        self.go_table = pd.DataFrame.from_dict(go_dict)
        self.go_table = self.go_table.fillna('na')
        self.go_table[self.go_table == ''] = 'na'
        
        # analyze and write report of GO terms
        if self.parent.sf_var.get() == 1:
            self.analyze_GOTerms()
        return
        
    def analyze_GOTerms(self):
        """Frequency anaylsis of GO terms and inverse mapping to proteins."""
        # write frequency report to GO Terms report file
        if not self.reports_folder:
            self.reports_folder = get_folder(self.default, 'Select a folder for reports')
            if not self.reports_folder:
                return
        report = open(os.path.join(self.reports_folder, 'GOTerms_report.txt'), 'w')
        print('GO Term Report generated on:', time.ctime(), file=report)

        for go_category in ['GO: Biological Process', 'GO: Cellular Component', 'GO: Molecular Function']:
            goterm_freq = {}
            for ident, goterms in zip(self.basic_table['Identifier'], self.go_table[go_category]):
                if ident == 'na':
                    continue
                for goterm in goterms.split('; '):
                    if goterm == 'na':
                        continue
                    if goterm in goterm_freq:
                        goterm_freq[goterm].append(ident)
                    else:
                        goterm_freq[goterm] = [ident]
            goterm_items = goterm_freq.items()
            goterm_rows = [(x, len(y), '; '.join(y)) for (x, y) in goterm_items]
            goterm_rows = sorted(goterm_rows, key=lambda x: x[1], reverse=True)
            print('\n\nTotal number of %s terms was: %s\n' % (go_category, len(goterm_rows)), file=report)
            columns = ['Identifier', 'Description', 'Link', 'Frequency', 'Proteins']
            print('\t'.join(columns), file=report)
            for go, freq, prots in goterm_rows:
                desc, acc = go.split('{')[0], go.split('{')[1][:-1]
                link = '=hyperlink("http://amigo.geneontology.org/amigo/term/' + acc + '", "' + acc + '")'
                row = [acc, desc, link, str(freq), prots]
                if freq > 1:
                    print('\t'.join(row), file=report)
        report.close()
        return
        
    def concatenate(self):
        """Merges all the annotation tables together."""
        # drop any unwanted columns before concatenating
        self.basic_table = self.basic_table.drop('Key Words', axis=1)    
        frames = [self.basic_table, self.mgi_table, self.kw_table, self.go_table, self.pw_table]
        self.table = pd.concat(frames, axis=1)
        self.table = self.table.set_index('Index')
        return

# MAIN program starts here
annotator = ProteinAnnotator()
# end when user hits quit button or closes window
