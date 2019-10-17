"""keywlist_sprot-dat_download.py - downloads key word list and Swiss-Prot DAT
   format records (all of the annotations) from UniProt FTP site. Parses out
   the three species of interest (human, mouse, and arabidopsis)

   20191006 - Phil Wilmarth, OHSU

The MIT License (MIT)

Copyright (c) 2019 Phillip A. Wilmarth, OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import sys
import ftplib
import datetime
import gzip
import io
import tkinter
from tkinter import filedialog

def get_folder(default_location, title_string=None):
    """Dialog box to browse to a folder.  Returns folder path.

    Usage: full_folder_name = get_folder(default_location, [title]),
        where "default_location" is a starting folder location,
        "title" is an optional message to list in the dialog box,
        and "full_folder_name" is the complete selected folder name.
    Written by Phil Wilmarth, 2008, 2016
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string and location if not passed
    if title_string is None:   
        title_string = 'Select a folder with desired files/dirs'
    if not default_location:
        default_location = os.getcwd()
    
    # create dialog box for folder selection
    root.update()   # helps make sure dialog box goes away after selection
    full_folder_name = filedialog.askdirectory(parent=root, initialdir=default_location, 
                                               title=title_string, mustexist=True)    
    # return full folder name
    return full_folder_name

def fetch_keywlist(location):
    """fetches keywlist.txt file from UniProt FTP site with error testing and retries.
    Has a hard failure if file cannot be downloaded.
    """
    # login
    ftp = ftplib.FTP('ftp.uniprot.org')
    ftp.login()
    
    # move into documents location
    ftp.cwd('/pub/databases/uniprot/current_release/knowledgebase/complete/docs/')

    # get file contents
    listing = []
    ftp.retrlines('RETR keywlist.txt', listing.append)
    
    with open(os.path.join(location, 'keywlist.txt'), 'wt') as fout:
        for line in listing:
            line = line.rstrip()
            print(line, file=fout)
        print('...key word file was downloaded and saved')

    # done
    ftp.quit()
    return

def fetch_sprot_dat(location):
    """Fetches uniprot_sprot.dat.gz file from UniProt FTP site.
    """
    # login or reconnect
    ftp = ftplib.FTP('ftp.uniprot.org')
    ftp.login() 
    
    # move into documents location
    ftp.cwd('/pub/databases/uniprot/current_release/knowledgebase/complete/')

    # get file contents
    print('...starting download of UniProt Sprot dat file')
    fname = 'uniprot_sprot.dat.gz'
    ftp.retrbinary('RETR {}'.format(fname), open(os.path.join(location, fname), 'wb').write)
    print('...{} was downloaded OK'.format(fname))

    # logout connection and return
    ftp.quit()
    return

def check_buffer(buff, out_lines):
    """Checks for human, mouse, or arabidopsis records
    and appends any found to out_lines.
    """
    species_list = ['3702', '9606', '10090']
    species_count = 0
    for line in buff:
        if line.startswith('OX   '):
            species = line.split('NCBI_TaxID=')[1].split(';')[0]
            if species in species_list:
                out_lines += buff
                return True
            else:
                return False

def parse_sprot_dat(location):
    """Parses sprot.dat file.
    """
    # make sure file exists
    fname = os.path.join(location, 'uniprot_sprot.dat.gz')
    if not os.path.exists(fname):
        print('%s file not found' % fname)
        return

    print('...parsing:', fname)
    buff = ['//']
    out_lines = []
    species_count = 0
    # read the gzipped file (https://pymotw.com/3/gzip/)
    with gzip.open(fname, 'rb') as fin:
        with io.TextIOWrapper(fin, encoding='utf-8') as dec:
            for line in dec:
                if (not line) or line.startswith('//'):
                    if check_buffer(buff, out_lines):
                        species_count += 1
                    buff = ['//']
                else:
                    buff.append(line.rstrip())
    
    print('...there were %d human/mouse/arabidopsis records' % species_count)
    out_lines.append('//')
    return out_lines[1:]    # skip first line ("\\")
    
################################################################################
# get loacation (folder) for downloads
location = get_folder(os.getcwd(), 'Select folder for downloads')
if not location:
    sys.exit()
            
# get key word list contents
fetch_keywlist(location)

# get UniProt full Swiss-Prot records
fetch_sprot_dat(location)

# filter out the human, mouse, or arabidopsis records
out_lines = parse_sprot_dat(location)

# save as gzipped file with appended date stamp (https://pymotw.com/3/gzip/)
today = datetime.date.today().strftime("%Y%m%d")
out_name = os.path.join(location, 'sprot-dat_3702-9606-10090_' + today + '.dat.gz')
print('...writing parsed records')
with gzip.open(out_name, 'wb') as output:
    with io.TextIOWrapper(output, encoding='utf-8') as enc:
        enc.write('\n'.join(out_lines) + '\n')

# remove the full Swiss-Prot DAT downloaded file
os.remove(os.path.join(location, 'uniprot_sprot.dat.gz'))

# fini

