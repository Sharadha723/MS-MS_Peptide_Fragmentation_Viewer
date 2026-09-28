# Write a program to display peptide fragmentation spectra from an mzXML file.
# The program will take an mzXML file, a scan number, and a peptide sequence as input.
# The peptide's b-ion and y-ion m/z values should be computed, and peaks matching these m/z values annotated with appropriate labels.
# The output figure/plot should aid the user in determining whether or not the peptide is a good match to the spectrum.

import sys
import gzip
from base64 import b64decode
import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from array import array 

# Error handling if not enough input given
if len(sys.argv)<4:
    print("Enter file name, scan number and peptide sequence", file=sys.stderr)
    sys.exit(1)

# unzipping file
mzxml_file = gzip.open(sys.argv[1],'rt')

# Make sure scan number provided is an integer
try:
    scan_num = int(sys.argv[2])
except ValueError:
    print("Scan number must be an integer")

# Define amino acid mass dictionary 
aa_mw = {
    'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07, 
    'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08, 
    'M': 131.04,'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10, 
    'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06
}

# Make sure only valid amino acid sequence is entered
peptide_seq = sys.argv[3]

for aa in peptide_seq.upper():
    if aa not in aa_mw:
        print("Invalid amino acid:", aa)
        sys.exit(1)

# Function to compute b-ion and y-ion m/z values 
def compute_b_y_ions(peptide_seq):
    b_ions = []
    y_ions = []
    peptide_length = len(peptide_seq)
    b_mass = 0
    y_mass = 0

    for aa_b, aa_y in zip(peptide_seq, reversed(peptide_seq)):
        b_mass += aa_mw[aa_b]
        b_ions.append(b_mass+1)

        y_mass += aa_mw[aa_y]
        y_ions.append(y_mass+19)

    # Print final b-ions and y-ions after calculation
    # print("Final b-ions:", b_ions)
    # print("Final y-ions:", y_ions)
    return b_ions, y_ions    
    

# Function to extract mz and intensity data from mzXML file
def parse_mzxml(mzxml_file, scan_number):
    ns = ''  
    scan_found = False
    mzs = None
    ints = None
    try:
        for event, elem in ET.iterparse(mzxml_file):
            if ns == '':
                p = elem.tag.find('}')
                if p>= 0:
                   ns = elem.tag[:(p+1)]

            if event == "end" and elem.tag == ns +'scan':
                if int(elem.attrib.get("num")) == scan_number:
                    scan_found = True
                    peaks_element = elem.find(ns + 'peaks')
                    if peaks_element.text:
                        peaks = array('f', b64decode(peaks_element.text))
                        if sys.byteorder != 'big':
                            peaks.byteswap()
                        mzs = peaks[::2]
                        ints = peaks[1::2]
                elem.clear()
        if not scan_found:
             print("Error: Scan number not found")
        mzxml_file.close()
    except ET.ParseError:
        print("Error parsing XML")
    return mzs, ints

# Funtion to find best peaks
def match_best_peak(mzs, ints, targets, tolerance = 0.05):
    matches = {}  #store best matches
    int_threshold = 0.05 * max(ints)   #set threshold
    
    for target_mz in targets:
        best_match = None
        best_intensity = 0

        for i, exp_mz in enumerate(mzs):
         # check if exp mz is within tolerance and above intensity threshold
            if abs(exp_mz - target_mz) <= tolerance and ints[i] >= int_threshold: 
                # if current intensity has higher intensity than previous best match, update best match
                if ints[i] > best_intensity:  
                    best_match = (exp_mz, ints[i])
                    best_intensity = ints[i]
        matches[target_mz] = best_match
    return matches

# Function to plot the spectrum and annotate with b/y ions
def plot_spectrum(mz_value, int_value, b_ions, y_ions, tolerance = 0.05):    
    if mz_value is None or int_value is None:
        print("No m/z or intensity value found")
        return
    
    # add labels for legends, but only once
    b_label_added = False
    y_label_added = False    

    unmatched_plot = plt.stem(mz_value, int_value, linefmt = 'gray', markerfmt = '', basefmt = '', 
                              bottom = 0, label =  'Unmatched peaks')
    unmatched_plot.baseline.set_visible(False) #Stem Container from matplotlib 

    # plot b_ions with one label
    for i, b_ion in enumerate(b_ions):
        matches = match_best_peak(mz_value, int_value, [b_ion], tolerance)
        best_match = matches.get(b_ion)
        if best_match:
            exp_mz, intensity = best_match
            if not b_label_added:
                b_ions_plot = plt.stem([exp_mz],[intensity],linefmt = 'blue', markerfmt= '',basefmt = '', 
                                       label = 'b-ions')
                b_ions_plot.baseline.set_visible(False)
                b_label_added = True
            else:
                plt.stem([exp_mz],[intensity],linefmt = 'blue', markerfmt = '', basefmt = '')
            plt.annotate("b" + str(i+1),
                         (exp_mz, intensity +10),
                         color = 'blue', 
                         fontsize = 10,
                         ha = 'center',
                         va = 'bottom')

    # plot y_ions with one label
    for i,y_ion in enumerate(y_ions):
        matches = match_best_peak(mz_value, int_value, [y_ion], tolerance)
        best_match = matches.get(y_ion)
        if best_match:
            exp_mz, intensity = best_match
            if not y_label_added:
                y_ions_plot = plt.stem([exp_mz],[intensity],linefmt = 'red',markerfmt = '',basefmt = '', 
                                       label = 'y-ions')
                y_ions_plot.baseline.set_visible(False) 
                y_label_added = True
            else:
                plt.stem([exp_mz],[intensity],linefmt = 'red',markerfmt = '',basefmt = '')
            plt.annotate("y" + str(i+1),
                        (exp_mz, intensity + 10),
                         color = 'red',
                         fontsize = 10,
                         ha = 'center',
                         va = 'bottom')

    # setting y axis limits
    plt.ylim(0, max(int_value) * 1.2) 
    
    # setting axis title and graph title                
    plt.title("Spectrum for Scan:"+str(scan_num)+" Peptide:"+peptide_seq, fontsize = 10)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")

    # adding a legend
    plt.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.show()

# Parse mzxml file and get mz and intensity values
mz_value, int_value = parse_mzxml(mzxml_file, scan_num)

# compute b and y ions
b_ions, y_ions = compute_b_y_ions(peptide_seq)

# plot peptide spectrum
plot_spectrum(mz_value, int_value, b_ions, y_ions, tolerance = 0.05)






