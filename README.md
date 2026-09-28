# Peptide Fragmentation Spectrum

A Python program for visualizing peptide fragmentation spectra from mzXML 
mass spectrometry data.

The program takes an mzXML file, a scan number, and a peptide sequence as
inputs. It calculates theoretical b-ion and y-ion m/z values for the peptide,
matches them to experimental peaks in the selected MS/MS spectrum, and
annotates matching peaks on the resulting spectrum.

## Project Overview

In tandem mass spectrometry (MS/MS), peptides can be fragmented into
characteristic b- and y-ions. Comparing theoretical fragment ion masses with
experimental spectra can help assess whether a proposed peptide sequence is
consistent with the observed spectrum.

This project performs this comparison by:

1. Reading a compressed mzXML file.
2. Selecting a specified scan.
3. Extracting experimental m/z and intensity values.
4. Calculating theoretical b-ion and y-ion m/z values from a peptide sequence.
5. Matching theoretical ions to experimental peaks within a specified m/z
   tolerance.
6. Annotating matched b- and y-ions on the experimental spectrum.
7. Displaying the resulting annotated spectrum.

The example commands in this repository were developed using the
`17mix_test2.mzxml` dataset provided for the course.

## Requirements

Python 3.x

Required Python packages:

- NumPy
- Matplotlib

The program also uses Python standard-library modules including:
- `sys`
- `gzip`
- `base64`
- `xml.etree.ElementTree`
- `array`

