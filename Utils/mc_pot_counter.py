#!/usr/bin/env python

# usage: python mc_pot_counter.py -f /path/to/ubxsec_output.root

import argparse
parser = argparse.ArgumentParser(description='Counts POT.')
parser.add_argument('-f', action='store', dest='infile',
                    help='File name (multiple files using * are ok.')

args = parser.parse_args()

import os,sys,string, time
import ROOT
from math import *
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain
#from ROOT import *
from array import array
from glob import glob

# Opening root file
fname = glob(args.infile)
print "Input file(s): ", fname

# Creating TChain
chainPOT = TChain("UBXSec/pottree")
for f in fname: 
  chainPOT.Add(f)
           
# Printing the number of entries
entriesPOT = chainPOT.GetEntries()
print "Number of entries in the POT tree: ", entriesPOT

# Getting POT number
NominalPOT = 6.0e20
TotalPOT = 0
for jentry in xrange( entriesPOT ):
  entry = chainPOT.GetEntry( jentry )
  #print '{:e}'.format(float(chainPOT.pot))
  TotalPOT += chainPOT.pot 
print "Accumulated POT: ", '{:e}'.format(float(TotalPOT))


raw_input("Please press enter to exit.")


