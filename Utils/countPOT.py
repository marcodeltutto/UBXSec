#!/usr/bin/env python

# usage: python countPOT.py -f /path/to/ubxsec_output.root

import argparse
parser = argparse.ArgumentParser(description='Counts POT.')
parser.add_argument('-f', action='store', dest='infile',
                    help='File name (multiple files using * are ok.')

args = parser.parse_args()
#print 'Input', args.infile

import os,sys,string, time
import ROOT
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain
from array import array
from glob import glob
import subprocess


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

file = open("runsubrun_list.txt","w") 

# Getting POT number
NominalPOT = 6.0e20
TotalPOT = 0
for jentry in xrange( entriesPOT ):
  entry = chainPOT.GetEntry( jentry )
  #print '{:e}'.format(float(chainPOT.pot))
  TotalPOT += chainPOT.pot
  #print "event:", chainPOT.run
  #print "event:", chainPOT.subrun
  file.write("%d" % chainPOT.run)
  file.write(" ")
  file.write("%d" % chainPOT.subrun)
  file.write("\n")
#print "Accumulated POT: ", '{:e}'.format(float(TotalPOT))
#print "Nominal POT:     ", '{:e}'.format(float(NominalPOT))

file.close()


os.system('/uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list runsubrun_list.txt')
os.system('rm runsubrun_list.txt')



#raw_input("Please press enter to exit.")


