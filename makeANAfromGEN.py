#!usr/bin/python
import os
import optparse 
import sys 
import subprocess
import datetime 
import re
from subprocess import check_output
import ROOT
import numpy as nn
from random import randint

basedir_path = os.path.dirname(os.path.realpath(__file__))
print basedir_path

usage = ""
parser = optparse.OptionParser(usage='\nExample: python %prog -c GenAna_cfg.py -i list.txt -t /tmp/santanas/ --outputDir `pwd`/TestOutput')
parser.add_option("-c","--pythonConfig",action="store",type="string",dest="PYTHONCONFIG",default="")
parser.add_option("-i","--inputGENList",action="store",type="string",dest="INPUTGENLIST",default="")
parser.add_option("-t","--tmpDir",action="store",type="string",dest="TMPDIR",default="")
parser.add_option("--numberOfevents",action="store",type="string",dest="NEVENTS",default="-1")
parser.add_option("--outputDir",action="store",type="string",dest="OUTPUTDIR",default="")

(options, args) = parser.parse_args()
PYTHONCONFIG = options.PYTHONCONFIG
INPUTGENLIST = options.INPUTGENLIST
TMPDIR = options.TMPDIR
NEVENTS = options.NEVENTS
OUTPUTDIR = options.OUTPUTDIR

if not options.PYTHONCONFIG:   
    parser.error('ERROR: Input python config is not given')
if not options.INPUTGENLIST:   
    parser.error('ERROR: Input GEN list is not given')
if not options.TMPDIR:   
    parser.error('ERROR: Tmp directory is not given')
if not options.OUTPUTDIR:   
    parser.error('ERROR: Output directory is not given')

#get list of lhe files
proc = subprocess.Popen(["less %s" % INPUTGENLIST], stdout=subprocess.PIPE, shell=True)
(genfilelist, err) = proc.communicate()
genfilelist = genfilelist.splitlines()
print genfilelist

#create output dir
IsEosDir = False
if ("/eos/" in OUTPUTDIR):
    print "IS EOS DIR"
    IsEosDir = True
else:
    print "IS NOT EOS DIR"
    IsEosDir = False

if IsEosDir:
    print("eos mkdir -p %s" % (OUTPUTDIR) )
    os.system("eos mkdir -p %s" % (OUTPUTDIR) )
else:
    print("mkdir -p %s" % (OUTPUTDIR) )
    os.system("mkdir -p %s" % (OUTPUTDIR) )

# loop over gripacks
for genfile in genfilelist:

    outputfilename = ((genfile.split("/")[-1]).split("."))[0]+"_ANALYSIS.root"
    #print outputfilename

    print("cmsRun %s files=%s output=%s/%s maxEvents=%s" % (PYTHONCONFIG,genfile,TMPDIR,outputfilename,NEVENTS) )
    os.system("cmsRun %s files=%s output=%s/%s maxEvents=%s" % (PYTHONCONFIG,genfile,TMPDIR,outputfilename,NEVENTS) )

    # move output in final directory
    if IsEosDir:
        TMPFILE = ("%s/%s" % (TMPDIR,outputfilename))
        TMPFILE = re.sub("//","/",TMPFILE) 
        OUTPUTFILE = ("%s/%s" % (OUTPUTDIR,outputfilename))
        OUTPUTFILE = re.sub("//","/",OUTPUTFILE) 

        print("eos cp %s %s" % (TMPFILE,OUTPUTFILE))
        os.system("eos cp %s %s" % (TMPFILE,OUTPUTFILE))
        print("rm -f %s" % (TMPFILE))
        os.system("rm -f %s" % (TMPFILE))
    else:
        print("mv %s/%s %s/%s" % (TMPDIR,outputfilename,OUTPUTDIR,outputfilename))
        os.system("mv %s/%s %s/%s" % (TMPDIR,outputfilename,OUTPUTDIR,outputfilename))


