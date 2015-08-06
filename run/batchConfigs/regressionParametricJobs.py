#! /usr/bin/env python

import sys
import os
import shutil
import optparse
import subprocess

sys.path.append(os.getcwd())

parser = optparse.OptionParser()
parser.add_option('--launch', action='store', type='string', dest='launch', help='Launch regressions')
parser.add_option('--status', action='store', type='string', dest='status', help='Chech status of regressions')
parser.add_option('--retrieve', action='store', type='string', dest='retrieve', help='Get output of regressions')
parser.add_option('--check', action='store', type='string', dest='check', help='Check output of regressions')
(options, args) = parser.parse_args()


if options.launch:
    parameterFile = options.launch
    b = __import__(os.path.splitext(parameterFile)[0]).batch
    b.parameterFile = parameterFile
    b.initialize()
    b.launch(simulate=False)
elif options.status:
    parameterFile = options.status
    b = __import__(os.path.splitext(parameterFile)[0]).batch
    b.initialize()
    print "Status version "+b.findLastVersion()
    outputDir = b.outputDirectory+"/"+b.findLastVersion()
    currentDir = os.getcwd()
    os.chdir(outputDir)
    subprocess.call(["/opt/exp_soft/vo.llr.in2p3.fr/tools/parametricJobs", "-status"])
    os.chdir(currentDir)
elif options.retrieve:
    parameterFile = options.retrieve
    b = __import__(os.path.splitext(parameterFile)[0]).batch
    b.initialize()
    print "Retrieve version "+b.findLastVersion()
    outputDir = b.outputDirectory+"/"+b.findLastVersion()
    currentDir = os.getcwd()
    nJobs = len(b.regressions.items())
    os.chdir(outputDir)
    subprocess.call(["/opt/exp_soft/vo.llr.in2p3.fr/tools/parametricJobs", "-getoutput", "-ids", "1-"+str(nJobs)])
    shutil.rmtree("include")
    shutil.rmtree("src")
    shutil.rmtree("obj")
    os.remove("regression.exe")
    os.chdir(currentDir)
elif options.check:
    parameterFile = options.check
    b = __import__(os.path.splitext(parameterFile)[0]).batch
    b.initialize()
    b.checkOutput()
else:
    parser.print_help()
