import os
import shutil
import subprocess
import copy



class ParametricJobsMP:
    def __init__(self, name):
        self.name = name
        self.exe = ""
        self.libs = []
        self.incl = []
        self.creatingDict = False
        self.src = []
        self.parameters = ""
        self.ce = "llrmpicream"
        self.nJobs = 1
        self.inputFiles = []
        self.parameterFileBaseName = "dummy.$jobid.par"
        self.parameterFiles = []
        self.outputDir = ""
        self.currentDir = os.getcwd()
        self.args = []
        self.cmsswRel   = os.environ["CMSSW_VERSION"]
        self.scram_arch = os.environ["SCRAM_ARCH"]

    def checkParams(self):
        if not os.path.isfile(self.exe):
            raise StandardError("ERROR: ParametricJobsMP: Non existing executable "+self.exe)
        for lib in self.libs:
            if not os.path.isfile(lib):
                raise StandardError("ERROR: ParametricJobsMP: Non existing library "+lib)
        for incl in self.incl:
            if not os.path.isfile(incl):
                raise StandardError("ERROR: ParametricJobsMP: Non existing include "+incl)
        for src in self.src:
            if not os.path.isfile(src):
                raise StandardError("ERROR: ParametricJobsMP: Non existing source "+src)
        for inp in self.inputFiles:
            if not os.path.isfile(self.outputDir+"/jobs/"+inp):
                raise StandardError("ERROR: ParametricJobsMP: Non existing input file "+inp+" in jobs/ directory")
        for par in self.parameterFiles:
            if not os.path.isfile(self.outputDir+"/jobs/"+par):
                raise StandardError("ERROR: ParametricJobsMP: Non existing parameter file "+par+" in jobs/ directory")
        if len(self.parameterFiles)!=self.nJobs:
            raise StandardError("ERROR: ParametricJobsMP: Number of parameter files inconsistent with the number of jobs")
        for i in range(1,self.nJobs+1):
            baseName = copy.copy(self.parameterFileBaseName)
            baseName = baseName.replace("$jobid", str(i))
            if baseName!=self.parameterFiles[i-1]:
                raise StandardError("ERROR: ParametricJobsMP: parameter file name "+self.parameterFiles[i-1]+" for job "+str(i)+" is not correctly formed. Should be "+baseName)


    def prepareArea(self):
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)
        if not os.path.isdir(self.outputDir+"/jobs"):
            os.makedirs(self.outputDir+"/jobs")
        if not os.path.isdir(self.outputDir+"/obj"):
            os.makedirs(self.outputDir+"/obj")
        if not os.path.isdir(self.outputDir+"/include"):
            os.makedirs(self.outputDir+"/include")
        if not os.path.isdir(self.outputDir+"/src"):
            os.makedirs(self.outputDir+"/src")
        shutil.copy(self.exe, self.outputDir)
        for lib in self.libs:
            shutil.copy(lib, self.outputDir+"/obj/")
        for incl in self.incl:
            shutil.copy(incl, self.outputDir+"/include/")
        for src in self.src:
            shutil.copy(src, self.outputDir+"/src/")



    def createScript(self):
        with open(self.outputDir+"/"+self.name+".sub", 'w') as script:
            print >>script, "#! /bin/sh"
            print >>script, "uname -a"
            print >>script, "\necho '\n-- Running parametric job $jobid'"
            print >>script, "\necho '\n-- Setting up cms environment'"
            print >>script, "export SCRAM_ARCH={}".format(self.scram_arch) 
            #print >>script, "source /opt/exp_soft/cms/cmsset_default.sh"
            print >>script, "source /cvmfs/cms.cern.ch/cmsset_default.sh "
            print >>script, "cmsrel {}".format(self.cmsswRel)
            print >>script, "cd {}/src".format(self.cmsswRel)
            print >>script, "cmsenv"
            print >>script, "cd -"
            print >>script, "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/obj/"
            if self.creatingDict:
                print >>script, "\necho '\n-- Creating ROOT dictionary'"
                print >>script, "cd include"
                print >>script, "root -b -q libDictionary.C+"
                print >>script, "mkdir ../obj"
                print >>script, "mv libDictionary_C.so ../obj/"
                print >>script, "cd -"
            print >>script, "\necho '\n-- Executing job'"
            print >>script, "chmod +x "+os.path.basename(self.exe)
            print >>script, "./"+os.path.basename(self.exe), "jobs/"+self.parameterFileBaseName
            print >>script, "rm -rf {}".format(self.cmsswRel)
            print >>script, "rm -rf obj"
            print >>script, "\necho '\n-- Scratch disk status'"
            print >>script, "ls -ltr"
            print >>script, "du -hs ."

    def prepareCommand(self):
        self.args.append("-submit")
        self.args.append("-ce")
        self.args.append("llrmpicream.in2p3.fr")
        self.args.append("-mpi")
        self.args.append("12345")
        self.args.append("-queue")
        self.args.append("gridq")
        self.args.append("-cmd")
        self.args.append("'./"+self.name+".sub'")
        self.args.append("-exe")
        self.args.append("'"+self.name+".sub'")
        self.args.append("-files")
        fileList = "'"
        fileList += os.path.basename(self.exe)+","
        for lib in self.libs:
            fileList += "obj/"+os.path.basename(lib)+","
        for incl in self.incl:
            fileList += "include/"+os.path.basename(incl)+","
        for src in self.src:
            fileList += "src/"+os.path.basename(src)+","
        fileList += self.name+".sub,"
        fileList += "jobs'"
        self.args.append(fileList)
        self.args.append("-ids")
        self.args.append("1-"+str(self.nJobs))

    def launch(self, simulate=False):
        os.chdir(self.outputDir)
        command = "/opt/exp_soft/vo.llr.in2p3.fr/tools/parametricJobs"
        for arg in self.args:
            command += " "+arg
        print command
        if not simulate:
            subprocess.call(["/opt/exp_soft/vo.llr.in2p3.fr/tools/parametricJobs"] + self.args)
        os.chdir(self.currentDir)



if __name__=="__main__":
    job = ParametricJobsMP("testJob")


