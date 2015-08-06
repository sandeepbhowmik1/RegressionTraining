from ParametricJobsMP import ParametricJobsMP
import copy
import os
import shutil
import stat
import sys
import datetime
import ROOT
from Regression import Regression



class RegressionParametricJobLauncher:
    def __init__(self):
        self.baseDir = "/home/llr/cms/sauvan/RegressionTraining/"
        self.exe = "regression.exe"
        self.libs = []
        self.includeDir = self.baseDir+"/include/"
        self.srcDir = self.baseDir+"/src/"
        self.trainerType = "GBRTrain"
        self.baseName = "BDT"
        self.tree = "regNtuplizer/simpleRootTree"
        self.method = "BDT"
        self.inputFiles = []
        self.tmvaFactoryOptions = ["!V","!Silent","!Color","!DrawProgressBar"]
        self.tmvaTrainingOptions = ["SplitMode=random","!V"]
        self.outputDirectory = "./"
        self.doErrors = False
        self.doCombine = False
        self.commonOptions = []
        self.commonVariablesEB = []
        self.commonVariablesEE = []
        self.commonVariablesComb = []
        self.allVariables = []
        self.commonCuts = []
        self.target = "el_genOverSCEnergy"
        self.targetError = "1.4826*abs(BDTresponse - TARGET)"
        self.targetComb = "(el_generatedEnergy-el_scl_rawEnergy*(1+el_scl_preshowerEOverRaw)*BDTresponse)/(el_gsftrk_pAtVtx-el_scl_rawEnergy*(1+el_scl_preshowerEOverRaw)*BDTresponse)"
        self.commonCuts = []
        self.commonCutsEB = []
        self.commonCutsEE = []
        self.commonCutsComb = []
        self.commonCutsError = []
        self.regressions = {}
        self.parameterFile = ""
        self.info = ""
        self.version = ""
        self.simulate = False

    def initialize(self):
        self.findVersion()

    def printInfo(self):
        with open(self.outputDirectory+"/"+self.version+"/"+self.baseName+".info", "w") as info:
            print >>info, datetime.datetime.today().ctime()
            print >>info, "Parameter file = "+self.parameterFile
            print >>info, "Base name = "+self.baseName
            print >>info, "Trainer = "+self.trainerType
            print >>info, "Regression target = "+self.target
            print >>info, "Regression target error = "+self.targetError
            print >>info, "Regression target comb = "+self.targetComb
            index = 1
            for name,reg in sorted(self.regressions.iteritems()):
                print >>info, "Regression "+str(index)+": "+name
                options = ""
                for opt in reg.options:
                    options += opt+","
                print >>info, "  Options = "+options[:-1]
                cuts = ""
                for cut in reg.cuts:
                    cuts += cut+","
                print >>info, "  Cuts = "+cuts[:-1]
                cutsEB = ""
                for cut in reg.cutsEB:
                    cutsEB += cut+","
                print >>info, "  Cuts EB = "+cutsEB[:-1]
                cutsEE = ""
                for cut in reg.cutsEE:
                    cutsEE += cut+","
                print >>info, "  Cuts EE = "+cutsEE[:-1]
                cutsComb = ""
                for cut in reg.cutsComb:
                    cutsComb += cut+","
                print >>info, "  Cuts combi = "+cutsComb[:-1]
                index += 1
            print >>info, "\nAdditional information:"
            print >>info, self.info


    def addRegression(self, name):
        if name in self.regressions.keys():
            raise StandardError("ERROR: regression "+name+" has already been registered")
        reg = Regression()
        reg.id = 1
        reg.name = self.baseName+"_"+name
        reg.inputFiles = self.inputFiles
        reg.tree = self.tree
        reg.method = self.method
        if self.trainerType=="TMVA":
            reg.tmvaTrainingOptions = copy.copy(self.tmvaTrainingOptions)
        reg.options = copy.copy(self.commonOptions)
        reg.doErrors = self.doErrors
        reg.doCombine = self.doCombine
        reg.variablesEB = copy.copy(self.commonVariablesEB)
        reg.variablesEE = copy.copy(self.commonVariablesEE)
        reg.variablesComb = copy.copy(self.commonVariablesComb)
        reg.target = self.target
        reg.targetError = self.targetError
        reg.targetComb = self.targetComb
        reg.cuts = copy.copy(self.commonCuts)
        reg.cutsEB = copy.copy(self.commonCutsEB)
        reg.cutsEE = copy.copy(self.commonCutsEE)
        reg.cutsError = copy.copy(self.commonCutsError)
        reg.cutsComb = copy.copy(self.commonCutsComb)
        self.regressions[name] = reg

    def doCombination(self, name, combine):
        self.regressions[name].doCombine = combine

    def addTrainingOptions(self, name, options):
        self.regressions[name].tmvaTrainingOptions.extend(options)

    def addOptions(self, name, options):
        self.regressions[name].options.extend(options)

    def addVariablesEB(self, name, variables):
        self.regressions[name].variablesEB.extend(variables)
        for var in variables:
            if not var in self.allVariables:
                self.allVariables.append(var)

    def addVariablesEE(self, name, variables):
        self.regressions[name].variablesEE.extend(variables)
        for var in variables:
            if not var in self.allVariables:
                self.allVariables.append(var)

    def addVariablesComb(self, name, variables):
        self.regressions[name].variablesComb.extend(variables)
        for var in variables:
            if not var in self.allVariables:
                self.allVariables.append(var)


    def addCutsEB(self, name, cuts):
        if len(cuts)>0:
            self.regressions[name].cutsEB.extend(cuts)

    def addCutsEE(self, name, cuts):
        if len(cuts)>0:
            self.regressions[name].cutsEE.extend(cuts)

    def addCutsComb(self, name, cuts):
        if len(cuts)>0:
            self.regressions[name].cutsComb.extend(cuts)

    def addCuts(self, name, cuts):
        if len(cuts)>0:
            self.regressions[name].cuts.extend(cuts)

    def addCutsError(self, name, cuts):
        if len(cuts)>0:
            self.regressions[name].cutsError.extend(cuts)

    def setInputTree(self, name, tree):
        self.regressions[name].tree = tree

    def setTarget(self, name, target):
        self.regressions[name].target = target

    def setTargetError(self, name, target):
        self.regressions[name].targetError = target

    def setTargetComb(self, name, target):
        self.regressions[name].targetComb = target

    def prepareArea(self):
        ### perform some checks
        fullPathExe = self.baseDir+"/"+self.exe
        if not os.path.isfile(fullPathExe):
            raise StandardError("ERROR: cannot find regression executable {0}".format(fullPathExe))
        for lib in self.libs:
            if not os.path.isfile(lib):
                raise StandardError("ERROR: cannot find library {}".format(lib))
        if not (stat.S_IXUSR & os.stat(fullPathExe)[stat.ST_MODE]):
            raise StandardError("ERROR: you don't have exe permission for "+fullPathExe)
        for inFile in self.inputFiles:
            if not os.path.isfile(inFile):
                raise StandardError("ERROR: cannot find input file "+inFile)
        ### Create running area
        if not os.path.isdir(self.outputDirectory):
            os.makedirs(self.outputDirectory)
        if not os.path.isdir(self.outputDirectory+"/"+self.version):
            os.makedirs(self.outputDirectory+"/"+self.version)
        if not os.path.isdir(self.outputDirectory+"/"+self.version+"/jobs"):
            os.makedirs(self.outputDirectory+"/"+self.version+"/jobs")
        if not os.path.isdir(self.outputDirectory+"/"+self.version+"/obj"):
            os.makedirs(self.outputDirectory+"/"+self.version+"/obj")


        ### Create configurations files
        factoryOptions = ""
        for opt in self.tmvaFactoryOptions:
            factoryOptions += opt+":"
        factoryOptions = factoryOptions[0:-1]
        jobid = 1
        for name,reg in sorted(self.regressions.iteritems()):
            with open(self.outputDirectory+"/"+self.version+"/jobs/"+self.baseName+"."+str(jobid)+".config", 'w') as fConfig:
                print >>fConfig, "## Configuration file for regression", reg.name
                print >>fConfig, "Trainer: ",self.trainerType
                print >>fConfig, "NumberOfRegressions: 1"
                if self.trainerType=="TMVA":
                    print >>fConfig, "TMVAFactoryOptions:",factoryOptions
                print >>fConfig, "OutputDirectory:", "./"
                print >>fConfig, "\n"
                print >>fConfig, reg.toString()
                print >>fConfig, "\n\n"
            jobid += 1

        ### Copy libraries
        for lib in self.libs:
            shutil.copy(lib, self.outputDirectory+"/"+self.version+"/obj/")

        ### Copy parameter file
        shutil.copy(self.parameterFile, self.outputDirectory+"/"+self.version+"/"+self.baseName+".py")

    def printLaunchedJobs(self):
        if self.trainerType!="TMVA" and self.trainerType!="GBRTrain" and self.trainerType!="GBRLikelihoodTrain":
            raise StandardError("ERROR: Unknown trainer typer "+self.trainerType)
        print "\n>>>>>> You are going to send", len(self.regressions), self.trainerType, "regressions on batch <<<<<<"
        number = 1
        print "  "+self.baseName+"  Version "+self.version
        for name,reg in sorted(self.regressions.iteritems()):
            print "   ",str(number)+")",name
            number += 1
        print ">>>>>> You are going to send", len(self.regressions), "regressions on batch <<<<<<\n\n"

    def findVersion(self):
        self.version = "v_1_"+str(datetime.date.today())
        if os.path.isdir(self.outputDirectory):
            listdirs= [f for f in os.listdir(self.outputDirectory) if os.path.isdir(os.path.join(self.outputDirectory,f))]
            numberMax = 0
            for d in listdirs:
                number = int(d.split("_")[1])
                if number > numberMax:
                    numberMax = number
            self.version = "v_"+str(numberMax+1)+"_"+str(datetime.date.today())

    def findLastVersion(self):
        lastVersion = ""
        if os.path.isdir(self.outputDirectory):
            listdirs= [f for f in os.listdir(self.outputDirectory) if os.path.isdir(os.path.join(self.outputDirectory,f))]
            numberMax = 0
            for d in listdirs:
                number = int(d.split("_")[1])
                if number > numberMax:
                    numberMax = number
                    lastVersion = d
        return lastVersion

    def checkPreviousVersion(self):
        ok = True
        diffJobs = ""
        if os.path.isdir(self.outputDirectory):
            listdirs= [f for f in os.listdir(self.outputDirectory) if os.path.isdir(os.path.join(self.outputDirectory,f))]
            numberMax = 0
            lastVersion = ""
            for d in listdirs:
                number = int(d.split("_")[1])
                if number > numberMax:
                    numberMax = number
                    lastVersion = d
            pyFiles = [f for f in os.listdir(self.outputDirectory+"/"+lastVersion) if os.path.splitext(f)[1]==".py"]
            if len(pyFiles)==0:
                print "WARNING: Cannot find .py configuration file for version "+lastVersion
            elif len(pyFiles)>1:
                print "WARNING: Found more than 1 .py file for version "+lastVersion
            else:
                sys.path.append(self.outputDirectory+"/"+lastVersion)
                lastJobs = __import__(os.path.splitext(pyFiles[0])[0]).batch
                diffJobs = self.diff(lastJobs)
                if diffJobs!="":
                    ok = False
        return ok, diffJobs

    def launch(self, simulate=False):
        ok, diffJobs = self.checkPreviousVersion()
        if not ok:
            print "WARNING: Previous job parameters are different from the ones of the current jobs"
            print diffJobs
            answer = ""
            while answer!='y' and answer!='n':
                sys.stdout.write("Send jobs anyway (y/n)? ")
                answer = sys.stdin.readline()
                answer = answer.strip()
            if answer =='n':
                print "-- STOP\n"
                return
        if self.info=="":
            print "No information given for these jobs "+self.baseName
            answer = ""
            while answer!='y' and answer!='n':
                sys.stdout.write("Send jobs anyway (y/n)? ")
                answer = sys.stdin.readline()
                answer = answer.strip()
            if answer =='n':
                print "-- STOP\n"
                return
        self.printLaunchedJobs()
        answer = ""
        go = False
        print "Did you init your proxy?"
        while answer!='y' and answer!='n':
            sys.stdout.write("Do you want to continue (y/n)? ")
            answer = sys.stdin.readline()
            answer = answer.strip()
        if answer == 'y':
            go = True
            print "-- GO\n"
        elif answer =='n':
            print "-- STOP\n"
        if go:
            self.prepareArea()
            self.printInfo()
            #for name,reg in sorted(self.regressions.iteritems()):
            job = ParametricJobsMP("job_"+self.baseName)
            job.nJobs = len(self.regressions.items())
            job.exe = self.baseDir+"/"+self.exe
            job.parameterFileBaseName = self.baseName+".$jobid.config"
            jobid = 1
            for name,reg in sorted(self.regressions.iteritems()):
                job.inputFiles.append(self.baseName+"."+str(jobid)+".config")
                job.parameterFiles.append(self.baseName+"."+str(jobid)+".config")
                jobid += 1
            for lib in self.libs:
                job.libs.append(lib)
            #job.incl.append(self.includeDir+"/GBRForest.h")
            #job.incl.append(self.includeDir+"/GBRTree.h")
            #job.incl.append(self.includeDir+"/libDictionary.C")
            #job.src.append(self.srcDir+"/GBRForest.cpp")
            #job.src.append(self.srcDir+"/GBRTree.cpp")
            job.creatingDict = False
            job.outputDir = self.outputDirectory+"/"+self.version
            job.checkParams()
            job.prepareArea()
            job.createScript()
            job.prepareCommand()
            job.launch(self.simulate)

    def checkOutput(self):
        lastVersion = self.findLastVersion()
        ROOT.gSystem.Load("/home/llr/cms/sauvan/TMVAregression/obj/libDictionnary_C.so")
        summary = ""
        # check if outputs have been retrieved
        print ">> Checking parametric job", self.baseName, lastVersion
        print "  > Checking parametricJobs_status file and output directory"
        if not os.path.exists(self.outputDirectory+"/"+lastVersion+"/parametricJobs_status"):
            print "  ERROR: non existing status file available"
            return
        with open(self.outputDirectory+"/"+lastVersion+"/parametricJobs_status") as statusFile:
            lines = statusFile.readlines()
            for line in lines:
                if not "RETRIEVED" in line:
                    print "  ERROR: non-retrieved output"
                    print "  "+line
                    return
        if not os.path.isdir(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs"):
            print "  ERROR: non existing output directory"
            return
        # check stderr and stdout
        print "  > Checking stderr and stdout"
        index = 1
        for name,reg in sorted(self.regressions.iteritems()):
            ## stderr
            if not os.path.exists(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs/stderr."+str(index)):
                print "  ERROR: non-existing stderr."+str(index)+" (regression "+name+")"
                return
            if os.stat(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs/stderr."+str(index))[6]!=0:
                print "  ERROR: non-empty stderr."+str(index)+" (regression "+name+")"
                print "  Details:"
                with open(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs/stderr."+str(index)) as stderr:
                    lines = stderr.readlines()
                    for line in lines:
                        print "   "+line
            ## stdout
            isFinish = False
            nTrees = 0
            nVars = 0
            nEvts = 0
            time = 0
            if not os.path.exists(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs/stdout."+str(index)):
                print "  ERROR: non-existing stdout."+str(index)+" (regression "+name+")"
                return
            if os.stat(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs/stdout."+str(index))[6]==0:
                print "  ERROR: empty stdout."+str(index)+" (regression "+name+")"
                return
            with open(self.outputDirectory+"/"+lastVersion+"/parametricJobs_outputs/stdout."+str(index)) as stdout:
                lines = stdout.readlines()
                for line in lines:
                    if "Finish - All good" in line:
                        isFinish = True
                    if "tree " in line:
                        n = int(line.split()[1])
                        if n>nTrees:
                            nTrees = n
                    if "nev " in line and "nvar " in line:
                        nev = int(line.split(", ")[0].split(" = ")[1])
                        nvar = int(line.split(", ")[1].split(" = ")[1])
                        nVars = nvar
                        nEvts = nev
                    if "RegressionManager::makeRegression(): Elapsed time" in line:
                        time = int(line.split(" = ")[1].split()[0])
            if not isFinish:
                print "  ERROR: not ended job "+str(index)+" (regression "+name+")"
                return
            #if nVars!=len(reg.variables):
            #    print "ERROR: Number of variables in log file ("+str(nVars)+") not compatible with the regression definition ("+str(len(reg.variables))+")"
            #    return
            print "  Regression "+name+":"
            print "   - Processed events = "+str(nEvts)
            print "   - Input variables  = "+str(nVars)
            print "   - Number of trees  = "+str(nTrees)
            print "   - Elapsed time     = "+str(datetime.timedelta(seconds=time))
            summary += "  Regression "+name+":\n"
            summary += "   - Processed events = "+str(nEvts)+"\n"
            summary += "   - Input variables  = "+str(nVars)+"\n"
            summary += "   - Number of trees  = "+str(nTrees)+"\n"
            summary += "   - Elapsed time     = "+str(datetime.timedelta(seconds=time))+"\n"
            index += 1
        ## check result files
        print "  > Checking result files"
        index = 1
        for name,reg in sorted(self.regressions.iteritems()):
            fileName = self.outputDirectory+"/"+lastVersion+"/"+self.baseName+"_"+name+"_results.root"
            results = ROOT.TFile.Open(fileName)
            if not results:
                print "ERROR: cannot open result file for regression "+name
                return
            if not results.FindKey("gbrtrainEB"):
                print "ERROR: cannot find forest in result file for regression "+name
                return
            if not results.FindKey("varlistEB"):
                print "ERROR: cannot find variable list in result file for regression "+name
                return
            varlist = ROOT.MakeNullPointer( ROOT.vector("string") )
            results.GetObject("varlistEB", varlist)
            nVars = len(varlist)
            #if nVars!=len(reg.variables):
            #    print "ERROR: Number of variables in result file ("+str(nVars)+") not compatible with the regression definition ("+str(len(reg.variables))+")"
            #    return
            results.Close()
            print "  "+self.baseName+"_"+name+"_results.root --> %7.2f Mb"%(os.path.getsize(fileName)/1.e6)
            summary += "  "+self.baseName+"_"+name+"_results.root --> %7.2f Mb\n"%(os.path.getsize(fileName)/1.e6)
        print "\n Everything's OK"
        with open(self.outputDirectory+"/"+lastVersion+"/"+self.baseName+".summary", 'w') as summaryFile:
            print >>summaryFile, datetime.datetime.today().ctime()
            print >>summaryFile, summary

    def diff(self,other):
        diff = ""
        if self.baseName!=other.baseName:
            diff += "Different base names: "+self.baseName+" != "+other.baseName+"\n"
        if len(self.regressions)!=len(other.regressions):
            diff += "Different numbers of regressions\n"
        diffReg1 = list(set(self.regressions.keys()) - set(other.regressions.keys()))
        diffReg2 = list(set(other.regressions.keys()) - set(self.regressions.keys()))
        if len(diffReg1)>0 or len(diffReg2)>0:
            diff += "Different regression names:\n"
            if len(diffReg1)>0:
                diff += "   New: "
                for reg in diffReg1:
                    diff += reg+" "
                diff += "\n"
            if len(diffReg2)>0:
                diff += "   Old: "
                for reg in diffReg2:
                    diff += reg+" "
                diff += "\n"
        else:
            for name,reg1 in sorted(self.regressions.iteritems()):
                reg2 = other.regressions[name]
                diffReg = reg1.diff(reg2)
                if diffReg!="":
                    diff += "------------------------------\n"
                    diff += "Differences in "+name+"\n"
                    diff += diffReg
        return diff



