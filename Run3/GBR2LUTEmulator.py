import ROOT
import array
import copy
import os
import math
from itertools import product
from ROOT import gStyle

class GBR2LUTEmulator:
    def __init__(self):
        self.name = ""
        self.inputFileName = "regression.root"
        self.inputFile = None
        self.forestName = "EBCorrection"
        self.forest = None
        self.inputVariables = []
        self.variablePoints = []
        self.sortedShapes = "../ShapeSorting/sortedShapes.txt"
        self.outputFileName = "corrections.txt"

    def retrieveForest(self):
        cmsswBase = os.environ["CMSSW_BASE"]
        scram_arch = os.environ["SCRAM_ARCH"]
        ROOT.gSystem.AddDynamicPath("{0}/lib/{1}/".format(cmsswBase, scram_arch))
        ROOT.gSystem.Load("libHiggsAnalysisGBRLikelihood.so") 
        #print "inputFileName",self.inputFileName
        self.inputFile = ROOT.TFile.Open(self.inputFileName)
        self.forest = ROOT.MakeNullPointer( "GBRForest" ) 
        #self.forest = ROOT.MakeNullPointer( "GBRForestD" ) 
        self.inputFile.GetObject(self.forestName, self.forest)
        #
        varlist = ROOT.MakeNullPointer( ROOT.vector("string") )
        self.inputFile.GetObject("varlistEB", varlist)
        #if len(varlist)!=2:
        #    raise StandardError("ERROR: Number of input variables != 2. Not implemented for the moment.")
        for i,name in enumerate(varlist):
            print "var ",i
            print name
            print self.variablePoints[i][0]
            if self.variablePoints[i][0]!=name:
                raise StandardError("ERROR: Input variables are not given in the correct order in self.variablePoints")
            self.inputVariables.append(name)

    def createLUT(self, header=""):
        sortedShapes = {}
        with open(self.sortedShapes, 'r') as f:
            lines = f.readlines()
            for line in lines:
                a,b = line.split()
                if not int(a) in sortedShapes:
                    sortedShapes[int(a)] = int(b)
                else:
                    print "WARNING: shape already filled"
        with open(self.outputFileName, 'w') as output:
            inputPoints = [ it[1] for it in self.variablePoints ]
            address = 0
            if header=="":
                header = """
# Calibration vs |ieta|,ET,hasEM,isMerged Derived from 13 TeV data (RunD) using mu+tau Tag & Probe selection, with semi-parametric regression
# The LUT output is (ET_off/ET_L1 - 0.5) between 0 and 2, encoded on 9 bits: corr = LUT_output/256. + 0.5 <-> pt_calib = (pt_raw>>1)+ ((LUT_output*pt_raw)>>8)
# Index is (BinEta<<6)+(BinEt<<2)+(BinHasEM<<1)+BinIsMerged
# anything after # is ignored with the exception of the header
# the header is first valid line starting with #<header> versionStr(unused but may be in future) </header>
#<header> V3 </header>
"""
# Calibration vs |ieta|,shape,E. Derived from Zee PU40bx25, with semi-parametric regression
# The LUT output is (ET_off/ET_L1 - 1) between O and 1, encoded on 9 bits
# Index is shape+E*128+(|ieta|-1)*256*128. 0 is |ieta|=1,shape=0,E=0 (0 GeV)
#anything after # is ignored with the exception of the header
#the header is first valid line starting with #<header> versionStr(unused but may be in future) nrBitsAddress nrBitsData </header>
#<header> V2 20 9 </header>
#"""
            print >>output, header
            for inputs in product(*inputPoints):
                # convert into sorted shapes for the BRT
                inputsCopy = list(copy.copy(inputs))
                for i,var in enumerate(self.inputVariables):
                    if "sorted" in var:
                        inputsCopy[i] = sortedShapes[inputs[i]]
                resp = 1.
                if self.forest:
                    #resp = self.forest.GetResponse(array.array('f',inputs))
                    resp = self.forest.GetResponse(array.array('f',inputsCopy))
                    #print "resp = ",resp
                    ## transformation for the semi-parametric regression
                    low  = 0.2;
                    high = 5.;
                    #high = 2.;
                    scale = 0.5*(high-low);
                    offset = low + 0.5*(high-low);
                    corr = offset + scale*math.sin(resp);
                    #print "corr = ",corr
                strInputs = ""
                for var,inp in zip(self.inputVariables,inputs):
                    strInputs += var+"="+str(inp)+" "
                #print "resp = ",resp
                #print "corr = ",corr
                #y = corr-1. # (ET_off/ET_L1) - 1 #was here JB
                #y = resp-1.#was here before re-encoding 04/01/16
                if resp>2.1: resp=2.1
                y = resp-0.5
                #print "y = ",y
                #yint = int(round(y*512)) # [0,1] -> [0, 512] (9 bits)
                yint = int(round(y/2.*512)) # [0,2] -> [0, 512] (9 bits)
                #print "yint = ",yint
                if yint>511: # threshold
                    yint=511
                print >>output, address, yint, "#", strInputs
                address += 1

    def createTH1(self):
        ## Creating binning for TH2F filling
        bins1 = []
        npoints1 = len(self.variablePoints[0][1])
        #
        binlow = self.variablePoints[0][1][0] - (self.variablePoints[0][1][1]-self.variablePoints[0][1][0])/2.
        bins1.append(binlow)
        for v1,v2 in zip(self.variablePoints[0][1], self.variablePoints[0][1][1:]):
            binlow = (v1+v2)/2.
            bins1.append(binlow)
        binhigh = self.variablePoints[0][1][npoints1-1] + (self.variablePoints[0][1][npoints1-1]-self.variablePoints[0][1][npoints1-2])/2.
        bins1.append(binhigh)
        #

        self.lut = ROOT.TH1F("LUT_"+self.name, "LUT", len(bins1)-1, array.array('f',bins1))

        with open(self.outputFileName, 'w') as output:
            inputPoints = [ it[1] for it in self.variablePoints ]
            for inputs in product(*inputPoints):
                value1 = inputs[0]
                if self.variablePoints[0][1]: value1 = self.variablePoints[0][1][value1] ## value mapping

                resp = self.forest.GetResponse(array.array('f',inputs))
                bin1 = self.lut.GetXaxis().FindBin(value1)

                #added by Olivier
                low  = 0.2;
                high = 5.;
                #high = 2.;
                scale = 0.5*(high-low);
                offset = low + 0.5*(high-low);
                corr = offset + scale*math.sin(resp);
                corr = resp
                self.lut.SetBinContent(bin1, corr)
                print "bin1 = ",bin1
                print "corr = ",corr
                strInputs = ""
                for inp in inputs:
                    strInputs += str(inp)+" "
                print >>output, strInputs, resp

        outputRooFile = ROOT.TFile(self.outputFileName,"RECREATE")
        outCanvasName = self.outputFileName
        #print  self.outputFileName
        outCanvasName = outCanvasName.replace(".root",".pdf")
        print outCanvasName
        canvas = ROOT.TCanvas("c_test","c_test",800,800)

        outputRooFile.cd()
        YaxisString = "L1_{E_{T}}/#tau_{E_{T}} Correction";
        XaxisString = "IEt [0.5 GeV]";         
        self.lut.GetXaxis().SetTitle(XaxisString)
        self.lut.GetXaxis().SetTitleOffset(1.2)
        self.lut.GetYaxis().SetTitle(YaxisString)
        self.lut.GetYaxis().SetTitleOffset(1.43)
        self.lut.GetYaxis().SetRangeUser(0.,200.);
        self.lut.SetMinimum(0.5)
        self.lut.SetMaximum(3.3)
        #if self.outputFileName.find("Cluster_iEta_Cluster_Shape"):
        #    self.lut.SetMinimum(1.3)
        #    self.lut.SetMaximum(2.2)

        self.lut.Write()

        gStyle.SetOptStat(0)
        canvas.SetLeftMargin(0.15);
        canvas.SetRightMargin(0.15);
        self.lut.Draw("")
        print outCanvasName
        canvas.SaveAs(outCanvasName)


    def createTH2(self):
        ## Creating binning for TH2F filling
        bins1 = []
        bins2 = []
        npoints1 = len(self.variablePoints[0][1])
        npoints2 = len(self.variablePoints[1][1])
        #
        binlow = self.variablePoints[0][1][0] - (self.variablePoints[0][1][1]-self.variablePoints[0][1][0])/2.
        bins1.append(binlow)
        for v1,v2 in zip(self.variablePoints[0][1], self.variablePoints[0][1][1:]):
            binlow = (v1+v2)/2.
            bins1.append(binlow)
        binhigh = self.variablePoints[0][1][npoints1-1] + (self.variablePoints[0][1][npoints1-1]-self.variablePoints[0][1][npoints1-2])/2.
        bins1.append(binhigh)
        #
        binlow = self.variablePoints[1][1][0] - (self.variablePoints[1][1][1]-self.variablePoints[1][1][0])/2.
        bins2.append(binlow)
        for v1,v2 in zip(self.variablePoints[1][1], self.variablePoints[1][1][1:]):
            binlow = (v1+v2)/2.
            bins2.append(binlow)
        binhigh = self.variablePoints[1][1][npoints2-1] + (self.variablePoints[1][1][npoints2-1]-self.variablePoints[1][1][npoints2-2])/2.
        bins2.append(binhigh)

        self.lut = ROOT.TH2F("LUT_"+self.name, "LUT", len(bins1)-1, array.array('f',bins1), len(bins2)-1, array.array('f',bins2))

        with open(self.outputFileName, 'w') as output:
            inputPoints = [ it[1] for it in self.variablePoints ]
            for inputs in product(*inputPoints):
                value1 = inputs[0]
                if self.variablePoints[0][1]: value1 = self.variablePoints[0][1][value1] ## value mapping
                #if self.variablePoints[0][2]: value1 = self.variablePoints[0][2][value1] ## value mapping
                value2 = inputs[1]
                if self.variablePoints[1][1]: value2 = self.variablePoints[1][1][value2] ## value mapping
                #if self.variablePoints[1][2]: value2 = self.variablePoints[1][2][value2] ## value mapping
                resp = self.forest.GetResponse(array.array('f',inputs))
                bin1 = self.lut.GetXaxis().FindBin(value1)
                bin2 = self.lut.GetYaxis().FindBin(value2)
                #added by Olivier
                low  = 0.2;
                high = 5.;
                #high = 2.;
                scale = 0.5*(high-low);
                offset = low + 0.5*(high-low);
                corr = offset + scale*math.sin(resp);
                corr = resp
                self.lut.SetBinContent(bin1, bin2, corr)
                #self.lut.SetBinContent(bin1, bin2, resp)
                strInputs = ""
                for inp in inputs:
                    strInputs += str(inp)+" "
                print >>output, strInputs, resp

        outputRooFile = ROOT.TFile(self.outputFileName,"RECREATE")
        outCanvasName = self.outputFileName
        print  self.outputFileName
        outCanvasName = outCanvasName.replace(".root",".pdf")
        #print outCanvasName
        canvas = ROOT.TCanvas("c_test","c_test",800,800)

        outputRooFile.cd()
        if self.outputFileName.find("compressedieta_compressedE"):
            XaxisString = "i#eta(L1 #tau)";
            YaxisString = "hwEt(L1 #tau)";
        if self.outputFileName.find("IEt_compressedshape"):
            XaxisString = "Shape index (compressed)";
            YaxisString = "hwEt(L1 #tau)";
        if self.outputFileName.find("IEt_Cluster_Shape"):
            XaxisString = "Shape index";
            YaxisString = "hwEt(L1 #tau)";
        if self.outputFileName.find("Cluster_iEta_Cluster_Shape"):
            YaxisString = "Shape index";
            XaxisString = "i#eta(L1 #tau)";
        if self.outputFileName.find("compressedieta_compressedshape"):
            YaxisString = "Compressed shape index";
            XaxisString = "Compressed i#eta(L1 #tau)";
        if self.outputFileName.find("compressedieta_compressedE"):
            YaxisString = "Compressed hwEt(L1 #tau)";
            XaxisString = "Compressed i#eta(L1 #tau)";         
        self.lut.GetXaxis().SetTitle(XaxisString)
        self.lut.GetXaxis().SetTitleOffset(1.2)
        self.lut.GetYaxis().SetTitle(YaxisString)
        self.lut.GetYaxis().SetTitleOffset(1.43)
        self.lut.GetYaxis().SetRangeUser(0.,200.);
        self.lut.SetMinimum(0.5)
        self.lut.SetMaximum(3.3)
        if self.outputFileName.find("Cluster_iEta_Cluster_Shape"):
            self.lut.SetMinimum(1.3)
            self.lut.SetMaximum(2.2)
        if self.outputFileName.find("Cluster_compressedieta_compressedshape"):
            self.lut.SetMinimum(0.8)
            self.lut.SetMaximum(2.5)
            self.lut.GetYaxis().SetRangeUser(0.,15.);
            self.lut.GetXaxis().SetRangeUser(0.,15.);
        if self.outputFileName.find("compressedieta_compressedE"):
            self.lut.SetMinimum(0.8)
            self.lut.SetMaximum(2.5)
            self.lut.GetYaxis().SetRangeUser(0.,15.);
            self.lut.GetXaxis().SetRangeUser(0.,15.);

        if self.outputFileName.find("pt18")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>18 GeV")
        if self.outputFileName.find("pt19")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>19 GeV")
        if self.outputFileName.find("pt20")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>20 GeV")
        if self.outputFileName.find("pt21")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>21 GeV")
        if self.outputFileName.find("pt22")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>22 GeV")
        if self.outputFileName.find("pt23")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>23 GeV")
        if self.outputFileName.find("pt24")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>24 GeV")
        if self.outputFileName.find("pt25")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>25 GeV")
        if self.outputFileName.find("pt26")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>26 GeV")
        if self.outputFileName.find("pt27")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>27 GeV")
        if self.outputFileName.find("pt28")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>28 GeV")
        if self.outputFileName.find("pt29")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>29 GeV")
        if self.outputFileName.find("pt30")!=-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>30 GeV")
        if self.outputFileName.find("pt")==-1:
            self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>18 GeV")
        self.lut.Write()

        gStyle.SetOptStat(0)
        canvas.SetLeftMargin(0.15);
        canvas.SetRightMargin(0.15);
        self.lut.Draw("colz")
        canvas.SaveAs(outCanvasName)


    def createTH3(self):
        ## Creating binning for TH2F filling
        bins1 = []
        bins2 = []
        bins3 = []
        npoints1 = len(self.variablePoints[0][1])
        npoints2 = len(self.variablePoints[1][1])
        npoints3 = len(self.variablePoints[2][1])

        print "npoints1",npoints1
        print "npoints2",npoints2
        print "npoints3",npoints3
        #
        binlow = self.variablePoints[0][1][0] - (self.variablePoints[0][1][1]-self.variablePoints[0][1][0])/2.
        bins1.append(binlow)
        for v1,v2 in zip(self.variablePoints[0][1], self.variablePoints[0][1][1:]):
            binlow = (v1+v2)/2.
            bins1.append(binlow)
        binhigh = self.variablePoints[0][1][npoints1-1] + (self.variablePoints[0][1][npoints1-1]-self.variablePoints[0][1][npoints1-2])/2.
        bins1.append(binhigh)
        #
        binlow = self.variablePoints[1][1][0] - (self.variablePoints[1][1][1]-self.variablePoints[1][1][0])/2.
        bins2.append(binlow)
        for v1,v2 in zip(self.variablePoints[1][1], self.variablePoints[1][1][1:]):
            binlow = (v1+v2)/2.
            bins2.append(binlow)
        binhigh = self.variablePoints[1][1][npoints2-1] + (self.variablePoints[1][1][npoints2-1]-self.variablePoints[1][1][npoints2-2])/2.
        bins2.append(binhigh)

        binlow = self.variablePoints[1][1][0] - (self.variablePoints[1][1][1]-self.variablePoints[1][1][0])/2.
        bins3.append(binlow)
        for v1,v2 in zip(self.variablePoints[1][1], self.variablePoints[1][1][1:]):
            binlow = (v1+v2)/2.
            bins3.append(binlow)
        binhigh = self.variablePoints[2][1][npoints3-1] + (self.variablePoints[2][1][npoints3-1]-self.variablePoints[2][1][npoints3-2])/2.
        bins3.append(binhigh)

        self.lut = ROOT.TH3F("LUT_"+self.name, "LUT", len(bins1)-1, array.array('f',bins1), len(bins2)-1, array.array('f',bins2), len(bins3)-1, array.array('f',bins3))

        with open(self.outputFileName, 'w') as output:
            inputPoints = [ it[1] for it in self.variablePoints ]
            for inputs in product(*inputPoints):
                value1 = inputs[0]
                if self.variablePoints[0][1]: value1 = self.variablePoints[0][1][value1] ## value mapping
                #if self.variablePoints[0][2]: value1 = self.variablePoints[0][2][value1] ## value mapping
                value2 = inputs[1]
                if self.variablePoints[1][1]: value2 = self.variablePoints[1][1][value2] ## value mapping
                #if self.variablePoints[1][2]: value2 = self.variablePoints[1][2][value2] ## value mapping
                value3 = inputs[2]
                if self.variablePoints[2][1]: value3 = self.variablePoints[2][1][value3] ## value mapping

                resp = self.forest.GetResponse(array.array('f',inputs))
                bin1 = self.lut.GetXaxis().FindBin(value1)
                bin2 = self.lut.GetYaxis().FindBin(value2)
                bin3 = self.lut.GetZaxis().FindBin(value3)
                print bin1,"",bin2,"",bin3," --> ",resp
                #added by Olivier
                low  = 0.2;
                high = 5.;
                #high = 2.;
                scale = 0.5*(high-low);
                offset = low + 0.5*(high-low);
                corr = offset + scale*math.sin(resp);
                corr = resp
                #print bin1,"",bin2,"",bin3," --> ",corr
                self.lut.SetBinContent(bin1, bin2, bin3, corr)
                #self.lut.SetBinContent(bin1, bin2, resp)
                strInputs = ""
                for inp in inputs:
                    strInputs += str(inp)+" "
                print >>output, strInputs, resp

        outputRooFile = ROOT.TFile(self.outputFileName,"RECREATE")
        outCanvasName = self.outputFileName
        print  self.outputFileName
        outCanvasName = outCanvasName.replace(".root",".pdf")
        #print outCanvasName
        canvas = ROOT.TCanvas("c_test","c_test",800,800)

        outputRooFile.cd()

        if self.outputFileName.find("compressedieta_compressedE"):
            XaxisString = "Compressed i#eta(L1 #tau)";
            YaxisString = "Comprssed hwEt(L1 #tau)";
            ZaxisString = "Compressed shape index";

        self.lut.GetXaxis().SetTitle(XaxisString) 
        self.lut.GetYaxis().SetTitle(YaxisString) 
        self.lut.GetZaxis().SetTitle(ZaxisString) 

        #if self.outputFileName.find("pt18")!=-1:
        #    self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>18")

        self.lut.Write()

        gStyle.SetOptStat(0)
        canvas.SetLeftMargin(0.15);
        canvas.SetRightMargin(0.15);
        self.lut.Draw("lego")
        canvas.SaveAs(outCanvasName)






    def createTH4(self):
        ## Creating binning for TH2F filling
        bins1 = []
        bins2 = []
        bins3 = []
        bins4 = []
        npoints1 = len(self.variablePoints[0][1])
        npoints2 = len(self.variablePoints[1][1])
        npoints3 = len(self.variablePoints[2][1])
        npoints4 = len(self.variablePoints[3][1])

        print "npoints1",npoints1
        print "npoints2",npoints2
        print "npoints3",npoints3
        print "npoints4",npoints4
        #
        binlow = self.variablePoints[0][1][0] - (self.variablePoints[0][1][1]-self.variablePoints[0][1][0])/2.
        bins1.append(binlow)
        for v1,v2 in zip(self.variablePoints[0][1], self.variablePoints[0][1][1:]):
            binlow = (v1+v2)/2.
            bins1.append(binlow)
        binhigh = self.variablePoints[0][1][npoints1-1] + (self.variablePoints[0][1][npoints1-1]-self.variablePoints[0][1][npoints1-2])/2.
        bins1.append(binhigh)
        #
        binlow = self.variablePoints[1][1][0] - (self.variablePoints[1][1][1]-self.variablePoints[1][1][0])/2.
        bins2.append(binlow)
        for v1,v2 in zip(self.variablePoints[1][1], self.variablePoints[1][1][1:]):
            binlow = (v1+v2)/2.
            bins2.append(binlow)
        binhigh = self.variablePoints[1][1][npoints2-1] + (self.variablePoints[1][1][npoints2-1]-self.variablePoints[1][1][npoints2-2])/2.
        bins2.append(binhigh)

        binlow = self.variablePoints[2][1][0] - (self.variablePoints[1][1][1]-self.variablePoints[1][1][0])/2.
        bins3.append(binlow)
        for v1,v2 in zip(self.variablePoints[2][1], self.variablePoints[2][1][1:]):
            binlow = (v1+v2)/2.
            bins3.append(binlow)
        binhigh = self.variablePoints[2][1][npoints3-1] + (self.variablePoints[2][1][npoints3-1]-self.variablePoints[2][1][npoints3-2])/2.
        bins3.append(binhigh)

        binlow = self.variablePoints[3][1][0] - (self.variablePoints[1][1][1]-self.variablePoints[1][1][0])/2.
        bins4.append(binlow)
        for v1,v2 in zip(self.variablePoints[3][1], self.variablePoints[3][1][1:]):
            binlow = (v1+v2)/2.
            bins4.append(binlow)
        binhigh = self.variablePoints[3][1][npoints3-1] + (self.variablePoints[3][1][npoints3-1]-self.variablePoints[3][1][npoints3-2])/2.
        bins4.append(binhigh)

        self.lut = ROOT.TH3F("LUT_isMerged0_"+self.name, "LUT_isMerged0", len(bins1)-1, array.array('f',bins1), len(bins2)-1, array.array('f',bins2), len(bins3)-1, array.array('f',bins3))
        self.lut2 = ROOT.TH3F("LUT_isMerged1_"+self.name, "LUT_isMerged1", len(bins1)-1, array.array('f',bins1), len(bins2)-1, array.array('f',bins2), len(bins3)-1, array.array('f',bins3))

        with open(self.outputFileName, 'w') as output:
            inputPoints = [ it[1] for it in self.variablePoints ]
            for inputs in product(*inputPoints):
                value1 = inputs[0]
                if self.variablePoints[0][1]: value1 = self.variablePoints[0][1][value1] ## value mapping
                #if self.variablePoints[0][2]: value1 = self.variablePoints[0][2][value1] ## value mapping
                value2 = inputs[1]
                if self.variablePoints[1][1]: value2 = self.variablePoints[1][1][value2] ## value mapping
                #if self.variablePoints[1][2]: value2 = self.variablePoints[1][2][value2] ## value mapping
                value3 = inputs[2]
                if self.variablePoints[2][1]: value3 = self.variablePoints[2][1][value3] ## value mapping

                value4 = inputs[3]
                if self.variablePoints[3][1]: value4 = self.variablePoints[3][1][value4] ## value mapping

                resp = self.forest.GetResponse(array.array('f',inputs))
                bin1 = self.lut.GetXaxis().FindBin(value1)
                bin2 = self.lut.GetYaxis().FindBin(value2)
                bin3 = self.lut.GetZaxis().FindBin(value3)
                bin4 = 0
                if(value4<0.5):
                    bin4 = 1
                else:
                    bin4 = 2

                print bin1,"",bin2,"",bin3,"",bin4," --> ",resp
                #added by Olivier
                low  = 0.2;
                high = 5.;
                #high = 2.;
                scale = 0.5*(high-low);
                offset = low + 0.5*(high-low);
                corr = offset + scale*math.sin(resp);
                corr = resp
                if corr>=2.1: corr=2.1
                #print bin1,"",bin2,"",bin3," --> ",corr
                if bin4==1:
                    self.lut.SetBinContent(bin1, bin2, bin3, corr)
                    print "filling histo1"
                else:
                    self.lut2.SetBinContent(bin1, bin2, bin3, corr)
                    print "filling histo2"
                #self.lut.SetBinContent(bin1, bin2, resp)
                strInputs = ""
                for inp in inputs:
                    strInputs += str(inp)+" "
                print >>output, strInputs, resp

        outputRooFile = ROOT.TFile(self.outputFileName,"RECREATE")
        outCanvasName = self.outputFileName
        print  self.outputFileName
        outCanvasName = outCanvasName.replace(".root",".pdf")
        #print outCanvasName

        outputRooFile.cd()

        if self.outputFileName.find("compressedieta_compressedE"):
            XaxisString = "Compressed i#eta(L1 #tau)";
            YaxisString = "Comprssed hwEt(L1 #tau)";
            ZaxisString = "Cluster has EM energy";

        self.lut.GetXaxis().SetTitle(XaxisString) 
        self.lut.GetYaxis().SetTitle(YaxisString) 
        self.lut.GetZaxis().SetTitle(ZaxisString) 

        self.lut2.GetXaxis().SetTitle(XaxisString) 
        self.lut2.GetYaxis().SetTitle(YaxisString) 
        self.lut2.GetZaxis().SetTitle(ZaxisString) 

        #if self.outputFileName.find("pt18")!=-1:
        #    self.lut.SetTitle("Calibration constants: p_{T}(Off. #tau)>18")

        self.lut.Write()
        self.lut2.Write()

        #canvas = ROOT.TCanvas("c_test","c_test",800,800)
        #gStyle.SetOptStat(0)
        #canvas.SetLeftMargin(0.15);
        #canvas.SetRightMargin(0.15);
        #self.lut.Draw("lego")
        #canvas.SaveAs(outCanvasName)

        
