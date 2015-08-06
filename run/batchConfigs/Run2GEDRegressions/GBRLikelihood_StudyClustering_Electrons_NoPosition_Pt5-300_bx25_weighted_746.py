from RegressionParametricJobLauncher import RegressionParametricJobLauncher
import glob
import os

info = """
GBR likelihood regressison
Based on 7_4_6
PU 1-50, bx 25
Version without position in input variables

SC regressions with weights
* PF-Mustache + Refined SC electron
* Barrel + Endcap

All cluster energies are non-corrected

RandomSeed=nEvent+123456*el_index
EventWeight=
min(1,exp(-(genPt-50)/50))

Correction + error + combine

"""

cmsswBase    = os.environ["CMSSW_BASE"]
cmsswRelBase = os.environ["CMSSW_RELEASE_BASE"]
scram_arch   = os.environ["SCRAM_ARCH"]

inputFiles = []
inputFiles.extend(glob.glob("/data_CMS/cms/sauvan/Regression/Ntuples/DoubleElectron_FlatPt-1To300/regression_ntuple_RunIISpring15DR74-AsymptFlat0to50bx25RawReco_MCRUN2_74_V9A-v3_2015_08_03/regression_ntuple_746_*.root"))


commonVariables = [
    "nVtx",
    "scRawEnergy",
    #"scEta",
    #"scPhi",
    "scEtaWidth",
    "scPhiWidth",
    "scSeedR9",
    "scSeedRawEnergy/scRawEnergy",
    "scSeedEmax/scRawEnergy",
    "scSeedE2nd/scRawEnergy",
    "scSeedLeftRightAsym",
    "scSeedTopBottomAsym",
    "scSeedSigmaIetaIeta",
    "scSeedSigmaIetaIphi",
    "scSeedSigmaIphiIphi",
    "N_ECALClusters",
    "clusterMaxDR",
    "clusterMaxDRDPhi",
    "clusterMaxDRDEta",
    "clusterMaxDRRawEnergy/scRawEnergy",
    "clusterRawEnergy[0]/scRawEnergy",
    "clusterRawEnergy[1]/scRawEnergy",
    "clusterRawEnergy[2]/scRawEnergy",
    "clusterDPhiToSeed[0]",
    "clusterDPhiToSeed[1]",
    "clusterDPhiToSeed[2]",
    "clusterDEtaToSeed[0]",
    "clusterDEtaToSeed[1]",
    "clusterDEtaToSeed[2]",
]


barrelVariables = [
    #"scSeedCryEta",
    #"scSeedCryPhi",
    #"scSeedCryIeta",
    #"scSeedCryIphi",
    "scSeedCryIetaV2",
    "scSeedCryIphiV2",
]

endcapVariables = [
    "scPreshowerEnergy/scRawEnergy",
    "scSeedCryIxV2",
    "scSeedCryIyV2",
]

combineVariables = [
    "(scRawEnergy+scPreshowerEnergy)*BDTresponse",
    "BDTerror/BDTresponse",
    "trkMomentum",
    "trkMomentumRelError",
    "BDTerror/BDTresponse/trkMomentumRelError",
    "(scRawEnergy+scPreshowerEnergy)*BDTresponse/trkMomentum",
    "(scRawEnergy+scPreshowerEnergy)*BDTresponse/trkMomentum*sqrt(BDTerror/BDTresponse*BDTerror/BDTresponse+trkMomentumRelError*trkMomentumRelError)", ## error on E/p ( E/p *sqrt((dE/E)^2 + (dp/p)^2) )
    "eleEcalDriven",
    "eleTrackerDriven",
    "eleClass",
    "scIsEB",
]



commonVariablesEB = commonVariables + barrelVariables
commonVariablesEE = commonVariables + endcapVariables
commonVariablesComb = combineVariables



batch = RegressionParametricJobLauncher()
batch.trainerType = "GBRLikelihoodTrain"
batch.baseDir = "/home/llr/cms/sauvan/Regression/RegressionTraining/"
batch.baseName = "GBRLikelihood_Clustering_746_bx25_Electrons_NoPosition"
batch.libs.append("{0}/lib/{1}/libHiggsAnalysisGBRLikelihood.so".format(cmsswBase, scram_arch))
batch.libs.append("{0}/lib/{1}/HiggsAnalysisGBRLikelihood_yr_rdict.pcm".format(cmsswBase, scram_arch))
batch.inputFiles = inputFiles
batch.outputDirectory = "/data_CMS/cms/sauvan/RegressionResults/StudyClustering/746/GBRLikelihood_Clustering_bx25_Electrons_NoPosition_PtCut50_PtSlope50_Sig5/"
batch.histoConfig = "../../dummy_Histo.config"
batch.commonOptions = ["MinEvents=200","Shrinkage=0.1","NTrees=1000","MinSignificance=5.0","RandomSeed=eventNumber+123456*scIndex","EventWeight=min(1,exp(-(genPt-50)/50))"]
batch.commonVariablesEB = commonVariablesEB
batch.commonVariablesEE = commonVariablesEE
batch.commonVariablesComb = commonVariablesComb
batch.commonCutsEB      = ["scIsEB"]
batch.commonCutsEE      = ["!scIsEB"]
batch.commonCuts        = ["(eventNumber%2==0)&&(isMatched==1)"]
batch.commonCutsError   = ["(eventNumber%2!=0)&&(((eventNumber-1)/2)%4==3)&&(isMatched==1)"]
#batch.commonCutsComb     = ["((eventNumber%2!=0)&&(((eventNumber-1)/2)%4!=3)&&(isMatched==1))"]
batch.commonCutsComb     = ["((eventNumber%2!=0)&&(((eventNumber-1)/2)%4>=2)&&(isMatched==1))"] ## not for production -- only testing
batch.doErrors = False
batch.doCombine = False

batch.target = "genEnergy/(scRawEnergy+scPreshowerEnergy)"
batch.targetError = "1.253*abs(BDTresponse - TARGET)"
batch.targetComb = "(genEnergy-(scRawEnergy+scPreshowerEnergy)*BDTresponse)/(trkMomentum-(scRawEnergy+scPreshowerEnergy)*BDTresponse)"

doCombine = {
    "GedGsfElectron":True,
    "PFMustache":False,
}


clusterings = [
    ("PFMustache","mustacheSCTree/SuperClusterTree"),
    ("GedGsfElectron","gedGsfElectronTree/RegressionTree"),
]

for tagClustering,treeClustering in clusterings:
    # define regressions
    batch.addRegression(tagClustering)
    batch.setInputTree(tagClustering, treeClustering)
    #batch.addVariablesEB(tagClustering, specificVariables[tagClustering], float)
    #batch.addVariablesEE(tagClustering, specificVariablesEE[tagClustering])
    batch.doCombination(tagClustering, doCombine[tagClustering])

#batch.simulate = True

batch.info = info

