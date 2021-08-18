# RegressionTraining

# The main package is Forked from 

```
https://github.com/jbsauvan/RegressionTraining

```

#  Environment setup
```
cmsrel CMSSW_8_2_0
cd CMSSW_8_2_0/src/
cmsenv
git clone https://github.com/cms-egamma/HiggsAnalysis
scram b -j8
```

# Download the RegressionTraining package

```
git clone https://github.com/sandeepbhowmik1/RegressionTraining
cd RegressionTraining
make clean
make

./regression.exe Run3/GBRFullLikelihood_Trigger_Stage2_2021_compressedieta_compressedE_hasEM_isMerged_20210727.config

```