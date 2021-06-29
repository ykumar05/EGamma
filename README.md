# EGamma
Samples used: DoubleElectron gun from 1.7 GeV to 300 GeV in pT

CMSSW version: CMSSW_11_0_1

1.cmsrel CMSSW_11_0_1

2.cd CMSSW_11_0_1/src

3.cmsenv

4.The gen-sim-raw samples with nominal configuration (please use config.Data.inputDBS = 'phys03'):

/SingleEle_V0_GENSIMRAW/shilpi-SingleEle_V0_miniAOD_nominal-28ea11753f803e89e7fd3e375196acc9/USER

Ntuples are already made. Can be found here: /eos/cms/store/group/phys_egamma/shilpi/improveTrackingHighEleEta/ntuples/SingleEle_V0_GENSIMRAW/crab_nominal/210317_004029/0000/

5.For modified config, the gen-sim-raw samples can be used from here: 
5.1: /SingleEle_V0_GENSIMRAW/shilpi-SingleEle_V2_miniAOD_modified_Kensaddition-28ea11753f803e89e7fd3e375196acc9/USER 

Parameters Modified:

process.TrajectoryFilterForElectrons.minimumNumberOfHits = cms.int32(2) ### 5(Original Value)                                                     
process.GsfElectronFittingSmoother.MinNumberOfHits = cms.int32(2) #### 5(Original Value)

5.2: /SingleEle_V0_GENSIMRAW/shilpi-SingleEle_V2_miniAOD_modified_Kensaddition_withMaxEtaAsWell-28ea11753f803e89e7fd3e375196acc9/USER

Parameters Modified:

process.TrajectoryFilterForElectrons.minimumNumberOfHits = cms.int32(2) ### 5(Original)                                                     
process.GsfElectronFittingSmoother.MinNumberOfHits = cms.int32(2) #### 5(Original)                                                         
process.trackerDrivenElectronSeeds.MaxEta = cms.double(3) ### 2.4(Original)

6.Once the above aod config is finalised, generate the aod samples.

To run over a large sample, submit CRAB Jobs. The configuration files used to submit CRAB Jobs are present in Configuration_CRAB_Files.

7.Generate the miniAOD using the above aod samples

8.Then run the ggNtuplizer.

9.Once the modified config ntuples are there, go to analysis directory, update files.list and do:

root -l -b -q runAll.C+
