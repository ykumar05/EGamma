#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
#import sys

#config = config()

from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException

from CRABClient.UserUtilities import config
config = config()

from multiprocessing import Process

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException, hte:
        print hte.headers

#**************************submit function***********************
#from CRABAPI.RawCommand import crabCommand
#from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException
#def submit(config):
#        try:
#                crabCommand('submit', config = config)
#        except HTTPException as hte:
#                print "Failed submitting task: %s" % (hte.headers)
#        except ClientException as cle:
#                print "Failed submitting task: %s" % (cle)
#****************************************************************

#submitVersion = "V1_RAWSIM"
#submitVersion = "V0_AOD_modified"
submitVersion = "AOD_Run3"
mainOutputDir = "/store/user/ykumar/%s" %submitVersion


#General
config.General.requestName = 'AOD_Run3_Modified'
config.General.transferLogs = True
config.General.workArea = 'crab_projects_%s' % submitVersion

#Site
config.Site.storageSite = 'T2_IN_TIFR'
config.Site.whitelist = ['T2_*']

#JobType
config.JobType.allowUndistributedCMSSW = True 
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName  = 'aod_cfg.py'
#config.JobType.maxJobRuntimeMin = '' 
config.JobType.numCores = 4
config.JobType.maxMemoryMB = 2500

#Data
config.Data.inputDataset = '/SingleEle_V0_GENSIMRAW/shilpi-crab_SingleEle_GENSIMRAW-2133ffc961f67f7c3fa0060c152fa455/USER'
config.Data.inputDBS = 'phys03'
config.Data.publication = True
#config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '%s' % (mainOutputDir)
#config.Data.splitting    = 'EventBased'
#config.Data.splitting    = 'Automatic'
config.Data.splitting     = 'FileBased'
config.Data.totalUnits    = 1000
config.Data.unitsPerJob   = 10
#config.Data.unitsPerJob  = 1000
#config.Data.totalUnits   = 10000000
config.Data.outputDatasetTag ='SingleEle_%s' %submitVersion
submit(config)
