from WMCore.Configuration import Configuration
config = Configuration()


### General ###
config.section_('General')
config.General.requestName = 'bJet2017G_qcdMC_trackTaggingTraining'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#### JobType ####
config.section_('JobType')
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "runForestAOD_pp_MC_94X.py"
#config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8

#### Data ####
config.section_("Data")

# qcdMC dataset
config.Data.inputDataset = "/QCD_pThat-15_Dijet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV_v1-v1/AODSIM"

# bJetMC dataset
#config.Data.inputDataset = "/QCD_pThat-15_bJet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV_v1-v1/AODSIM"

#config.Data.inputDataset ='/QCD_pThat-15_Mujet_TuneCP5_5p02TeV_pythia8/RunIIpp5Spring18DR-94X_mc2017_realistic_forppRef5TeV_v1-v3/AODSIM'
#config.Data.inputDBS = "global"
#config.Data.inputDBS = "phys03"
config.Data.splitting = "FileBased"
#config.Data.splitting = "LumiBased"
#config.Data.splitting = "Automatic"
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/lkalipol/bJet2022'
config.Data.outputDatasetTag = config.General.requestName
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True

#### Site ####
config.section_('Site')
config.Site.whitelist = ['T2_FR_*', "T2_CH_*", "T2_US_MIT"]
config.Site.storageSite = 'T2_FR_GRIF_LLR'
