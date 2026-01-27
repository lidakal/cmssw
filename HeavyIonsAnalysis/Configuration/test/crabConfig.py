from WMCore.Configuration import Configuration
config = Configuration()


### General ###s
config.section_('General')
config.General.requestName = 'GNucleus_QCD_Pthat5_jtPtMin_5_genPtMin_3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#### JobType ####
config.section_('JobType')
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "forest_miniAOD_run3_ppRECO_MC_upcJet.py"
# config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
# config.JobType.numCores = 8
# config.JobType.sendExternalFolder = True # To load libxgboost.so
# config.JobType.outputFiles = ['outfileName.root']

#### Data ####
config.section_("Data")
config.Data.inputDataset = '/GNucleus-QCD_Pthat5_5p36TeV_pythia8/HINPbPbSpring23MiniAOD-NoPU_UPC_UPC_132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM'
# config.Data.lumiMask = 
config.Data.inputDBS = "global"
# config.Data.inputDBS = "phys03"
config.Data.splitting = "FileBased"
# config.Data.splitting = "LumiBased"
# config.Data.splitting = "Automatic"
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
# config.Data.outLFNDirBase = '/store/user/lkalipol/bJet2023'
# config.Data.outputDatasetTag = config.General.requestName
# config.Data.ignoreLocality = True
# config.Data.allowNonValidInputDataset = True

config.section_('User')
config.section_('Site')
# config.Site.whitelist = ['T3_KR_KNU']
config.Site.storageSite = 'T3_KR_KNU'