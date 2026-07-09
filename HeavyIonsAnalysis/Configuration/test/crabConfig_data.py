from WMCore.Configuration import Configuration
config = Configuration()


### General ###s
config.section_('General')
config.General.requestName = 'HiForward19_HIRun2023A-16Jan2024-v1_ZDCjson_wCCfilter'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#### JobType ####
config.section_('JobType')
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "forest_miniAOD_run3_forward_rereco_DATA_upcJet.py"
# config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
# config.JobType.numCores = 8
config.JobType.sendExternalFolder = True # To load libxgboost.so
# config.JobType.outputFiles = ['outfileName.root']

#### Data ####
config.section_("Data")
config.Data.inputDataset = '/HIForward19/HIRun2023A-16Jan2024-v1/MINIAOD'
config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Good_ZDC_Golden.json'
config.Data.inputDBS = "global"
# config.Data.inputDBS = "phys03"
# config.Data.splitting = "FileBased"
config.Data.splitting = "LumiBased"
# config.Data.splitting = "Automatic"
config.Data.unitsPerJob = 3
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