import FWCore.ParameterSet.Config as cms

process = cms.Process("EcalSelectiveReadoutValid")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
###process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

#Number of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#Input files
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/1cea400a-8fc0-42f5-92b8-5c8463fb5774.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/33e42cb1-7143-45ff-bb63-58debfc87d2d.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/462aa48f-1067-4c71-97ea-3a1ce447d1b7.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/1362c4e5-71cd-4f5f-bf5e-20c7f14e6783.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/878cb5ee-47ba-4c14-8ee9-25b44b4237cc.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/5c3df509-7ee5-4db6-87e5-c8387c6cb7b2.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/b7797f91-2a45-4d23-bf35-05f1a62caeb5.root",
          "file:/eos/user/c/cmcginn/HIEGReadoutTesting/RunData/362315/72279191-6fb7-42a5-a72a-1b9cd245d5df.root"
        )
)

# initialize  MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50

# initialize magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

# Added by Stathes
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

# DQM services
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
#process.dqmSaver.workflow = "/HIRun2018A/HIHardProbes/RAW"
process.dqmSaver.workflow = "/HIRun2022A/HITestRaw0/RAW"
#process.load("CalibCalorimetry.Configuration.Ecal_FakeConditions_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##process.load("DQM.Integration.config.FrontierCondition_GT_Offline_cfi")
process.load("DQM.Integration.config.FrontierCondition_GT_cfi")

# Local database for emulator setup (to change ZS thresholds!)
#----------------------------------------------------------------------
#To overwrite Selective readout settings with settings from a local DB (sqlite file): 
#process.load("CondCore.CondDB.CondDB_cfi")
#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("EcalSRSettingsRcd"),
#             tag = cms.string("EcalSRSettings_EB8_EE8"),
#             connect = cms.string("sqlite_file:EcalSRSettings_EB8_EE8.db")
#         )
#)

#----------------------------------------------------------------------

# Sequence info:
# 0. ecalEBunpacker runs: 
#    - it will provide the Trigger Prims (i.e. LI, HI trigger towers, their Et etc).
#    - the name of the output is "EcalTriggerPrimitives"
#    - if input file is data, then the TPs are data
# 1. ecalDigiSequence defined but only "simEcalDigis" runs:
# ecalDigiSequence = cms.Sequence(simEcalTriggerPrimitiveDigis*simEcalDigis*simEcalPreshowerDigis)
#    - Here the simEcalDigis, reads EB and EE digis and TPs from the unpacker
#    - Output names are "egDigis" and "eeDigis"
#    - simEcalDigis also output the "ebSrFlags" and "eeSrFlags"
#    - The user may invoke emulation (EcalSelectiveReadoutSuppressor.cc) uncommenting
#    the flags at the end of this macro.
#    The difference is that the simEcalDigis re-does zero suppression from scratch.
#    It uses the input TTPrimitives getting from there the suppressed Et(), but then
#    it redefines LI,MI,HI etc. Of course I would hope nothing changes as far as the 
#    TTPrims are concerned, but the ZS being repeated should cut more data and 
#    reduce the DCC load!
# 2. ecalSelectiveReadoutValidation runs:
#    Its input can be the output of the unpacker, or the output of the simEcalDigis.
#    - It simply makes a number of SR validation histograms based on the input digis.
#    - If we want to plot the emulated ADC FIR histograms, these histograms are calculated 
#      in the validation code by applying the FIR weights on ADC samples coming from the 
#      input digis. So the ADC FIR histograms depend on the input ebDigis eeDigis.
#

# ECAL Unpacker:
process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi")
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")
process.ecalEBunpacker.silentMode = cms.untracked.bool(True)
process.ecalEBunpacker.InputLabel = cms.InputTag("rawDataRepacker")

#
# ECAL digitization sequence
#
process.load("SimCalorimetry.Configuration.ecalDigiSequence_cff")
process.simEcalDigis.trigPrimProducer = cms.string('ecalEBunpacker')# Label name of input ECAL trigger primitive collection
process.simEcalDigis.trigPrimCollection =  cms.string('EcalTriggerPrimitives')# Input trigger prim collection
process.simEcalDigis.digiProducer = cms.string('ecalEBunpacker')#Input Label of EB EE digi collections
process.simEcalDigis.EBdigiCollection = cms.string('ebDigis')#Input Digis
process.simEcalDigis.EEdigiCollection = cms.string('eeDigis')#Input Digis
process.simEcalDigis.EBSRPdigiCollection = cms.string('SRebDigis')#output Digis (emulated)
process.simEcalDigis.EESRPdigiCollection = cms.string('SReeDigis')#output Digis (emulated)
#Note that the output flags have names (default): ebSrFlags,eeSrFlags


#
# Ecal selective readout validation module, ecalSelectiveReadoutValidation:
#
process.load("Validation.EcalDigis.ecalSelectiveReadoutValidation_cfi")
process.ecalSelectiveReadoutValidation.outputFile = 'srvalid_Run362315_DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0.root'
process.ecalSelectiveReadoutValidation.verbose = cms.untracked.bool(False);

#INPUTS:
process.ecalSelectiveReadoutValidation.EbUnsuppressedDigiCollection = cms.InputTag("simEcalUnsuppressedDigis")
process.ecalSelectiveReadoutValidation.EeUnsuppressedDigiCollection = cms.InputTag("simEcalUnsuppressedDigis")

#Stathes changes:
#DQM: here we use straight the data (no emulator)
#     the output of the unpacker will be used by the ecalSelectiveReadoutValidation.
#process.ecalSelectiveReadoutValidation.EbDigiCollection = cms.InputTag("ecalEBunpacker","ebDigis")#input digis to SRValid
#process.ecalSelectiveReadoutValidation.EeDigiCollection = cms.InputTag("ecalEBunpacker","eeDigis")#input digis to SRValid
#process.ecalSelectiveReadoutValidation.EbSrFlagCollection = cms.InputTag("ecalEBunpacker","")#input flags
#process.ecalSelectiveReadoutValidation.EeSrFlagCollection = cms.InputTag("ecalEBunpacker","")#input flags
#process.ecalSelectiveReadoutValidation.TrigPrimCollection = cms.InputTag("ecalEBunpacker", "EcalTriggerPrimitives")
#
# My MC script: (if you turn on the emulator uncomment the following).
# Here we force the simEcalDigis to do the zero suppression and their output will go 
# to the ecalSelectiveReadoutValidation.
process.ecalSelectiveReadoutValidation.EbDigiCollection = cms.InputTag("simEcalDigis","SRebDigis")#input digis to SRValid
process.ecalSelectiveReadoutValidation.EeDigiCollection = cms.InputTag("simEcalDigis","SReeDigis")#input digis to SRValid
process.ecalSelectiveReadoutValidation.EbSrFlagCollection = cms.InputTag("simEcalDigis","ebSrFlags")#input is the output of simEcalDigis
process.ecalSelectiveReadoutValidation.EeSrFlagCollection = cms.InputTag("simEcalDigis","eeSrFlags")#input
process.ecalSelectiveReadoutValidation.TrigPrimCollection = cms.InputTag("simEcalDigis", "simEcalTriggerPrimitives")
# end

#process.ecalSelectiveReadoutValidation.ecalDccZs1stSample = 2
#process.ecalSelectiveReadoutValidation.dccWeights = [ -0.374, -0.374, -0.3629, 0.2721, 0.4681, 0.3707 ]
process.ecalSelectiveReadoutValidation.histDir = ''
process.ecalSelectiveReadoutValidation.histograms = [ 'all' ]
#process.ecalSelectiveReadoutValidation.useEventRate = False
process.ecalSelectiveReadoutValidation.LocalReco = cms.bool(False) # local pulse Ampl reco
process.ecalSelectiveReadoutValidation.IsDATA    = cms.bool(True)
process.ecalSelectiveReadoutValidation.ebZsThrADCCount = 4.5
process.ecalSelectiveReadoutValidation.eeZsThrADCCount = 6.5

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

process.tpparams12 = cms.ESSource("EmptyESSource",
    recordName = cms.string('EcalTPGPhysicsConstRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

#process.load("pgras.ListCollection.ListCollection_cfi")

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.Trigger = hltHighLevel.clone(
    HLTPaths = ["HLT_HIMinimumBias_v*"]
)
process.p1 = cms.Path(process.Trigger*process.ecalEBunpacker*process.simEcalDigis*process.ecalSelectiveReadoutValidation*process.dqmSaver)

#
# switch Modify ECAL SR
#
# to work you need (1) to uncomment the bypass flags below,
# (2) you need to uncomment the 4 lines above (see "turn on the emulator" above)
from SimCalorimetry.EcalSimProducers.ecaldigi_cfi import *
from SimCalorimetry.EcalSelectiveReadoutProducers.ecalDigis_cfi import *
#simEcalDigis.trigPrimBypass = cms.bool(True) # uncomment to bypass
#simEcalDigis.trigPrimBypassMode = cms.int32(1) #bypass mode (uncomment)
#simEcalDigis.trigPrimBypassLTH = cms.double(4.0)# 2xGeV C1 DEFAULT
#simEcalDigis.trigPrimBypassHTH = cms.double(8.0)# 2xGeV C1 DEFAULT
#simEcalDigis.trigPrimBypassLTH = cms.double(4.0)# 2xGeV C1 //2GeV Et 2018 Default
#simEcalDigis.trigPrimBypassHTH = cms.double(8.0)# 2xGeV C1 //4GeV Et 2018 Default
#simEcalDigis.trigPrimBypassLTH = cms.double(0.0)# 2xGeV C1
#simEcalDigis.trigPrimBypassHTH = cms.double(0.0)# 2xGeV C1
#switch to apply selective readout decision on the digis and produce
simEcalDigis.produceDigis = cms.untracked.bool(True)


process.CaloGeometryBuilder.SelectedCalos = cms.vstring(
        'HCAL',
        'EcalBarrel',
        'EcalEndcap',
        'EcalPreshower',
        'TOWER'
    )

process.TFileService = cms.Service("TFileService",
                                fileName = cms.string("srvalid_Run362315_DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0.root"))

process.dqmSaver.saveByLumiSection = cms.untracked.int32(1)
process.dqmSaver.saveByRun = cms.untracked.int32(1)

simEcalDigis.trigPrimBypass = cms.bool(True)
simEcalDigis.trigPrimBypassMode = cms.int32(1)
simEcalDigis.trigPrimBypassLTH = cms.double(4.0)
simEcalDigis.trigPrimBypassHTH = cms.double(8.0)


process.GlobalTag.toGet = cms.VPSet(
 cms.PSet(record = cms.string("EcalSRSettingsRcd"),
   tag = cms.string("EcalSRSettings_DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0"),
   connect = cms.string("sqlite_file:EcalSRSettings_DoSR_LTH4p0_HTH8p0_NoZS_MIEB8p0_MIEE8p0.db")
 )
)

