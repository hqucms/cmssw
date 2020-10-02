import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.inputFiles = '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/SinglePhoton_PT2to200/FEVT/PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/250000/FFE433A9-4EA7-A14C-9D9B-ACB9D1E9097F.root'
options.maxEvents = 10
options.parseArguments()

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('TEST',Phase2C9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents),
)

## Source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles)
)

## Options and Output Report
process.options = cms.untracked.PSet(
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4),
    wantSummary = cms.untracked.bool(True),
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

# load LC PU weights producer
process.load('RecoHGCal.TICL.hgcalLayerClusterPileupWeights_cfi')
process.load('RecoHGCal.Configuration.recoHGCAL_cff')
process.iterTICLTask.add(process.hgcalLayerClusterPileupWeights)
process.filteredLayerClustersTrk.LayerClustersInputMask = 'hgcalLayerClusterPileupWeights'

## Output Module Configuration (expects a path 'p')
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('test_ticl_with_lc_puId.root'),
                               outputCommands = cms.untracked.vstring(
                                   'drop *', 
                                   'keep *_hgcalLayerClusters*_*_*',
                                   'keep *_hgcalLayerClusterPileupWeights*_*_*',
                                   'keep *_ticlTrackstersEM_*_*',
                                   'keep *_ticlTrackstersHAD_*_*',
                                   'keep *_ticlTrackstersTrk_*_*',
                                   'keep *_ticlTrackstersMIP_*_*',
                                   'keep *_ticlTrackstersMerge_*_*',
                                   'keep *_ticlCandidateFromTracksters_*_*',
                                   'keep *_pfTICL_*_*',
                                   )
                               )

process.outpath = cms.EndPath(process.out, process.iterTICLTask)
