import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.inputFiles = ''
options.maxEvents = 10
options.parseArguments()

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('PU',Phase2C9)
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles)
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.load('RecoHGCal.TICL.layerClusterPileUpWeightsProducer_cfi')
process.hgcalLayerClusterPileupWeights.debugMode = True
process.p = cms.Path(process.hgcalLayerClusterPileupWeights)

## Output Module Configuration (expects a path 'p')
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('test_lc_pu_id.root'),
                               outputCommands = cms.untracked.vstring(
                                   'drop *', 
                                   'keep *_hgcalLayerClusters*_*_*',
                                   'keep *_hgcalLayerClusterPileupWeights*_*_*',
                                   )
                               )

process.outpath = cms.EndPath(process.out)
