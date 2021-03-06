import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('*'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG')
    ),
    destinations = cms.untracked.vstring('cout')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource")

process.load("DQM.SiStripCommon.TkHistoMap_cff")
# load TrackerTopology (needed for TkDetMap and TkHistoMap)
process.load("Geometry.TrackerGeometryBuilder.trackerParameters_cfi")
process.trackerTopology = cms.ESProducer("TrackerTopologyEP")

process.tester = cms.EDAnalyzer("testTkHistoMap",
                              readFromFile = cms.bool(False)
                              )
process.p = cms.Path(process.tester)


