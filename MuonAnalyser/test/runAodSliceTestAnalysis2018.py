import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('SliceTestAnalysis',eras.Run2_2018,eras.run3_GEM)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v1', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#process.maxEvents.input = cms.untracked.int32(10)
# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.skipEvents = cms.untracked.uint32(0)

# filename = "FCAE5CCD-C29B-E811-8CD3-FA163EF4E3E9"
# process.source.fileNames.append('file:step3_%s.root' % filename)

process.source.fileNames.append("file:/xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv2/step3_139.root")
#process.source.fileNames.append("file:FCAE5CCD-C29B-E811-8CD3-FA163EF4E3E9.root")

process.options = cms.untracked.PSet()

process.TFileService = cms.Service("TFileService",fileName = cms.string("histo_mini.root"))

process.SliceTestAnalysis = cms.EDAnalyzer('AodSliceTestAnalysis',
    process.MuonServiceProxy,
    gemRecHits = cms.InputTag("gemRecHits"),
    muons = cms.InputTag("muons"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
)
process.p = cms.Path(process.SliceTestAnalysis)
