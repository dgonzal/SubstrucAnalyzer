import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")


## check the event content 
process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.demo = cms.EDAnalyzer('SimHitResponse',
                              ProducerModule = cms.string("g4SimHits"), #options at the moment "g4SimHits" or "famosSimHits"
                              #HitLabel = cms.InputTag("HcalHits"),#HcalHits,EcalHitsEB
                              OutputName = cms.string("0T.root") 
                              
)

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring( ('keep *')),
    fileName = cms.untracked.string('fastsim.root'),
    eventAutoFlushCompressedSize = cms.untracked.int32(1048576),
    splitLevel = cms.untracked.int32(0),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RECO')
    )
)    



fastsim = cms.untracked.vstring()
fastsim.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/fast/gensim.root'])


fullsim = cms.untracked.vstring()
fullsim.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/full/gensim.root'])

Mag_0T = cms.untracked.vstring()
Mag_0T.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/full/0T_gensim.root'])

Mag_0T_eta = cms.untracked.vstring()
Mag_0T_eta.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/full/0T_eta_gensim.root'])

Mag_4T = cms.untracked.vstring()
Mag_4T.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/full/4T_gensim.root'])

Mag_4T_eta = cms.untracked.vstring()
Mag_4T_eta.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/full/4T_eta_gensim.root'])

test =  cms.untracked.vstring()
test.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_2_0_pre4/src/GeneratorInterface/Pythia8Interface/full/gensim.root'])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source ("PoolSource",fileNames = Mag_0T)

process.p = cms.Path(process.demo)


#write out all the trees in one file, was used for testing purposes
#process.output_step = cms.EndPath(process.output)
#process.schedule = cms.Schedule(*[process.p,process.output_step])


