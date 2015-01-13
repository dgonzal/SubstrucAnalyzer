#using grid-control to run. information can be found on: 
#https://ekptrac.physik.uni-karlsruhe.de/trac/grid-control
#
import os
import math
import glob

# 1. Invent some configuration variables and set defaults for all of them here; use values suitable for local testing.
# They will be overwritten below in case this is run on the grid via grid-control.


NEvents = 10000
pdgId = '211'
energy = '20'
minEta = '1.0'
deltaEta = '0.1'

histoBinning = [150,0,40]

string = "newPlotsTimeEMenergy_"

directory="/nfs/dust/cms/user/gonvaq/fastsim/fullsim/no_magneticField/"

# 2. Now the grid-control interface: let grid-control replace some special variables and transform them into the global ones we just defined above.
# Put them in quotes such that if not using grid-control, it is not a python syntax error.

gc_nickname = os.getenv('MY_JOBID')

if gc_nickname is not None: # This means that we run on the grid with grid-control. 
    NEvents = int('__NEvents__')
    pdgId = '__pdgId__'
    energy = '__Energy__'
    minEta = '__minEta__'
    directory = '__directory__'

if float(energy)>999:
    histoBinning = [5000,500,4500]


maxEta = str(float(deltaEta) + float(minEta)) 
# 3. now continue with the usual configuration. At some point, you will make use of the variables, e.g. like this:

print "******************"
print "number of events",NEvents,"PdgId", pdgId, "energy",energy, "min eta",minEta,"max eta", maxEta
print "response binning", histoBinning
print "******************"

source = cms.untracked.vstring()

pattern = directory+"gensim_pdgId_"+str(pdgId)+"_eta_"+str(float(minEta))+"_"+str(float(maxEta))+"_E_"+str(energy)+"."+"*"
patHigh = directory+"high_energy/"+"gensim_pdgId_"+str(pdgId)+"_eta_"+str(float(minEta))+"_"+str(float(maxEta))+"_E_"+str(energy)+"."+"*"

print 'pattern: "%s"' % pattern
print 'pattern high energies: "%s"' % patHigh
print glob.glob(pattern) 

for name in glob.glob(pattern):
    source.extend( ['file:'+name])

for name in glob.glob(patHigh):
    source.extend( ['file:'+name])

print source

import FWCore.ParameterSet.Config as cms
    
process = cms.Process('RESPONSE')

        
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

#process.load('Configuration.Geometry.GeometryRecoECALHCAL_cff')
#process.load('Configuration.Geometry.GeometrySimECALHCAL_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')

#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEvents)
    )

# Input source
#process.source = cms.Source("EmptySource")
process.source = cms.Source ("PoolSource",fileNames = source)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.options = cms.untracked.PSet(
    
    )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('GeneratorInterface/Pythia8Interface/Py8EGun_cfi nevts:'+str(NEvents)),
    name = cms.untracked.string('Applications')
    )

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
                                           splitLevel = cms.untracked.int32(0),
                                           eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                           outputCommands = process.FEVTDEBUGEventContent.outputCommands,
                                           fileName = cms.untracked.string(string+"gensim_pdgId_"+str(pdgId)+"_eta_"+str(minEta)+"_"+str(maxEta)+"_E_"+str(energy)+".root"),
                                           dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
        ),
                                           SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
        )
                                           )

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.responseAnalyzer = cms.EDAnalyzer('SimHitFit',
                                          ProducerModule = cms.string("g4SimHits"), #options at the moment "g4SimHits" or "famosSimHits"
                                          OutputName = cms.string(string+"response_pdgId_"+str(pdgId)+"_eta_"+str(minEta)+"_"+str(maxEta)+"_E_"+str(energy)+".root"),
                                          ResponseBinning = cms.vdouble(histoBinning[0],histoBinning[1],histoBinning[2])
                                          )

process.response_step = cms.Path(process.responseAnalyzer)


