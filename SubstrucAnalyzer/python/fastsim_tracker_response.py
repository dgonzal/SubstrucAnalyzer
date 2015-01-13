#using grid-control to run. information can be found on: 
#https://ekptrac.physik.uni-karlsruhe.de/trac/grid-control
#
import os
import math

# 1. Invent some configuration variables and set defaults for all of them here; use values suitable for local testing.
# They will be overwritten below in case this is run on the grid via grid-control.


NEvents = 10000
pdgId = 211
energy = 9
minEta = 0.0
deltaEta = 0.1

histoBinning = [5000,0,500]

string = ""


# 2. Now the grid-control interface: let grid-control replace some special variables and transform them into the global ones we just defined above.
# Put them in quotes such that if not using grid-control, it is not a python syntax error.

gc_nickname = os.getenv('MY_JOBID')

if gc_nickname is not None: # This means that we run on the grid with grid-control. 
    NEvents = int('__NEvents__')
    pdgId = int('__pdgId__')
    energy = float('__Energy__')
    minEta = float('__minEta__')

if energy>999:
	histoBinning = [5000,500,4500] 

maxEta = deltaEta + minEta 
# 3. now continue with the usual configuration. At some point, you will make use of the variables, e.g. like this:

print "******************"
print "number of events",NEvents,"PdgId", pdgId, "energy",energy, "min eta",minEta,"max eta", maxEta
print "response binning", histoBinning
print "******************"

import FWCore.ParameterSet.Config as cms
    
process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
        

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEvents)
    )

# Input source
process.source = cms.Source("EmptySource")

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
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = False
process.Realistic8TeVCollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Realistic8TeVCollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Realistic8TeVCollisionVtxSmearingParameters
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')



"""
# Make the tracker transparent (unuseful?)
process.famosSimHits.MaterialEffects.PairProduction = False
process.famosSimHits.MaterialEffects.Bremsstrahlung = False
process.famosSimHits.MaterialEffects.EnergyLoss = False
process.famosSimHits.MaterialEffects.MultipleScattering = False
process.famosSimHits.MaterialEffects.NuclearInteraction = False
"""


process.generator = cms.EDProducer("FlatRandomEGunProducer",
                                   PGunParameters = cms.PSet(
        MinE = cms.double(energy),
        MaxE = cms.double(energy),
        MinEta = cms.double(minEta),
        MaxEta = cms.double(maxEta),
        MinPhi = cms.double(-3.14159265359),#-3.14159265359
        MaxPhi = cms.double(3.14159265359),
        PartID = cms.vint32(pdgId)
        ),
                                   #pythiaHepMCVerbosity = cms.untracked.bool(True),
                                   maxEventsToPrint = cms.untracked.int32(1),
                                   AddAntiParticle = cms.bool(False)
                                   #pythiaPylistVerbosity = cms.untracked.int32(1),
                                   #PythiaParameters = cms.PSet(
                                   #    parameterSets = cms.vstring()
                                   #)
                                   )

process.responseAnalyzer = cms.EDAnalyzer('SimHitFit',
                                          ProducerModule = cms.string("famosSimHits"), #options at the moment "g4SimHits" or "famosSimHits"
                                          OutputName = cms.string(string+"response_pdgId_"+str(pdgId)+"_eta_"+str(minEta)+"_"+str(maxEta)+"_E_"+str(energy)+".root"),
                                          ResponseBinning = cms.vdouble(histoBinning[0],histoBinning[1],histoBinning[2])
                                          )

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("ekin_fast_histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen_genonly)
process.simulation_step = cms.Path(process.simulationWithFamos)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

process.response_step = cms.Path(process.responseAnalyzer)
        
# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.response_step,process.endjob_step,process.FEVTDEBUGoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
    
# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs

"""
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()
"""

process = customisePostLS1(process)

# End of customisation functions
        
