import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

# setup 'standard'  options
options = VarParsing.VarParsing ('standard')

# setup any defaults you want
options.output = "GenAnalq_M2000_Lambda0p1.root"
#options.files= 'file:/afs/cern.ch/work/s/santanas/Workspace/CMS/CMSSW_10_6_28_LQSim/src/TT_13TeV_Pow_Herwig7_cff_py_LHE_GEN_inRAWSIM.root'
options.files= 'root://eoscms///eos/cms/store/group/phys_exotica/lq-LQ-lq/test_GEN_1/singleLQ_13TeV_Pow_Herwig7_M2000_Lambda0p1_GEN.root'
#options.files= 'root://eoscms///eos/cms/store/group/phys_exotica/lq-LQ-lq/test_GEN_1/singleLQ_13TeV_Pow_Herwig7_M2000_Lambda0p1_GEN.root'
#options.files= 'root://eoscms///eos/cms/store/cmst3/user/santanas/MCsamples/Res1ToRes2GluTo3Glu/Res1ToRes2GluTo3Glu_M1-500_R-0p1_slc6_amd64_gcc481_CMSSW_7_1_30_GEN.root'
#options.files= 'root://eoscms///eos/cms/store/cmst3/user/santanas/MCsamples/TrijetRes_g_ggg_BP2_testV2/TrijetRes_g_ggg_BP2_testV2_MGKK2000R0p7_slc6_amd64_gcc481_CMSSW_7_1_30_GEN.root'
options.maxEvents = -1
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        options.files        
    ),
    #                         firstEvent = cms.untracked.uint32(908)
)

process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string(options.output),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

process.genAnalyzer = cms.EDAnalyzer('GenAnalq',
                                     genParticlesInputTag = cms.InputTag('genParticles'),
                                     genJetsAK4InputTag = cms.InputTag('ak4GenJets'),
                                     genJetsAK8InputTag = cms.InputTag('ak8GenJets')
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(3),
                                   printVertex = cms.untracked.bool(False),
                                   printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.  
                                   src = cms.InputTag("genParticles")
)

process.p = cms.Path(process.printTree+process.genAnalyzer)
