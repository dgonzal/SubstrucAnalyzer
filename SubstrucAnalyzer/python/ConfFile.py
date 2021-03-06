import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#/store/relval/CMSSW_7_1_0_pre7/RelValQCD_FlatPt_15_3000HS/GEN-SIM-RECO/PRE_STA71_V3-v1/00000/0C74A8E9-40D1-E311-A1D8-0025905A6056.root
#/store/relval/CMSSW_7_1_0_pre7/RelValQCD_FlatPt_15_3000HS/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_STA71_V3-v1/00000/1C32DED0-2DD1-E311-989D-003048678FFA.root
#/store/relval/CMSSW_7_1_0_pre7/RelValTTbar/GEN-SIM-RECO/PRE_STA71_V3-v1/00000/887DEA5B-5CD1-E311-BB97-002618943923.root
#'/store/relval/CMSSW_7_1_0_pre7/RelValTTbar/GEN-SIM-RECO/PRE_STA71_V3-v1/00000/887DEA5B-5CD1-E311-BB97-002618943923.root',
#'/store/relval/CMSSW_7_1_0_pre7/RelValTTbar/GEN-SIM-RECO/PRE_STA71_V3-v1/00000/8EEC0F1F-9FD1-E311-966F-003048FFD76E.root' ] );

ttbar_File = cms.untracked.vstring()
ttbar_File.extend( [
        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_7_1_0_pre7/RelValTTbar/GEN-SIM-RECO/PRE_STA71_V3-v1/00000/887DEA5B-5CD1-E311-BB97-002618943923.root',
        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_7_1_0_pre7/RelValTTbar/GEN-SIM-RECO/PRE_STA71_V3-v1/00000/8EEC0F1F-9FD1-E311-966F-003048FFD76E.root'
        ])

ttbar_fs = cms.untracked.vstring()

ttbar_fs.extend( [
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/000E45D8-AFD0-E311-B051-02163E00EA10.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/069135E9-B1D0-E311-8DE3-02163E00EA8F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/0CCF32A6-9DD2-E311-B64C-02163E00F535.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/10E0AB48-B0D0-E311-831B-02163E00E74E.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/1214DB57-B4D0-E311-90C4-02163E00B786.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/123EE7D4-B0D0-E311-8F32-02163E00B1AF.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/1414C636-B1D0-E311-84AD-02163E00E868.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/1666EC72-B4D0-E311-8B05-02163E00ADE5.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/16C641A9-AFD0-E311-AD55-02163E00EB06.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/18F4093C-B1D0-E311-A35D-02163E00E7BE.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/2402F011-B3D0-E311-853F-02163E00B55D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/2624E65F-AFD0-E311-B352-02163E00EA93.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/26E04D1B-C1D0-E311-A8F8-02163E00E7E8.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/28160F85-FED0-E311-AF56-02163E00CDFC.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/28D1143B-B0D0-E311-9F2F-02163E00E999.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/2A47BA7F-B2D0-E311-92AF-02163E00AE81.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/2E3DD0D4-B1D0-E311-8CA7-02163E00BF9F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/343671C8-B4D0-E311-A337-02163E00E69D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/362C21A2-B5D0-E311-A6B2-02163E00CDFF.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/3AD77981-B5D0-E311-966E-02163E00EA05.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/3C2DD378-B2D0-E311-8184-02163E00B6F2.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/3E2B9D29-D8D0-E311-9557-02163E00E9BA.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/3E386891-CFD0-E311-A0C5-02163E00E7C7.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/3EA91D0B-B2D0-E311-8DF0-02163E00B507.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/4A38EFED-B2D0-E311-A4C5-02163E00AD8F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/4EE82F8E-AFD0-E311-B5C4-02163E00F32F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/52C6E7AD-AFD0-E311-BED6-02163E00EB68.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/562BC562-AFD0-E311-B330-02163E00E9E8.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/564215BF-AFD0-E311-98B3-02163E00E90A.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/5C24AE95-B7D0-E311-B990-02163E00BA35.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/5CD44015-D2D0-E311-9D75-02163E00F500.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/62E855F7-B1D0-E311-82BC-02163E00AE04.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/66A7CC43-B3D0-E311-84E6-02163E00B7AF.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/6E48E292-AFD0-E311-B624-02163E00EA93.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/7294D111-B5D0-E311-9DFA-02163E00CDDB.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/74038C2D-B4D0-E311-9972-02163E00EA11.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/746FA29B-AFD0-E311-B62C-02163E00EBBF.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/764BDAA9-B3D0-E311-991F-02163E00C399.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/76F85AA3-B4D0-E311-9E2C-02163E00C76F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/7A94F07C-AFD0-E311-A006-02163E00E8B0.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/7E427187-B3D0-E311-BA97-02163E00ADA4.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/7E4D4F93-B4D0-E311-A2E0-02163E00EA47.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/8252DF7A-B1D0-E311-8816-02163E00BA13.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/86D1E022-B2D0-E311-9F43-02163E00E92F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/8C93E65C-B4D0-E311-B5C0-02163E00C5B4.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/8E5FBEC5-AFD0-E311-9F2B-02163E00E8A9.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/907AFAA1-B1D0-E311-96DB-02163E00EB35.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/9422AAC9-B0D0-E311-ABCE-02163E00B575.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/9840D037-B4D0-E311-AEC6-02163E00C61F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/9A678969-AFD0-E311-BAE3-02163E00E705.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/9AA62DCC-AFD0-E311-A3D0-02163E00E770.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/9E40F32B-B5D0-E311-98A5-02163E00BC2E.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/9E91D5C2-C2D0-E311-A8BD-02163E00E760.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/A20D5C8C-B1D0-E311-9593-02163E00BF80.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/A4F66922-D7D0-E311-95B3-02163E00E78D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/A61C7ACF-B0D0-E311-80F6-02163E00BD4D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/A85C231B-B1D0-E311-9779-02163E00EB42.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/AA36069E-B5D0-E311-8364-02163E00BFB7.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/AE44995C-B1D0-E311-9542-02163E00B4F3.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/B019F0D4-B0D0-E311-90AD-02163E00B1AF.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/B01B278E-AFD0-E311-A68C-02163E00E7B6.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/B07FE854-B4D0-E311-9C06-02163E00E79D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/B8D1FAC4-AFD0-E311-8C91-02163E00EA3F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/BCD1147C-B2D0-E311-B4FE-02163E00C569.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/C228FED7-B0D0-E311-B792-02163E00B772.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/C44B2EB6-B4D0-E311-A652-02163E00E88A.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/C492F4AD-AFD0-E311-9B52-02163E00E795.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/C6F6361D-C0D0-E311-9FC6-02163E00EA65.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/CE10FE15-B1D0-E311-B64A-02163E00E72D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D0447BC8-B4D0-E311-9CA2-02163E00E792.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D0CB5FF2-B0D0-E311-BD0F-02163E00BFFD.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D2984E1C-B0D0-E311-9E40-02163E00F4D2.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D2D5FCF2-DED0-E311-9052-02163E00CDA2.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D468C1A8-B1D0-E311-94C8-02163E00E734.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D4A35BBB-BFD0-E311-9774-02163E00EA3F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D691977E-B2D0-E311-B68D-02163E00AE81.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D6C95008-B5D0-E311-AE35-02163E00EB5F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/D83429B7-B1D0-E311-BA48-02163E00CE27.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/DC1943A7-AFD0-E311-A7B2-02163E00E84A.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/DCE830AD-B5D0-E311-9CD7-02163E00B77A.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/DE1E62AA-D6D0-E311-9389-02163E00E7BF.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/DEFA0EF2-B0D0-E311-A6F5-02163E00B742.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/E063E216-B2D0-E311-9114-02163E00C8A7.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/E2D5C93E-B1D0-E311-AFAE-02163E00CAC0.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/E4E83FFE-C0D0-E311-9E0C-02163E00BA1D.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/EACD8374-B2D0-E311-97B7-02163E00C4DA.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/F2141290-AFD0-E311-8417-02163E00E7B6.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/F2390D61-AFD0-E311-9783-02163E00F32F.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/F24812D3-B0D0-E311-927C-02163E00EB68.root',
       '/store/relval/CMSSW_7_1_0_pre7/RelValTTbar_13/GEN-SIM-DIGI-RECO/PRE_LS171_V7_FastSim-v1/00000/F8BE024F-AFD0-E311-8EB3-02163E00EBBC.root' ] );

ttbar_fs_File = cms.untracked.vstring()
for s in ttbar_fs:
  ttbar_fs_File.extend(['root://cms-xrd-global.cern.ch/'+ s])

process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("genParticles"),
    ignoreParticleIDs = cms.vuint32(
         1000022,
         1000012, 1000014, 1000016,
         2000012, 2000014, 2000016,
         1000039, 5100039,
         4000012, 4000014, 4000016,
         9900012, 9900014, 9900016,
         39,
         12, 14, 16),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(True),
    excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
    tausAsJets = cms.bool(False)
)

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak8GenJets = ak5GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJets"))


process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string('!isFake && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0')
)

process.particleFlowPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
    src = cms.InputTag("particleFlow")
)

process.pfPileUpJME = cms.EDProducer("PFPileUp",
    checkClosestZVertex = cms.bool(False),
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlowPtrs"),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("goodOfflinePrimaryVertices")
)

process.pfNoPileUpJME = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUpJME"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)

process.ak8PFJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('PFJet'),
    src = cms.InputTag("particleFlow"),
    doPUOffsetCorr = cms.bool(False),
    radiusPU = cms.double(0.5),
    doAreaDiskApprox = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8)
)

process.ak8PFJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('PFJet'),
    src = cms.InputTag("particleFlow"),
    doPUOffsetCorr = cms.bool(False),
    radiusPU = cms.double(0.5),
    doAreaDiskApprox = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8)
)



process.ak8PFJetsConstituents = cms.EDFilter("PFJetConstituentSelector",
    src = cms.InputTag("ak8PFJets"),
    cut = cms.string('pt > 100.0 && abs(rapidity()) < 2.4')
)

process.ak8PFJetsCHSConstituents = cms.EDFilter("PFJetConstituentSelector",
    src = cms.InputTag("ak8PFJetsCHS"),
    cut = cms.string('pt > 100.0 && abs(rapidity()) < 2.4')
)

process.cmsTopTagPFJets = cms.EDProducer("CATopJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    voronoiRfact = cms.double(-0.9),
    Ghost_EtaMax = cms.double(5.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doPUOffsetCorr = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    src = cms.InputTag("ak8PFJetsConstituents","constituents"),
    inputEMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    jetPtMin = cms.double(100.0),
    doAreaDiskApprox = cms.bool(False),
    radiusPU = cms.double(0.5),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    sumEtEtaCut = cms.double(3.0),
    ptFracBins = cms.vdouble(0.05, 0.05, 0.05),
    rBins = cms.vdouble(0.8, 0.8, 0.8),
    tagAlgo = cms.int32(1),
    etFrac = cms.double(0.7),
    useMaxTower = cms.bool(False),
    deltarBins = cms.vdouble(0.19, 0.19, 0.19),
    nCellBins = cms.vdouble(1.9, 1.9, 1.9),
    debugLevel = cms.untracked.int32(0),
    sumEtBins = cms.vdouble(0, 1600, 2600),
    centralEtaCut = cms.double(2.5),
    useAdjacency = cms.int32(2),
    algorithm = cms.int32(1),
    jetCollInstanceName = cms.string('cmsTopSubJets'),
    verbose = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    rParam = cms.double(0.8),
    ptFrac = cms.double(0.05),
    rFrac = cms.double(0.01),
    adjacencyParam  = cms.double(0.00001),
    writeCompound = cms.bool(True)
)

process.cmsTopTaggerPFJets = cms.EDProducer("CATopJetTagger",
                                         src  = cms.InputTag("cmsTopTagPFJets"),
                                         TopMass = cms.double(171.),
                                         WMass = cms.double(80.4),
                                         TopMassMin = cms.double(0.),
                                         TopMassMax = cms.double(250.),
                                         WMassMin = cms.double(0.0),
                                         WMassMax = cms.double(200.0),
                                         MinMassMin = cms.double(0.0),
                                         MinMassMax = cms.double(200.0),
                                         verbose = cms.bool(False)
)

process.cmsTopTagPFJetsCHS = cms.EDProducer("CATopJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    voronoiRfact = cms.double(-0.9),
    Ghost_EtaMax = cms.double(5.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doPUOffsetCorr = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    src = cms.InputTag("ak8PFJetsCHSConstituents","constituents"),
    inputEMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    jetPtMin = cms.double(100.0),
    doAreaDiskApprox = cms.bool(False),
    radiusPU = cms.double(0.5),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    sumEtEtaCut = cms.double(3.0),
    ptFracBins = cms.vdouble(0.05, 0.05, 0.05),
    rBins = cms.vdouble(0.8, 0.8, 0.8),
    tagAlgo = cms.int32(1),
    etFrac = cms.double(0.7),
    useMaxTower = cms.bool(False),
    deltarBins = cms.vdouble(0.19, 0.19, 0.19),
    nCellBins = cms.vdouble(1.9, 1.9, 1.9),
    debugLevel = cms.untracked.int32(0),
    sumEtBins = cms.vdouble(0, 1600, 2600),
    centralEtaCut = cms.double(2.5),
    useAdjacency = cms.int32(2),
    algorithm = cms.int32(1),
    jetCollInstanceName = cms.string('cmsTopSubJetsCHS'),
    verbose = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    rParam = cms.double(0.8),
    ptFrac = cms.double(0.05),
    rFrac = cms.double(0.01),
    adjacencyParam  = cms.double(0.00001),
    writeCompound = cms.bool(True)
)

process.NSubjettinesPFJets = cms.EDProducer("CATopJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    voronoiRfact = cms.double(-0.9),
    Ghost_EtaMax = cms.double(5.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doPUOffsetCorr = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    src = cms.InputTag("ak8PFJetsConstituents","constituents"),
    inputEMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    jetPtMin = cms.double(100.0),
    doAreaDiskApprox = cms.bool(False),
    radiusPU = cms.double(0.5),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    sumEtEtaCut = cms.double(3.0),
    ptFracBins = cms.vdouble(0.05, 0.05, 0.05),
    rBins = cms.vdouble(0.8, 0.8, 0.8),
    tagAlgo = cms.int32(4),
    etFrac = cms.double(0.7),
    useMaxTower = cms.bool(False),
    deltarBins = cms.vdouble(0.19, 0.19, 0.19),
    nCellBins = cms.vdouble(1.9, 1.9, 1.9),
    debugLevel = cms.untracked.int32(0),
    sumEtBins = cms.vdouble(0, 1600, 2600),
    centralEtaCut = cms.double(2.5),
    useAdjacency = cms.int32(2),
    jetCollInstanceName = cms.string('NSubJets'),
    verbose = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True),
    tau2Cut = cms.double(0.08),
    cosThetaSCut = cms.double(0.8),
    useExclusive = cms.bool(False)
)

process.hepTopTagPFJets = cms.EDProducer("CATopJetProducer",
    muCut = cms.double(10),
    maxSubjetMass = cms.double(50),
    useSubjetMass = cms.bool(True),
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    voronoiRfact = cms.double(-0.9),
    Ghost_EtaMax = cms.double(5.0),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doPUOffsetCorr = cms.bool(False),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    src = cms.InputTag("ak8PFJetsConstituents","constituents"),
    inputEMin = cms.double(0.0),
    jetType = cms.string('PFJet'),
    jetPtMin = cms.double(100.0),
    doAreaDiskApprox = cms.bool(False),
    radiusPU = cms.double(0.5),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    sumEtEtaCut = cms.double(3.0),
    ptFracBins = cms.vdouble(0.05, 0.05, 0.05),
    rBins = cms.vdouble(0.8, 0.8, 0.8),
    tagAlgo = cms.int32(2),
    etFrac = cms.double(0.7),
    useMaxTower = cms.bool(False),
    deltarBins = cms.vdouble(0.19, 0.19, 0.19),
    nCellBins = cms.vdouble(1.9, 1.9, 1.9),
    debugLevel = cms.untracked.int32(0),
    sumEtBins = cms.vdouble(0, 1600, 2600),
    centralEtaCut = cms.double(2.5),
    useAdjacency = cms.int32(2),
    jetCollInstanceName = cms.string('hepTopSubJets'),
    verbose = cms.bool(False),
    jetAlgorithm = cms.string('CambridgeAachen'),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
)

process.ak4PFJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('PFJet'),
    src = cms.InputTag("particleFlow"),
    doPUOffsetCorr = cms.bool(False),
    radiusPU = cms.double(0.5),
    doAreaDiskApprox = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.4)
)
process.ak4PFJetsCHS = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(False),
    srcPVs = cms.InputTag(""),
    doAreaFastjet = cms.bool(True),
    maxBadHcalCells = cms.uint32(9999999),
    doAreaDiskApprox = cms.bool(False),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('PFJet'),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    maxBadEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.4),
    maxProblematicHcalCells = cms.uint32(9999999),
    src = cms.InputTag("pfNoPileUpJME"),
    inputEtMin = cms.double(0.0),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)

process.demo = cms.EDAnalyzer('SubstrucAnalyzer',
                              jetPtMins = cms.vdouble(100.0, 100.0, 100.0, 100.0, 100.0,100.0,100.0,100.0),
                              matching_radius = cms.double(1.),
                              GenParticle = cms.InputTag("genParticles"),
                              OutputName = cms.string("test.root"),
                              jetLabels = cms.VInputTag("ak8PFJets","ak8GenJets","ak8PFJetsConstituents","cmsTopTagPFJets","hepTopTagPFJets","NSubjettinesPFJets")
)

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring( ('keep *')),
    fileName = cms.untracked.string('result_all.root'),
    eventAutoFlushCompressedSize = cms.untracked.int32(1048576),
    splitLevel = cms.untracked.int32(0),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RECO')
    )
)    


fastsim_low = cms.untracked.vstring()
fastsim_low.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/fast/pgun_tt_0_300.root',
                     'file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/fast/pgun_tt_300_800.root',
                     'file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/fast/pgun_tt_800_1500.root',
                     'file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/fast/pgun_tt_1500_2000.root'])
#test.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/full/reco.root'])

fastsim1 = cms.untracked.vstring()
fastsim1.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/fast/pgun_tt.root'])


fullsim1 = cms.untracked.vstring()
fullsim1.extend( ['file:/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_7_1_0_pre9/src/GeneratorInterface/Pythia8Interface/full/reco.root'])


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source ("PoolSource",fileNames = fullsim1)

process.p = cms.Path(process.genParticlesForJets+process.ak8GenJets+process.ak8PFJets+process.ak8PFJetsConstituents+process.NSubjettinesPFJets+process.hepTopTagPFJets+process.cmsTopTagPFJets+process.cmsTopTaggerPFJets+process.demo)


#write out all the trees in one file, was used for testing purposes
#process.output_step = cms.EndPath(process.output)
#process.schedule = cms.Schedule(*[process.p,process.output_step])


