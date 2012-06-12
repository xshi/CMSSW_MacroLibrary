#include "toolsCMG.h"
#include "event.h"

LeptonVariables::LeptonVariables(const Event & ev) {
	l1_px = ev.findVariableIndex("l1_px");
	l1_py = ev.findVariableIndex("l1_py");
	l1_pz = ev.findVariableIndex("l1_pz");
	l1_en = ev.findVariableIndex("l1_en");
	l1_ptErr = ev.findVariableIndex("l1_ptErr");
	l1_ecalIso = ev.findVariableIndex("l1_ecalIso");
	l1_hcalIso = ev.findVariableIndex("l1_hcalIso");
	l1_trkIso = ev.findVariableIndex("l1_trkIso");
	l1_gIso = ev.findVariableIndex("l1_gIso");
	l1_chIso = ev.findVariableIndex("l1_chIso");
	l1_puchIso = ev.findVariableIndex("l1_puchIso");
	l1_nhIso = ev.findVariableIndex("l1_nhIso");
	l1_id = ev.findVariableIndex("l1_id");
	l1_genid = ev.findVariableIndex("l1_genid");
	l1_ensf = ev.findVariableIndex("l1_ensf");
	l1_ensferr = ev.findVariableIndex("l1_ensferr");
	l1_d0 = ev.findVariableIndex("l1_d0");
	l1_dZ = ev.findVariableIndex("l1_dZ");
	l1_ip3d = ev.findVariableIndex("l1_ip3d");
	l1_trkpt = ev.findVariableIndex("l1_trkpt");
	l1_trketa = ev.findVariableIndex("l1_trketa");
	l1_trkphi = ev.findVariableIndex("l1_trkphi");
	l1_trkchi2 = ev.findVariableIndex("l1_trkchi2");
	l1_trkValidPixelHits = ev.findVariableIndex("l1_trkValidPixelHits");
	l1_trkValidTrackerHits = ev.findVariableIndex("l1_trkValidTrackerHits");
	l1_trkLostInnerHits = ev.findVariableIndex("l1_trkLostInnerHits");
	l1_pid = ev.findVariableIndex("l1_pid");

	l2_px = ev.findVariableIndex("l2_px");
	l2_py = ev.findVariableIndex("l2_py");
	l2_pz = ev.findVariableIndex("l2_pz");
	l2_en = ev.findVariableIndex("l2_en");
	l2_ptErr = ev.findVariableIndex("l2_ptErr");
	l2_ecalIso = ev.findVariableIndex("l2_ecalIso");
	l2_hcalIso = ev.findVariableIndex("l2_hcalIso");
	l2_trkIso = ev.findVariableIndex("l2_trkIso");
	l2_gIso = ev.findVariableIndex("l2_gIso");
	l2_chIso = ev.findVariableIndex("l2_chIso");
	l2_puchIso = ev.findVariableIndex("l2_puchIso");
	l2_nhIso = ev.findVariableIndex("l2_nhIso");
	l2_id = ev.findVariableIndex("l2_id");
	l2_genid = ev.findVariableIndex("l2_genid");
	l2_ensf = ev.findVariableIndex("l2_ensf");
	l2_ensferr = ev.findVariableIndex("l2_ensferr");
	l2_d0 = ev.findVariableIndex("l2_d0");
	l2_dZ = ev.findVariableIndex("l2_dZ");
	l2_ip3d = ev.findVariableIndex("l2_ip3d");
	l2_trkpt = ev.findVariableIndex("l2_trkpt");
	l2_trketa = ev.findVariableIndex("l2_trketa");
	l2_trkphi = ev.findVariableIndex("l2_trkphi");
	l2_trkchi2 = ev.findVariableIndex("l2_trkchi2");
	l2_trkValidPixelHits = ev.findVariableIndex("l2_trkValidPixelHits");
	l2_trkValidTrackerHits = ev.findVariableIndex("l2_trkValidTrackerHits");
	l2_trkLostInnerHits = ev.findVariableIndex("l2_trkLostInnerHits");
	l2_pid = ev.findVariableIndex("l2_pid");

	ln = ev.findVariableIndex("ln");
	ln_px = ev.findVariableIndex("ln_px");
	ln_py = ev.findVariableIndex("ln_py");
	ln_pz = ev.findVariableIndex("ln_pz");
	ln_en = ev.findVariableIndex("ln_en");
	ln_ptErr = ev.findVariableIndex("ln_ptErr");
	ln_ecalIso = ev.findVariableIndex("ln_ecalIso");
	ln_hcalIso = ev.findVariableIndex("ln_hcalIso");
	ln_trkIso = ev.findVariableIndex("ln_trkIso");
	ln_gIso = ev.findVariableIndex("ln_gIso");
	ln_chIso = ev.findVariableIndex("ln_chIso");
	ln_puchIso = ev.findVariableIndex("ln_puchIso");
	ln_nhIso = ev.findVariableIndex("ln_nhIso");
	ln_id = ev.findVariableIndex("ln_id");
	ln_genid = ev.findVariableIndex("ln_genid");
	ln_ensf = ev.findVariableIndex("ln_ensf");
	ln_ensferr = ev.findVariableIndex("ln_ensferr");
	ln_d0 = ev.findVariableIndex("ln_d0");
	ln_dZ = ev.findVariableIndex("ln_dZ");
	ln_ip3d = ev.findVariableIndex("ln_ip3d");
	ln_trkpt = ev.findVariableIndex("ln_trkpt");
	ln_trketa = ev.findVariableIndex("ln_trketa");
	ln_trkphi = ev.findVariableIndex("ln_trkphi");
	ln_trkchi2 = ev.findVariableIndex("ln_trkchi2");
	ln_trkValidPixelHits = ev.findVariableIndex("ln_trkValidPixelHits");
	ln_trkValidTrackerHits = ev.findVariableIndex("ln_trkValidTrackerHits");
	ln_trkLostInnerHits = ev.findVariableIndex("ln_trkLostInnerHits");
	ln_pid = ev.findVariableIndex("ln_pid");
}

ElectronVariables::ElectronVariables(const Event & ev) {
	e_idbits = ev.findVariableIndex("en_idbits");
	e_hoe = ev.findVariableIndex("en_hoe");
	e_dphiin = ev.findVariableIndex("en_dphiin");
	e_detain = ev.findVariableIndex("en_detain");
	e_sihih = ev.findVariableIndex("en_sihih");
	e_sipip = ev.findVariableIndex("en_sipip");
	e_r9 = ev.findVariableIndex("en_r9");
	e_sce = ev.findVariableIndex("en_sce");
	e_sceta = ev.findVariableIndex("en_sceta");
	e_scphi = ev.findVariableIndex("en_scphi");
	e_e2x5max = ev.findVariableIndex("en_e2x5max");
	e_e1x5 = ev.findVariableIndex("en_e1x5");
	e_e5x5 = ev.findVariableIndex("en_e5x5");
	e_h2te = ev.findVariableIndex("en_h2te");
	e_h2tebc = ev.findVariableIndex("en_h2tebc");
	e_ooemoop = ev.findVariableIndex("en_ooemoop");
	e_fbrem = ev.findVariableIndex("en_fbrem");
	e_eopin = ev.findVariableIndex("en_eopin");
	e_dEtaCalo = ev.findVariableIndex("en_dEtaCalo");
	e_kfchi2 = ev.findVariableIndex("en_kfchi2");
	e_kfhits = ev.findVariableIndex("en_kfhits");
	e_etawidth = ev.findVariableIndex("en_etawidth");
	e_phiwidth = ev.findVariableIndex("en_phiwidth");
	e_e1x5e5x5 = ev.findVariableIndex("en_e1x5e5x5");
	e_preShowerOverRaw = ev.findVariableIndex("en_preShowerOverRaw");
	e_eopout = ev.findVariableIndex("en_eopout");
}

MuonVariables::MuonVariables(const Event & ev) {
	m_idbits = ev.findVariableIndex("mn_idbits");
	m_nMatches = ev.findVariableIndex("mn_nMatches");
	m_validMuonHits = ev.findVariableIndex("mn_validMuonHits");
	m_innerTrackChi2 = ev.findVariableIndex("mn_innerTrackChi2");
	m_trkLayersWithMeasurement = ev.findVariableIndex("mn_trkLayersWithMeasurement");
	m_pixelLayersWithMeasurement = ev.findVariableIndex("mn_pixelLayersWithMeasurement");
}

JetVariables::JetVariables(const Event & ev) {
	jn = ev.findVariableIndex("ajn");
	j_px = ev.findVariableIndex("ajn_px");
	j_py = ev.findVariableIndex("ajn_py");
	j_pz = ev.findVariableIndex("ajn_pz");
	j_en = ev.findVariableIndex("ajn_en");
	j_btag = ev.findVariableIndex("ajn_btag3");
	j_genpt = ev.findVariableIndex("ajn_genpt");
}

MCVariables::MCVariables(const Event & ev) {
	mcparticles = ev.findVariableIndex("nmcparticles");
	mc_px = ev.findVariableIndex("mc_px");
	mc_py = ev.findVariableIndex("mc_py");
	mc_pz = ev.findVariableIndex("mc_pz");
	mc_en = ev.findVariableIndex("mc_en");
	mc_id = ev.findVariableIndex("mc_id");
}
