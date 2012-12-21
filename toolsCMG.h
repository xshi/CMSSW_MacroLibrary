#ifndef TOOLSCMG_H
#define TOOLSCMG_H

class Event;

struct LeptonVariables {
	unsigned l1_px;
	unsigned l1_py;
	unsigned l1_pz;
	unsigned l1_en;
	unsigned l1_ptErr;
	unsigned l1_ecalIso;
	unsigned l1_hcalIso;
	unsigned l1_trkIso;
	unsigned l1_gIso;
	unsigned l1_chIso;
	unsigned l1_puchIso;
	unsigned l1_nhIso;
	unsigned l1_id;
	unsigned l1_genid;
	unsigned l1_ensf;
	unsigned l1_ensferr;
	unsigned l1_d0;
	unsigned l1_dZ;
	unsigned l1_ip3d;
	unsigned l1_trkpt;
	unsigned l1_trketa;
	unsigned l1_trkphi;
	unsigned l1_trkchi2;
	unsigned l1_trkValidPixelHits;
	unsigned l1_trkValidTrackerHits;
	unsigned l1_trkLostInnerHits;
	unsigned l1_pid;

	unsigned l2_px;
	unsigned l2_py;
	unsigned l2_pz;
	unsigned l2_en;
	unsigned l2_ptErr;
	unsigned l2_ecalIso;
	unsigned l2_hcalIso;
	unsigned l2_trkIso;
	unsigned l2_gIso;
	unsigned l2_chIso;
	unsigned l2_puchIso;
	unsigned l2_nhIso;
	unsigned l2_id;
	unsigned l2_genid;
	unsigned l2_ensf;
	unsigned l2_ensferr;
	unsigned l2_d0;
	unsigned l2_dZ;
	unsigned l2_ip3d;
	unsigned l2_trkpt;
	unsigned l2_trketa;
	unsigned l2_trkphi;
	unsigned l2_trkchi2;
	unsigned l2_trkValidPixelHits;
	unsigned l2_trkValidTrackerHits;
	unsigned l2_trkLostInnerHits;
	unsigned l2_pid;

	unsigned ln;
	unsigned ln_px;
	unsigned ln_py;
	unsigned ln_pz;
	unsigned ln_en;
	unsigned ln_ptErr;
	unsigned ln_ecalIso;
	unsigned ln_hcalIso;
	unsigned ln_trkIso;
	unsigned ln_gIso;
	unsigned ln_chIso;
	unsigned ln_puchIso;
	unsigned ln_nhIso;
	unsigned ln_id;
	unsigned ln_genid;
	unsigned ln_ensf;
	unsigned ln_ensferr;
	unsigned ln_d0;
	unsigned ln_dZ;
	unsigned ln_ip3d;
	unsigned ln_trkpt;
	unsigned ln_trketa;
	unsigned ln_trkphi;
	unsigned ln_trkchi2;
	unsigned ln_trkValidPixelHits;
	unsigned ln_trkValidTrackerHits;
	unsigned ln_trkLostInnerHits;
	unsigned ln_pid;

	LeptonVariables(const Event & ev);
};

struct ElectronVariables {
	unsigned e_idbits;
	unsigned e_hoe;
	unsigned e_dphiin;
	unsigned e_detain;
	unsigned e_sihih;
	unsigned e_sipip;
	unsigned e_r9;
	unsigned e_sce;
	unsigned e_sceta;
	unsigned e_scphi;
	unsigned e_e2x5max;
	unsigned e_e1x5;
	unsigned e_e5x5;
	unsigned e_h2te;
	unsigned e_h2tebc;
	unsigned e_ooemoop;
	unsigned e_fbrem;
	unsigned e_eopin;
	unsigned e_dEtaCalo;
	unsigned e_kfchi2;
	unsigned e_kfhits;
	unsigned e_etawidth;
	unsigned e_phiwidth;
	unsigned e_e1x5e5x5;
	unsigned e_preShowerOverRaw;
	unsigned e_eopout;

	ElectronVariables(const Event & ev);
};

struct MuonVariables {
	unsigned m_idbits;
	unsigned m_nMatches;
	unsigned m_nMatchedStations;
	unsigned m_validMuonHits;
	unsigned m_innerTrackChi2;
	unsigned m_trkLayersWithMeasurement;
	unsigned m_pixelLayersWithMeasurement;

	MuonVariables(const Event & ev);
};

struct JetVariables {
	unsigned jn;
	unsigned j_px;
	unsigned j_py;
	unsigned j_pz;
	unsigned j_en;
	unsigned j_btag;
	unsigned j_genpt;
	unsigned j_idbits;

	JetVariables(const Event & ev);
};

struct PhotonVariables {
	unsigned gn;
	unsigned g_px;
	unsigned g_py;
	unsigned g_pz;
	unsigned g_en;
	unsigned g_iso1;
	unsigned g_iso2;
	unsigned g_iso3;
	unsigned g_sihih;
	unsigned g_sipip;
	unsigned g_r9;
	unsigned g_hoe;
	unsigned g_htoe;
	unsigned g_corren;
	unsigned g_correnerr;
	unsigned g_idbits;
	
	PhotonVariables(const Event & ev);
};

struct MCVariables {
	unsigned mcparticles;
	unsigned mc_px;
	unsigned mc_py;
	unsigned mc_pz;
	unsigned mc_en;
	unsigned mc_id;

	MCVariables(const Event & ev);
};

#endif // TOOLSCMS_H
