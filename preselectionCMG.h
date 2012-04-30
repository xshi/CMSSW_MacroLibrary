#ifndef PRESELECTIONCMG_H
#define PRESELECTIONCMG_H

// Standard Libraries
#include <vector>
#include <string>

enum PreselType { ELE, MU, PHOT }; 

class Event;
//class TLorentzVector;
class Options;
//class TriggerInfo;
class Muon;
class Electron;
class LeptonVariables;
class ElectronVariables;
class MuonVariables;

void LeptonPreselectionCMG( const Options & opt, PreselType type );
std::vector<Muon> buildMuonCollection(const Event & ev, const LeptonVariables & leptonVars, const MuonVariables & muonVars);
std::vector<Electron> buildElectronCollection(const Event & ev, const LeptonVariables & leptonVars, const ElectronVariables & electronVars);
void selectElectronsCMG(const Event & ev, std::vector<unsigned> & electrons, double ptMin = 10);
void selectMuonsCMG(const Event & ev, std::vector<unsigned> & muons, double ptMin = 10);
void selectSoftMuonsCMG(const Event & ev, std::vector<unsigned> & softmuons, const std::vector<unsigned> & muons20);
void selectJetsCMG(const Event & ev, std::vector<unsigned> & jets, double ptMin = 30, double etaMax = 5);

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

	ElectronVariables(const Event & ev);
};

struct MuonVariables {
	unsigned m_idbits;
	unsigned m_nMatches;
	unsigned m_validMuonHits;

	MuonVariables(const Event & ev);
};
/*
void selectPhotonsCMG(const Event & ev, std::vector<unsigned> & photons);
bool triggerAcceptCMG( const Event & ev, double pt, double & weight );
bool triggerAcceptCMG( const Event & ev, PreselType type );
*/
#endif
