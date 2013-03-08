#ifndef PRESELECTIONCMG_H
#define PRESELECTIONCMG_H

// Standard Libraries
#include <vector>
#include <string>

enum PreselType { ELE, MU, EMU, PHOT }; 

class Electron;
class ElectronVariables;
class Event;
class Jet;
class JetVariables;
class LeptonVariables;
class Muon;
class MuonVariables;
class Photon;
class PhotonVariables;
class TLorentzVector;
class RooWorkspace;
class JetCorrectionUncertainty;

void LeptonPreselectionCMG( PreselType type, RooWorkspace * massPeak = nullptr );
std::vector<Muon> buildMuonCollection(const Event & ev, const LeptonVariables & leptonVars, const MuonVariables & muonVars);
std::vector<Electron> buildElectronCollection(const Event & ev, const LeptonVariables & leptonVars, const ElectronVariables & electronVars);
std::vector<Jet> selectJetsCMG(const Event & ev, const JetVariables & jetVars, JetCorrectionUncertainty  & jecUnc, TLorentzVector * diff = 0, unsigned mode = 0, double ptMin = 10, double etaMax = 4.7);
Jet smearedJet(const Jet & origJet, unsigned mode);
TLorentzVector smearJets(std::vector<Jet> & jets, unsigned mode = 0);
std::vector<Photon> selectPhotonsCMG(const Event & ev, const PhotonVariables & photonVars);
	
/*
bool triggerAcceptCMG( const Event & ev, double pt, double & weight );
bool triggerAcceptCMG( const Event & ev, PreselType type );
*/
#endif
