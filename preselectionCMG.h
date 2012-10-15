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
class Options;
class Photon;
class PhotonVariables;
class TLorentzVector;
class RooWorkspace;

void LeptonPreselectionCMG( const Options & opt, PreselType type, RooWorkspace * massPeak = nullptr );
Jet smearedJet(const Jet & origJet);
TLorentzVector smearJets(std::vector<Jet> & jets);
std::vector<Muon> buildMuonCollection(const Event & ev, const LeptonVariables & leptonVars, const MuonVariables & muonVars);
std::vector<Electron> buildElectronCollection(const Event & ev, const LeptonVariables & leptonVars, const ElectronVariables & electronVars);
void selectElectronsCMG(const Event & ev, std::vector<unsigned> & electrons, double ptMin = 10);
void selectMuonsCMG(const Event & ev, std::vector<unsigned> & muons, double ptMin = 10);
void selectSoftMuonsCMG(const Event & ev, std::vector<unsigned> & softmuons, const std::vector<unsigned> & muons20);
std::vector<Jet> selectJetsCMG(const Event & ev, const JetVariables & jetVars, double ptMin = 10, double etaMax = 5);
std::vector<Photon> selectPhotonsCMG(const Event & ev, const PhotonVariables & photonVars);
	
/*
bool triggerAcceptCMG( const Event & ev, double pt, double & weight );
bool triggerAcceptCMG( const Event & ev, PreselType type );
*/
#endif
