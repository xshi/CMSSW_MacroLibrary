#ifndef PRESELECTIONCMG_H
#define PRESELECTIONCMG_H

// Standard Libraries
#include <vector>
#include <string>

enum PreselType { ELE, MU, PHOT }; 

class Event;
class TLorentzVector;
class Options;
//class TriggerInfo;
class Muon;
class Electron;
class LeptonVariables;
class ElectronVariables;
class MuonVariables;
class JetVariables;
class Jet;

void LeptonPreselectionCMG( const Options & opt, PreselType type, bool isData );
Jet smearedJet(const Jet & origJet);
TLorentzVector smearJets(std::vector<Jet> & jets);
std::vector<Muon> buildMuonCollection(const Event & ev, const LeptonVariables & leptonVars, const MuonVariables & muonVars);
std::vector<Electron> buildElectronCollection(const Event & ev, const LeptonVariables & leptonVars, const ElectronVariables & electronVars);
void selectElectronsCMG(const Event & ev, std::vector<unsigned> & electrons, double ptMin = 10);
void selectMuonsCMG(const Event & ev, std::vector<unsigned> & muons, double ptMin = 10);
void selectSoftMuonsCMG(const Event & ev, std::vector<unsigned> & softmuons, const std::vector<unsigned> & muons20);
std::vector<Jet> selectJetsCMG(const Event & ev, const JetVariables & jetVars, double ptMin = 10, double etaMax = 5);
	
/*
void selectPhotonsCMG(const Event & ev, std::vector<unsigned> & photons);
bool triggerAcceptCMG( const Event & ev, double pt, double & weight );
bool triggerAcceptCMG( const Event & ev, PreselType type );
*/
#endif
