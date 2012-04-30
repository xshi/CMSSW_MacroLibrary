#ifndef PRESELECTION_H
#define PRESELECTION_H

// Standard Libraries
#include <vector>
#include <string>

enum PreselType { ELE, MU, PHOT }; 

class Event;
class TLorentzVector;
class Options;
class TriggerInfo;

void LeptonPreselection( const Options & opt, PreselType type );
void selectElectrons(const Event & ev, std::vector<unsigned> & electrons, double ptMin = 10);
void selectMuons(const Event & ev, std::vector<unsigned> & muons, double ptMin = 10);
void selectSoftMuons(const Event & ev, std::vector<unsigned> & softmuons, const std::vector<unsigned> & muons20);
void selectJets(const Event & ev, std::vector<unsigned> & jets, const std::vector<TLorentzVector> & leptons, double ptMin = 30, double etaMax = 5);
void selectPhotons(const Event & ev, std::vector<unsigned> & photons);
void addTriggerInfo(const Event & ev, const std::string & trigName, TriggerInfo & trig);
bool triggerAccept( const Event & ev, double pt, double & weight );
bool triggerAccept( const Event & ev, PreselType type );
//TLorentzVector calculateTrackMet(const Event & ev, const std::vector<TLorentzVector> & leptons);
void calculatePFIsol_ele(const Event & ev, 
			 int idx,
			 float & iso_nopu,
 			 float & iso_nopu_eacorr, 
 			 float & iso_nopu_dbcorr, 
 			 float & iso_nopu_pt1GeV,
			 float & iso_std,
			 float & iso_std_rhocorr); 

void calculatePFIsol_mu(const Event & ev, 
			int idx,
			float & iso_nopu,
			float & iso_nopu_eacorr, 
			float & iso_nopu_dbcorr, 
			float & iso_nopu_pt1GeV,
			float & iso_std,
			float & iso_std_rhocorr);

#endif
