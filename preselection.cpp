// ROOT Libraries
#include <TFile.h>
#include <TLorentzVector.h>
// Standard Libraries
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
// Other
#include "preselection.h"
#include "event.h"
#include "options.h"
#include "triggerinfo.h"
#include "toolbox.h"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;
using std::setw;

void LeptonPreselection( const Options & opt, PreselType type ) {
	if (type == ELE)
		cout << "Entering ElectronPreselection() ..." << endl;
	else if (type == MU)
		cout << "Entering MuonPreselection() ..." << endl;

	TFile * file = new TFile( opt.checkStringOption("inpath").c_str() );
	TTree * tree = ( TTree * ) file->Get( "HZZ2l2nuAnalysis" );
	Event ev( tree );

	string outputFile( opt.checkStringOption("outpath") );
	size_t pos = outputFile.find(".root");
	if ( pos != string::npos )
		outputFile.erase(pos);

	if ( opt.checkBoolOption("addTrigger") )
		outputFile += "_triggers";

	if (type == ELE)
		outputFile += "_elePresel.root";
	else if (type == MU)
		outputFile += "_muPresel.root";
	
	cout << outputFile << endl;
	TFile * out = new TFile( outputFile.c_str(), "recreate" );

	unsigned run;
	unsigned lumi;
	unsigned event;
	double pfmet;
	int nele;
	int nmu;
	int nsoftmu;
	double l1pt;
	double l1eta;
	double l1phi;
	double l2pt;
	double l2eta;
	double l2phi;
	double zmass;
	double zpt;
	double zeta;
	double mt;
	int njet;
	int nsoftjet;
	double maxJetBTag;
	double minDeltaPhiJetMet;
	double minDeltaPhiSoftJetMet;
	int nvtx;
	int ni;
	double pftrkmet;

	// PF Isolation
	float pfisol1,pfisol1_eacorr,pfisol1_dbcorr,pfisol1_pt1GeV;
	float pfisol2,pfisol2_eacorr,pfisol2_dbcorr,pfisol2_pt1GeV;
	float stdisol1,stdisol1_rhocorr;
	float stdisol2,stdisol2_rhocorr;

	// Trigger Info
	TriggerInfo HLT_DoubleMu7;
	TriggerInfo HLT_Ele17_CaloIdL_CaloIsoVL;
	TriggerInfo HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
	TriggerInfo HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
	TriggerInfo HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30;
	TriggerInfo HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30;
	TriggerInfo HLT_Ele32_CaloIdL_CaloIsoVL_SC17;
	TriggerInfo HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
	TriggerInfo HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT;
	TriggerInfo HLT_IsoMu17;
	TriggerInfo HLT_IsoMu24;
	TriggerInfo HLT_IsoMu30;
	TriggerInfo HLT_IsoMu24_eta2p1;
	TriggerInfo HLT_IsoMu30_eta2p1;
	TriggerInfo HLT_Mu13_Mu8;
	TriggerInfo HLT_Mu17_Mu8;
	TriggerInfo HLT_Mu17_TkMu8;

	TTree * smallTree = new TTree("HZZ2l2nuAnalysis", "HZZ2l2nu Analysis Tree");
	smallTree->Branch( "Run", &run, "Run/i" );
	smallTree->Branch( "Lumi", &lumi, "Lumi/i" );
	smallTree->Branch( "Event", &event, "Event/i" );
	smallTree->Branch( "PFMET", &pfmet, "PFMET/D" );
	smallTree->Branch( "NELE", &nele, "NELE/I" );
	smallTree->Branch( "NMU", &nmu, "NMU/I" );
	smallTree->Branch( "NSOFTMU", &nsoftmu, "NSOFTMU/I" );
	smallTree->Branch( "L1PT", &l1pt, "L1PT/D" );
	smallTree->Branch( "L1ETA", &l1eta, "L1ETA/D" );
	smallTree->Branch( "L1PHI", &l1phi, "L1PHI/D" );
	smallTree->Branch( "L2PT", &l2pt, "L2PT/D" );
	smallTree->Branch( "L2ETA", &l2eta, "L2ETA/D" );
	smallTree->Branch( "L2PHI", &l2phi, "L2PHI/D" );
	smallTree->Branch( "ZMASS", &zmass, "ZMASS/D" );
	smallTree->Branch( "ZPT", &zpt, "ZPT/D" );
	smallTree->Branch( "ZETA", &zeta, "ZETA/D" );
	smallTree->Branch( "MT", &mt, "MT/D" );
	smallTree->Branch( "NJET", &njet, "NJET/I" );
	smallTree->Branch( "NSOFTJET", &nsoftjet, "NSOFTJET/I" );
	smallTree->Branch( "MAXJETBTAG", &maxJetBTag, "MAXJETBTAG/D" );
	smallTree->Branch( "MINDPJETMET", &minDeltaPhiJetMet, "MINDPJETMET/D" );
	smallTree->Branch( "MINDPSOFTJETMET", &minDeltaPhiSoftJetMet, "MINDPSOFTJETMET/D" );
	smallTree->Branch( "NVTX", &nvtx, "NVTX/I" );
	smallTree->Branch( "nInter" , &ni, "nInter/I" );
	smallTree->Branch( "PFTRKMET", &pftrkmet, "PFTRKMET/D" );

	// Isolation
	if ( opt.checkBoolOption("addIsolation") ) {
		smallTree->Branch( "PFISOL1", &pfisol1, "PFISOL1/F" );
		smallTree->Branch( "PFISOL1_EACORR", &pfisol1_eacorr, "PFISOL1_EACORR/F" );
		smallTree->Branch( "PFISOL1_DBCORR", &pfisol1_dbcorr, "PFISOL1_DBCORR/F" );
		smallTree->Branch( "PFISOL1_PT1GEV", &pfisol1_pt1GeV, "PFISOL1_PF1GEV/F" );
		smallTree->Branch( "STDISOL1", &stdisol1, "STDISOL1/F" );
		smallTree->Branch( "STDISOL1_RHOCORR", &stdisol1_rhocorr, "STDISOL1_RHOCORR/F" );

		smallTree->Branch( "PFISOL2", &pfisol2, "PFISOL2/F" );
		smallTree->Branch( "PFISOL2_EACORR", &pfisol2_eacorr, "PFISOL2_EACORR/F" );
		smallTree->Branch( "PFISOL2_DBCORR", &pfisol2_dbcorr, "PFISOL2_DBCORR/F" );
		smallTree->Branch( "PFISOL2_PT1GEV", &pfisol2_pt1GeV, "PFISOL2_PF1GEV/F" );
		smallTree->Branch( "STDISOL2", &stdisol2, "STDISOL2/F" );
		smallTree->Branch( "STDISOL2_RHOCORR", &stdisol2_rhocorr, "STDISOL2_RHOCORR/F" );
	}
	if ( opt.checkBoolOption("addTrigger") ) {
		smallTree->Branch( "HLT_DoubleMu7", &HLT_DoubleMu7 );
		smallTree->Branch( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", &HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL );
		smallTree->Branch( "HLT_Ele17_CaloIdL_CaloIsoVL", &HLT_Ele17_CaloIdL_CaloIsoVL );
		smallTree->Branch( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL );
		smallTree->Branch( "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30", &HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30 );
		smallTree->Branch( "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30", &HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 );
		smallTree->Branch( "HLT_Ele32_CaloIdL_CaloIsoVL_SC17", &HLT_Ele32_CaloIdL_CaloIsoVL_SC17 );
		smallTree->Branch( "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17", &HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17 );
		smallTree->Branch( "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT", &HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT );
		smallTree->Branch( "HLT_IsoMu17", &HLT_IsoMu17 );
		smallTree->Branch( "HLT_IsoMu24", &HLT_IsoMu24 );
		smallTree->Branch( "HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1 );
		smallTree->Branch( "HLT_IsoMu30_eta2p1", &HLT_IsoMu30_eta2p1 );
		smallTree->Branch( "HLT_IsoMu30", &HLT_IsoMu30 );
		smallTree->Branch( "HLT_Mu13_Mu8", &HLT_Mu13_Mu8 );
		smallTree->Branch( "HLT_Mu17_Mu8", &HLT_Mu17_Mu8 );
		smallTree->Branch( "HLT_Mu17_TkMu8", &HLT_Mu17_TkMu8 );
	}

	bool isData = opt.checkBoolOption("isData");

	unsigned long nentries = tree->GetEntries();
	for ( unsigned long i = 0; i < nentries; i++ ) {
		tree->GetEntry( i );

		run = ev.getSingleVariableValue<unsigned>("Run");
		lumi = ev.getSingleVariableValue<unsigned>("LumiSection");
		event = ev.getSingleVariableValue<unsigned>("Event");

		if (isData) {
			if (!triggerAccept(ev, type))
				continue;
		}

		pfmet = ev.getSingleVariableValue<float>("PFMET");

		vector<unsigned> leptons20;
		if (type == ELE)
			selectElectrons( ev, leptons20, 20 );
		else if (type == MU)
			selectMuons( ev, leptons20, 20 );

		if ( leptons20.size() != 2 )
			continue;

		vector<unsigned> electrons10;
		selectElectrons( ev, electrons10 );
		nele = electrons10.size();

		vector<unsigned> muons10;
		selectMuons( ev, muons10 );
		nmu = muons10.size();
		
		if (nmu + nele > 2)
			continue;

		vector<unsigned> softMuons;
		if (type == ELE) {
			vector<unsigned> empty;
			empty.clear();
			selectSoftMuons( ev, softMuons, empty );
		} else if (type == MU)
			selectSoftMuons( ev, softMuons, leptons20 );

		nsoftmu = softMuons.size();
		if (nsoftmu > 0)
			continue;

		string leptonsType;
		if (type == ELE)
			leptonsType = "Electrons";
		else if (type == MU)
			leptonsType = "Muons";

		const vector<int> & lepCharge = *ev.getVectorIntAdr( (leptonsType + "_CHARGE").c_str() );
		if ( lepCharge[ leptons20[0] ] == lepCharge[ leptons20[1] ] )
			continue;

		const vector<float> & lepE = *ev.getVectorFloatAdr( (leptonsType + "_E").c_str() );
		const vector<float> & lepPt = *ev.getVectorFloatAdr( (leptonsType + "_PT").c_str() );
		const vector<float> & lepEta = *ev.getVectorFloatAdr( (leptonsType + "_ETA").c_str() );
		const vector<float> & lepPhi = *ev.getVectorFloatAdr( (leptonsType + "_PHI").c_str() );
		TLorentzVector lep1;
		lep1.SetPtEtaPhiE( lepPt[ leptons20[0] ], lepEta[ leptons20[0] ], lepPhi[ leptons20[0] ], lepE[ leptons20[0] ]);
		TLorentzVector lep2;
		lep2.SetPtEtaPhiE( lepPt[ leptons20[1] ], lepEta[ leptons20[1] ], lepPhi[ leptons20[1] ], lepE[ leptons20[1] ]);

		if (lep2.Pt() > lep1.Pt()) {
			TLorentzVector temp = lep1;
			lep1 = lep2;
			lep2 = temp;
		}

		l1pt = lep1.Pt();
		l1eta = lep1.Eta();
		l1phi = lep1.Phi();

		l2pt = lep2.Pt();
		l2eta = lep2.Eta();
		l2phi = lep2.Phi();

		TLorentzVector Zcand = lep1 + lep2;
		zpt = Zcand.Pt();
		zeta = Zcand.Eta();
		zmass = Zcand.M();
		if (zmass < 76.1876 || zmass > 106.1876 || zpt < 55)
			continue;

		double metphi = ev.getSingleVariableValue<float>("PFMET_PHI");
		double metPx = pfmet * cos( metphi );
		double metPy = pfmet * sin( metphi );
		double px = metPx + Zcand.Px();
		double py = metPy + Zcand.Py();
		double pt2 = px * px + py * py;
		double e = sqrt(zpt * zpt + zmass * zmass) + sqrt(pfmet * pfmet + zmass * zmass);
		double mt2 = e * e - pt2;
		mt = (mt2 > 0) ? sqrt(mt2) : 0;

		vector<TLorentzVector> leptons;
		leptons.push_back( lep1 );
		leptons.push_back( lep2 );

		vector<unsigned> jets;
		selectJets( ev, jets, leptons );
		njet = jets.size();

		vector<unsigned> softjets;
		selectJets( ev, softjets, leptons, 15 );
		nsoftjet = softjets.size();

		const vector<float> & btag = *ev.getVectorFloatAdr("Jets_BTAG_TrackCountingHighEff");
		const vector<float> & jetEta = *ev.getVectorFloatAdr("Jets_ETA");
		maxJetBTag = -999;
		for ( int j = 0; j < njet; ++j ) {
			if ( btag[ jets[j] ] > maxJetBTag && fabs(jetEta[ jets[j] ]) < 2.4 )
				maxJetBTag = btag[ jets[j] ];
		}
		/*
		unsigned temprun = 160956;
		unsigned templumi = 53;
		unsigned tempevent = 27115287;
		if (run == temprun && lumi == templumi  && event == tempevent) {
			cout << "Found!" << endl;
			cout << "Trigger: " << triggerAccept(ev, type) << endl;
			cout << "Nlep20: " << leptons20.size() << endl;
			cout << "Nlep: " << nmu + nele << endl;
			cout << "Nsoftmu: " << nsoftmu << endl;
			cout << "Same charge: " << (lepCharge[ leptons20[0] ] == lepCharge[ leptons20[1] ]) << endl;
			cout << "zmass: " << zmass << endl;
			cout << "zpt: " << zpt << endl;
			cout << "maxJetBTag: " << maxJetBTag << endl;
			cin.get();
		}
		*/
		if (maxJetBTag > 2.0)
			continue;

		const vector<float> & delPhiJetMet = *ev.getVectorFloatAdr("Jets_DELTAPHIMET");
		minDeltaPhiJetMet = 999;
		for ( int j = 0; j < njet; ++j ) {
			double tempDelPhiJetMet = fabs(delPhiJetMet[ jets[j] ]);
			if ( tempDelPhiJetMet < minDeltaPhiJetMet )
				minDeltaPhiJetMet = tempDelPhiJetMet;
		}
		minDeltaPhiSoftJetMet = 999;
		for ( int j = 0; j < nsoftjet; ++j ) {
			double tempDelPhiSoftJetMet = fabs(delPhiJetMet[ softjets[j] ]);
			if ( tempDelPhiSoftJetMet < minDeltaPhiSoftJetMet )
				minDeltaPhiSoftJetMet = tempDelPhiSoftJetMet;
		}

		nvtx = ev.getSingleVariableValue<int>("NVTX");

		ni = nvtx;
		const vector<int> & bx = *ev.getVectorIntAdr("PUP_BunchCrossing");
		const vector<int> & ninter = *ev.getVectorIntAdr("PUP_NumInteractions");
		for ( unsigned j = 0; j < bx.size(); ++j )
			if ( bx[j] == 0 )
				ni = ninter[j];
		
		pftrkmet = ev.getSingleVariableValue<float>("PFTRACKMET");

		if ( opt.checkBoolOption("addIsolation") ) {
			float iso_nopu = -1.;
			float iso_nopu_eacorr = -1.;
			float iso_nopu_dbcorr = -1.;
			float iso_nopu_pt1GeV = -1.;
			float iso_std = -1.;
			float iso_std_rhocorr = -1.;

			if (type == ELE)
				calculatePFIsol_ele(ev, leptons20[0], iso_nopu, iso_nopu_eacorr, iso_nopu_dbcorr, iso_nopu_pt1GeV, iso_std, iso_std_rhocorr);
			else if (type == MU)
				calculatePFIsol_mu(ev, leptons20[0], iso_nopu, iso_nopu_eacorr, iso_nopu_dbcorr, iso_nopu_pt1GeV, iso_std, iso_std_rhocorr);

			pfisol1 = iso_nopu / l1pt;
			pfisol1_eacorr = iso_nopu_eacorr / l1pt;
			pfisol1_dbcorr = iso_nopu_dbcorr / l1pt;
			pfisol1_pt1GeV = iso_nopu_pt1GeV / l1pt;
			stdisol1 = iso_std / l1pt;
			stdisol1_rhocorr = iso_std_rhocorr / l1pt;

			if (type == ELE)
				calculatePFIsol_ele(ev, leptons20[1], iso_nopu, iso_nopu_eacorr, iso_nopu_dbcorr, iso_nopu_pt1GeV, iso_std, iso_std_rhocorr);
			else if (type == MU)
				calculatePFIsol_mu(ev, leptons20[1], iso_nopu, iso_nopu_eacorr, iso_nopu_dbcorr, iso_nopu_pt1GeV, iso_std, iso_std_rhocorr);

			pfisol2 = iso_nopu / l2pt;
			pfisol2_eacorr = iso_nopu_eacorr / l2pt;
			pfisol2_dbcorr = iso_nopu_dbcorr / l2pt;
			pfisol2_pt1GeV = iso_nopu_pt1GeV / l2pt;
			stdisol2 = iso_std / l2pt;
			stdisol2_rhocorr = iso_std_rhocorr / l2pt;
		}
		
		if ( opt.checkBoolOption("addTrigger") ) {
			addTriggerInfo(ev, "HLT_DoubleMu7", HLT_DoubleMu7);
			addTriggerInfo(ev, "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
			addTriggerInfo(ev, "HLT_Ele17_CaloIdL_CaloIsoVL", HLT_Ele17_CaloIdL_CaloIsoVL);
			addTriggerInfo(ev, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
			addTriggerInfo(ev, "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30", HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30);
			addTriggerInfo(ev, "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30", HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
			addTriggerInfo(ev, "HLT_Ele32_CaloIdL_CaloIsoVL_SC17", HLT_Ele32_CaloIdL_CaloIsoVL_SC17);
			addTriggerInfo(ev, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17", HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17);
			addTriggerInfo(ev, "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT", HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
			addTriggerInfo(ev, "HLT_IsoMu17", HLT_IsoMu17);
			addTriggerInfo(ev, "HLT_IsoMu24", HLT_IsoMu24);
			addTriggerInfo(ev, "HLT_IsoMu24_eta2p1", HLT_IsoMu24_eta2p1);
			addTriggerInfo(ev, "HLT_IsoMu30_eta2p1", HLT_IsoMu30_eta2p1);
			addTriggerInfo(ev, "HLT_IsoMu30", HLT_IsoMu30);
			addTriggerInfo(ev, "HLT_Mu13_Mu8", HLT_Mu13_Mu8);
			addTriggerInfo(ev, "HLT_Mu17_Mu8", HLT_Mu17_Mu8);
			addTriggerInfo(ev, "HLT_Mu17_TkMu8", HLT_Mu17_TkMu8);
		}

		smallTree->Fill();
	}
	delete file;
	smallTree->Write("", TObject::kOverwrite);
	delete smallTree;
	delete out;
}

void selectElectrons(const Event & ev, vector<unsigned> & electrons, double ptMin) {
	const vector<float> & pt = *ev.getVectorFloatAdr("Electrons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Electrons_ETA");
	const vector<float> & clEta = *ev.getVectorFloatAdr("Electrons_superClusterETA");
	const vector<int> & eleID = *ev.getVectorIntAdr("Electrons_simpleEleId80cIso");
	const vector<float> & dxy = *ev.getVectorFloatAdr("Electrons_dxy");
	const vector<float> & dz = *ev.getVectorFloatAdr("Electrons_dz");
	const vector<float> & sigIeIe = *ev.getVectorFloatAdr("Electrons_sigmaIetaIeta");
	const vector<float> & delPhi = *ev.getVectorFloatAdr("Electrons_deltaPhiSuperClusterTrackAtVtx");
	const vector<float> & delEta = *ev.getVectorFloatAdr("Electrons_deltaEtaSuperClusterTrackAtVtx");
	const vector<float> & hOe = *ev.getVectorFloatAdr("Electrons_hadronicOverEm");
	const vector<float> & trackIso = *ev.getVectorFloatAdr("Electrons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Electrons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Electrons_HCALISO");
	double rho = ev.getSingleVariableValue<float>("RhoCorrectionFactor");

	electrons.clear();
	const int nEle = pt.size();
	for ( int i = 0; i < nEle; ++i ) {
		if ( pt[i] > ptMin  &&
				fabs(eta[i]) < 2.5 &&
				!( fabs(clEta[i]) > 1.4442 && fabs(clEta[i]) < 1.566) &&
				eleID[i] > 3.5 &&
				fabs(dxy[i]) < 0.02 &&
				fabs(dz[i]) < 0.1
		   ) {
			if ( fabs(clEta[i]) < 1.4442 ) {
				double isol = ( trackIso[i] + max(ecalIso[i] - 1, 0) + hcalIso[i] - (rho * M_PI * 0.3 * 0.3) ) / pt[i];
				if ( sigIeIe[i] < 0.01 &&
						fabs(delPhi[i]) < 0.06 &&
						fabs(delEta[i]) < 0.004 &&
						hOe[i] < 0.04 &&
						isol < 0.1
				   ) {
					electrons.push_back( i );
				}
			} else if ( fabs(clEta[i]) > 1.566 ) {
				double isol = ( trackIso[i] + ecalIso[i] + hcalIso[i] - (rho * M_PI * 0.3 * 0.3) ) / pt[i];
				if ( sigIeIe[i] < 0.03 &&
						fabs(delPhi[i]) < 0.03 &&
						fabs(delEta[i]) < 0.007 &&
						hOe[i] < 0.1 &&
						isol < 0.1
				   ) {
					electrons.push_back( i );
				}
			}
		}
	}
}

void selectMuons(const Event & ev, vector<unsigned> & muons, double ptMin) {
	const vector<int> & global = *ev.getVectorIntAdr("Muons_GLOBAL");
	const vector<int> & tracker = *ev.getVectorIntAdr("Muons_TRACKER");
	const vector<float> & pt = *ev.getVectorFloatAdr("Muons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Muons_ETA");
	const vector<float> & ptErr = *ev.getVectorFloatAdr("Muons_innerTrackPtError");
	const vector<float> & chi2 = *ev.getVectorFloatAdr("Muons_normChi2");
	const vector<float> & dxy = *ev.getVectorFloatAdr("Muons_dxy");
	const vector<float> & dz = *ev.getVectorFloatAdr("Muons_dz");
	const vector<int> & nTrackHit = *ev.getVectorIntAdr("Muons_TrackerHits");
	const vector<int> & nPixelHit = *ev.getVectorIntAdr("Muons_PixelHits");
	const vector<int> & nMuonHit = *ev.getVectorIntAdr("Muons_ValidMuonHits");
	const vector<int> & nStation = *ev.getVectorIntAdr("Muons_numberOfMatchedStations");
	const vector<float> & trackIso = *ev.getVectorFloatAdr("Muons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Muons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Muons_HCALISO");
	double rho = ev.getSingleVariableValue<float>("RhoCorrectionFactor");

	muons.clear();
	const int nMu = pt.size();
	for ( int i = 0; i < nMu; ++i ) {
		if ( global[i] == 1 &&
				tracker[i] == 1 &&
				pt[i] > ptMin &&
				fabs(eta[i]) < 2.4 &&
				ptErr[i] / pt[i] < 0.1 &&
				chi2[i] < 10 &&
				fabs(dxy[i]) < 0.02 &&
				fabs(dz[i]) < 0.1 &&
				nTrackHit[i] > 10 &&
				nPixelHit[i] >= 1 &&
				nMuonHit[i] >= 1 &&
				nStation[i] >= 2
		   ) {
			double isol = ( trackIso[i] + ecalIso[i] + hcalIso[i] - (rho * M_PI * 0.3 * 0.3) ) / pt[i];
			if ( isol < 0.15 )
				muons.push_back( i );
			//cout << isol << endl;
		}
	}
}

void selectSoftMuons(const Event & ev, vector<unsigned> & softmuons, const vector<unsigned> & muons20) {
	const vector<int> & tracker = *ev.getVectorIntAdr("Muons_TRACKER");
	const vector<float> & pt = *ev.getVectorFloatAdr("Muons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Muons_ETA");
	const vector<float> & dxy = *ev.getVectorFloatAdr("Muons_dxy");
	const vector<float> & dz = *ev.getVectorFloatAdr("Muons_dz");
	const vector<int> & nTrackHit = *ev.getVectorIntAdr("Muons_TrackerHits");
	const vector<int> & muID = *ev.getVectorIntAdr("Muons_TMLastStationAngTight");

	const vector<float> & trackIso = *ev.getVectorFloatAdr("Muons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Muons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Muons_HCALISO");
	double rho = ev.getSingleVariableValue<float>("RhoCorrectionFactor");

	softmuons.clear();
	const unsigned nMu = pt.size();
	for ( unsigned i = 0; i < nMu; ++i ) {
		for ( unsigned kk = 0; kk < muons20.size(); ++kk ) {
			if ( muons20[kk] == i )
				continue;
		}
		if ( tracker[i] == 1 &&
				pt[i] > 3 &&
				fabs(eta[i]) < 2.4 &&
				fabs(dxy[i]) < 0.2 &&
				fabs(dz[i]) < 0.2 &&
				nTrackHit[i] > 10 &&
				muID[i] == 1
		   ) {
			if ( pt[i] > 20 ) {
				double isol = ( trackIso[i] + ecalIso[i] + hcalIso[i] - (rho * M_PI * 0.3 * 0.3) ) / pt[i];
				if ( isol > 0.1 )
					softmuons.push_back( i );
			} else
				softmuons.push_back( i );
		}
	}
}

void selectJets(const Event & ev, vector<unsigned> & jets, const vector<TLorentzVector> & leptons, double ptMin, double etaMax) {
	const vector<float> & pt = *ev.getVectorFloatAdr("Jets_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Jets_ETA");
	const vector<float> & phi = *ev.getVectorFloatAdr("Jets_PHI");
	const vector<float> & nHEF = *ev.getVectorFloatAdr("Jets_neutralHadronEnergyFraction");
	const vector<float> & nEEF = *ev.getVectorFloatAdr("Jets_neutralEmEnergyFraction");
	const vector<int> & nC = *ev.getVectorIntAdr("Jets_nConstituents");
	const vector<float> & cHEF = *ev.getVectorFloatAdr("Jets_chargedHadronEnergyFraction");
	const vector<int> & cM = *ev.getVectorIntAdr("Jets_chargedMultiplicity");
	const vector<float> & cEEF = *ev.getVectorFloatAdr("Jets_chargedEmEnergyFraction");


	jets.clear();
	const int nJets = pt.size();
	for ( int i = 0; i < nJets; ++i ) {
		bool passedID = false;
		if ( pt[i] > ptMin &&
				fabs(eta[i]) < etaMax
				&& nHEF[i] < 0.99 &&
				nEEF[i] < 0.99 &&
				nC[i] > 1
		   ) {
			if ( fabs(eta[i]) < 2.4 ) {
				if ( cHEF[i] > 0.0 &&
						cM[i] > 0 &&
						cEEF[i] < 0.99
				   )
					passedID = true;
			} else
				passedID = true;
		}
		bool nearLepton = false;
		for (unsigned j = 0; j < leptons.size(); ++j) {
			if ( deltaR( eta[i], phi[i], leptons[j].Eta(), leptons[j].Phi() ) < 0.5 )
				nearLepton = true;
		}
		if ( passedID && !nearLepton )
			jets.push_back( i );
	}
}

void selectPhotons(const Event & ev, vector<unsigned> & photons) {
	/*
	const vector<float> & pt = *ev.getVectorFloatAdr("Photons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Photons_superClusterETA");
	const vector<int> & pixSeed = *ev.getVectorIntAdr("Photons_hasPixelSeed");
	const vector<float> & sigIeIe = *ev.getVectorFloatAdr("Photons_sigmaIetaIeta");
	const vector<float> & hOe = *ev.getVectorFloatAdr("Photons_hadronicOverEm");
	const vector<float> & trkIso = *ev.getVectorFloatAdr("Photons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Photons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Photons_HCALISO");

	photons.clear();
	const unsigned nPhot = pt.size();
	for ( unsigned i = 0; i < nPhot; ++i ) {
		if ( pt[i] > 55 &&
				fabs(eta[i]) < 1.4442 &&
				pixSeed[i] == 0 &&
				sigIeIe[i] > 0.001 &&
				sigIeIe[i] < 0.013 &&
				hOe[i] < 0.05 &&
				trkIso[i] < (2.0 + 0.001 * pt[i]) &&
				ecalIso[i] < (4.2 + 0.006 * pt[i]) &&
				hcalIso[i] < (2.2 + 0.0025 * pt[i])
		) {
			photons.push_back( i );
		}
	}
	*/
	//*
	const vector<float> & pt = *ev.getVectorFloatAdr("Photons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Photons_superClusterETA");
	const vector<float> & r9 = *ev.getVectorFloatAdr("Photons_R9");
	const vector<float> & IsoSumOverEt = *ev.getVectorFloatAdr("Photons_IsoSumOverEt");
	const vector<float> & IsoSumOverEtWorst = *ev.getVectorFloatAdr("Photons_IsoSumOverEtWorst");
	const vector<float> & TrkIsoOverEt = *ev.getVectorFloatAdr("Photons_TrkIsoOverEt");
	const vector<float> & sigIeIe = *ev.getVectorFloatAdr("Photons_sigmaIetaIeta");
	const vector<float> & hOe = *ev.getVectorFloatAdr("Photons_hadronicOverEm");
	const vector<float> & delR = *ev.getVectorFloatAdr("Photons_dREleTrack");
	const vector<float> & nLostHits = *ev.getVectorFloatAdr("Photons_nLostHitsEleTrack");

	photons.clear();
	const unsigned nPhot = pt.size();
	for ( unsigned i = 0; i < nPhot; ++i ) {
		bool cat1 = (fabs(eta[i]) < 1.4442) && (r9[i] > 0.94);
		bool cat2 = (fabs(eta[i]) < 1.4442) && (r9[i] <= 0.94);
		bool cat3 = (fabs(eta[i]) > 1.566 && fabs(eta[i]) < 2.5) && (r9[i] > 0.94);
		bool cat4 = (fabs(eta[i]) > 1.566 && fabs(eta[i]) < 2.5) && (r9[i] <= 0.94);
		
		bool passID = false;
		if ( pt[i] > 55 ) {
			if ( cat1 ) {
				if ( IsoSumOverEt[i] < 3.8 &&
						IsoSumOverEtWorst[i] < 11.7 &&
						TrkIsoOverEt[i] < 3.5 &&
						sigIeIe[i] < 0.0106 &&
						hOe[i] < 0.082 &&
						r9[i] > 0.94 )
					passID = true;
			}
			if ( cat2 ) {
				if ( IsoSumOverEt[i] < 2.2 &&
						IsoSumOverEtWorst[i] < 3.4 &&
						TrkIsoOverEt[i] < 2.2 &&
						sigIeIe[i] < 0.0097 &&
						hOe[i] < 0.062 &&
						r9[i] > 0.36
						// && delR[i] < 0.062
						)
					passID = true;
			}
			if ( cat3 ) {
				if ( IsoSumOverEt[i] < 1.77 &&
						IsoSumOverEtWorst[i] < 3.9 &&
						TrkIsoOverEt[i] < 2.3 &&
						sigIeIe[i] < 0.028 &&
						hOe[i] < 0.065 &&
						r9[i] > 0.94 )
					passID = true;
			}
			if ( cat4 ) {
				if ( IsoSumOverEt[i] < 1.29 &&
						IsoSumOverEtWorst[i] < 1.84 &&
						TrkIsoOverEt[i] < 1.45 &&
						sigIeIe[i] < 0.027 &&
						hOe[i] < 0.048 &&
						r9[i] > 0.32 )
					passID = true;
			}
		}
		bool passEleVeto = true;
		if ( delR[i] < 0.5 && nLostHits[i] == 0 )
			passEleVeto = false;
		if (passID && passEleVeto)
			photons.push_back( i );
	}
	//*/
}

void addTriggerInfo(const Event & ev, const std::string & trigName, TriggerInfo & trig) {
	trig.name = trigName;
	trig.version = ev.getSingleVariableValue<int>( trigName + "_version" );
	trig.accept = ev.getSingleVariableValue<int>( trigName + "_accept" );
	trig.prescale = ev.getSingleVariableValue<int>( trigName + "_prescale" );
	trig.TO_PT = *ev.getVectorFloatAdr( trigName + "_TO_PT" );
	trig.TO_ETA = *ev.getVectorFloatAdr( trigName + "_TO_ETA" );
	trig.TO_PHI = *ev.getVectorFloatAdr( trigName + "_TO_PHI" );
}

bool triggerAccept( const Event & ev, double pt, double & weight ) {
	bool hlt50 = (ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_accept") == 1);
	bool hlt50i = (ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_IsoL_accept") == 1);
	bool hlt75 = (ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_accept") == 1);
	bool hlt75i = (ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_IsoL_accept") == 1);
	bool hlt90 = (ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_accept") == 1);
	bool hlt90i = (ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_IsoL_accept") == 1);
	int prescale50 = ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_prescale");
	int prescale50i = ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_IsoL_prescale");
	int prescale75 = ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_prescale");
	int prescale75i = ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_IsoL_prescale");
	int prescale90 = ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_prescale");
	int prescale90i = ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_IsoL_prescale");

	if ( pt >= 55 && pt < 80 ) {
		if ( hlt50 || hlt50i ) {
			vector<int> prescales;
			if (hlt50)
				prescales.push_back( prescale50 );
			if (hlt50i)
				prescales.push_back( prescale50i );
			weight = min<int>( prescales );
			return true;
		}
	} else if ( pt >= 80 && pt < 95 ) {
		if ( hlt50 || hlt50i || hlt75 || hlt75i ) {
			vector<int> prescales;
			if (hlt50)
				prescales.push_back( prescale50 );
			if (hlt50i)
				prescales.push_back( prescale50i );
			if (hlt75)
				prescales.push_back( prescale75 );
			if (hlt75i)
				prescales.push_back( prescale75i );
			weight = min<int>( prescales );
			return true;
		}
	} else if ( pt >= 95 ) {
		if ( hlt50 || hlt50i || hlt75 || hlt75i || hlt90 || hlt90i ) {
			vector<int> prescales;
			if (hlt50)
				prescales.push_back( prescale50 );
			if (hlt50i)
				prescales.push_back( prescale50i );
			if (hlt75)
				prescales.push_back( prescale75 );
			if (hlt75i)
				prescales.push_back( prescale75i );
			if (hlt90)
				prescales.push_back( prescale90 );
			if (hlt90i)
				prescales.push_back( prescale90i );
			weight = min<int>( prescales );
			return true;
		}
	} else {
		cout << "ERROR: Unknown preselection type!" << endl;
		exit( EXIT_FAILURE );
	}
	return false;
}

bool triggerAccept( const Event & ev, PreselType type ) {
	int run = ev.getSingleVariableValue<unsigned>("Run");
	if ( type == ELE ) {
		if ( run >= 160329 && run <= 170064 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltEle17cEle8c = (ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_accept") == 1);
			if ( hltEle17cEle8c )
				return true;
		} else if ( run >= 170065 && run <= 180296 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltEle17ctEle8ct = (ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_accept") == 1);
			if ( hltEle17ctEle8ct )
				return true;
		} else {
			cout << "ERROR: Unknown run: " << run << "!" << endl;
			exit( EXIT_FAILURE );
		}
	} else if ( type == MU ) {
		if ( run >= 160329 && run <= 163869 ) {
			if ( ev.getSingleVariableValue<int>("HLT_DoubleMu7_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltDoubleMu7 = (ev.getSingleVariableValue<int>("HLT_DoubleMu7_accept") == 1);
			if ( hltDoubleMu7 )
				return true;
		} else if ( run >= 165071 && run <= 178410 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Mu13_Mu8_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltMu13Mu8 = (ev.getSingleVariableValue<int>("HLT_Mu13_Mu8_accept") == 1);
			if ( hltMu13Mu8 )
				return true;
		} else if ( run >= 178411 && run <= 180296 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Mu17_Mu8_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltMu17Mu8 = (ev.getSingleVariableValue<int>("HLT_Mu17_Mu8_accept") == 1);
			if ( hltMu17Mu8 )
				return true;
		} else {
			cout << "ERROR: Unknown run: " << run << "!" << endl;
			exit( EXIT_FAILURE );
		}
	} else {
		cout << "ERROR: Unknown preselection type!" << endl;
		exit( EXIT_FAILURE );
	}
	return false;
}

/*
TLorentzVector calculateTrackMet(const Event & ev, const std::vector<TLorentzVector> & leptons) {
	const vector<float> * PFCen = ev.getVectorFloatAdr("PFCandidate_E");
	const vector<float> * PFCpt = ev.getVectorFloatAdr("PFCandidate_PT");
	const vector<float> * PFCeta = ev.getVectorFloatAdr("PFCandidate_ETA");
	const vector<float> * PFCphi = ev.getVectorFloatAdr("PFCandidate_PHI");
	const vector<float> * PFCtrackEta = ev.getVectorFloatAdr("PFCandidate_Track_ETA");
	const vector<float> * PFCtrackPhi = ev.getVectorFloatAdr("PFCandidate_Track_PHI");

	double trkMet_x = 0;
	double trkMet_y = 0;

	int nPFC = PFCen->size();
	for ( int i = 0; i < nPFC; ++i ) {
		double delR1 = deltaR( leptons[0].Eta(), leptons[0].Phi(), (*PFCtrackEta)[i], (*PFCtrackPhi)[i]);
		double delR2 = deltaR( leptons[1].Eta(), leptons[1].Phi(), (*PFCtrackEta)[i], (*PFCtrackPhi)[i]);
		if ( delR1 > 0.1 && delR2 > 0.1 ) {
			TLorentzVector pfCand;
			pfCand.SetPtEtaPhiE( (*PFCpt)[i], (*PFCeta)[i], (*PFCphi)[i], (*PFCen)[i]);
			trkMet_x += pfCand.Px();
			trkMet_y += pfCand.Py();
		}
	}

	trkMet_x = trkMet_x + leptons[0].Px() + leptons[1].Px();
	trkMet_y = trkMet_y + leptons[0].Py() + leptons[1].Py();
	
	trkMet_x = - trkMet_x;
	trkMet_y = - trkMet_y;

	double trkMet = sqrt( trkMet_x * trkMet_x + trkMet_y * trkMet_y );
	
	TLorentzVector trackMet;
	trackMet.SetPxPyPzE( trkMet_x, trkMet_y, 0, trkMet );

	return trackMet;
}
*/

void calculatePFIsol_ele(const Event & ev, 
		int idx,
		float & iso_nopu,
		float & iso_nopu_eacorr,
		float & iso_nopu_dbcorr,
		float & iso_nopu_pt1GeV,
		float & iso_std,
		float & iso_std_rhocorr) {

	float rho = ev.getSingleVariableValue<float>("RhoCorrectionFactor");
	const vector<float> & eleEtaVec = *ev.getVectorFloatAdr("Electrons_ETA");
	const vector<float> & elePhiVec = *ev.getVectorFloatAdr("Electrons_PHI");

	float eleEta = eleEtaVec[idx];
	float elePhi = elePhiVec[idx];


	const vector<float> & pfCandPt = *ev.getVectorFloatAdr("CandIso_Pt");
	const vector<float> & pfCandEta = *ev.getVectorFloatAdr("CandIso_Eta");
	const vector<float> & pfCandPhi = *ev.getVectorFloatAdr("CandIso_Phi");
	const vector<int> & pfCandType = *ev.getVectorIntAdr("CandIso_Type");
	const vector<int> &   pfCandPU = *ev.getVectorIntAdr("CandIso_PU");


	float ChNoPUSumPtDr03 = 0.;
	float ChPUSumPtDr03 = 0.;
	float NhSumPtDr03 = 0.;
	float NhSumPtDr03_Pt1GeV = 0.;

	unsigned nCand = pfCandPt.size();
	for ( unsigned i = 0; i < nCand; ++i ) {
		float dR = deltaR( eleEta, elePhi, pfCandEta[i], pfCandPhi[i] );
		if ( dR < 0.3 ) {

			// charged hadrons
			if(pfCandType[i] > 0.99 && pfCandType[i] < 1.01) {
				// get only candidates from PFNoPU
				// barrel no veto
				if(fabs(eleEta)  < 1.4442) {
					if(pfCandPU[i] == 0) 
						ChNoPUSumPtDr03 += pfCandPt[i];
					else
						ChPUSumPtDr03 += pfCandPt[i];
				}
				else {
					// veto on the endcap
					if(dR > 0.015) {
						if(pfCandPU[i] == 0) 
							ChNoPUSumPtDr03 += pfCandPt[i];
						else
							ChPUSumPtDr03 += pfCandPt[i];
					}
				}		  
			}


			// photons
			if(pfCandType[i] > 3.99 && pfCandType[i] < 4.01) {
				if(fabs(eleEta)  < 1.4442) {
					NhSumPtDr03 += pfCandPt[i];
					if(pfCandPt[i] > 1.) 
						NhSumPtDr03_Pt1GeV += pfCandPt[i];
				}
				else {
					if(dR > 0.08) {
						NhSumPtDr03 += pfCandPt[i];
						if(pfCandPt[i] > 1.) 
							NhSumPtDr03_Pt1GeV += pfCandPt[i];
					}
				}
			}

			// neutral hadrons
			if(pfCandType[i] > 4.99 && pfCandType[i] < 5.01) {
				NhSumPtDr03 +=  pfCandPt[i];
				if(pfCandPt[i] > 1.) 
					NhSumPtDr03_Pt1GeV += pfCandPt[i];
			}

			// to decide if add electrons and muons... probably to add electrons
			// some extra condition is needed.

		}
	}


	// corrections got from DATA

	float Rho03 = rho * M_PI * 0.3 * 0.3;
	float EA_barrel =  (0.054667  / 0.212706) * Rho03;
	float EA_endcap =  (0.0558634 / 0.215981) * Rho03;

	// 	float EA_barrel =  (0.056234 / 0.751912) * RhoCorrectionFactor;
	// 	float EA_endcap =  (0.051635 / 0.758839) * RhoCorrectionFactor;

	float DB_barrel =  (0.054667 / 0.151236) * ChPUSumPtDr03;
	float DB_endcap =  (0.0558634 / 0.157923) * ChPUSumPtDr03;

	float NhSumPtDr03_EACorr = 0.;
	float NhSumPtDr03_DBCorr = 0.;

	if(fabs(eleEta) < 1.4442) {	  
		NhSumPtDr03_EACorr = NhSumPtDr03 - EA_barrel;
		NhSumPtDr03_DBCorr = NhSumPtDr03 - DB_barrel;
	}
	else {
		NhSumPtDr03_EACorr = NhSumPtDr03 - EA_endcap;
		NhSumPtDr03_DBCorr = NhSumPtDr03 - DB_endcap;
	}

	// this is under investation 
	if(NhSumPtDr03_EACorr < 0.)
		NhSumPtDr03_EACorr = 0.;
	if(NhSumPtDr03_DBCorr < 0.)
		NhSumPtDr03_DBCorr = 0.;


	iso_nopu        = ChNoPUSumPtDr03 + NhSumPtDr03;
	iso_nopu_eacorr = ChNoPUSumPtDr03 + NhSumPtDr03_EACorr;
	iso_nopu_dbcorr = ChNoPUSumPtDr03 + NhSumPtDr03_DBCorr;
	iso_nopu_pt1GeV = ChNoPUSumPtDr03 + NhSumPtDr03_Pt1GeV;


	const vector<float> & trackIso = *ev.getVectorFloatAdr("Electrons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Electrons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Electrons_HCALISO");

	iso_std = trackIso[idx] + max(ecalIso[idx] - 1, 0) + hcalIso[idx]; 
	iso_std_rhocorr =  iso_std - Rho03;

	return;
}
void calculatePFIsol_mu(const Event & ev, 
		int idx,
		float & iso_nopu,
		float & iso_nopu_eacorr,
		float & iso_nopu_dbcorr,
		float & iso_nopu_pt1GeV,
		float & iso_std,
		float & iso_std_rhocorr) {

	float rho = ev.getSingleVariableValue<float>("RhoCorrectionFactor");
	const vector<float> & muEtaVec = *ev.getVectorFloatAdr("Muons_ETA");
	const vector<float> & muPhiVec = *ev.getVectorFloatAdr("Muons_PHI");

	float muEta = muEtaVec[idx];
	float muPhi = muPhiVec[idx];


	const vector<float> & pfCandPt = *ev.getVectorFloatAdr("CandIso_Pt");
	const vector<float> & pfCandEta = *ev.getVectorFloatAdr("CandIso_Eta");
	const vector<float> & pfCandPhi = *ev.getVectorFloatAdr("CandIso_Phi");
	const vector<int> & pfCandType = *ev.getVectorIntAdr("CandIso_Type");
	const vector<int> &   pfCandPU = *ev.getVectorIntAdr("CandIso_PU");



	float ChNoPUSumPtDr03 = 0.;
	float ChPUSumPtDr03 = 0.;
	float NhSumPtDr03 = 0.;
	float NhSumPtDr03_Pt1GeV = 0.;

	unsigned nCand = pfCandPt.size();
	for ( unsigned i = 0; i < nCand; ++i ) {
		float dR = deltaR( muEta, muPhi, pfCandEta[i], pfCandPhi[i] );
		if ( dR < 0.3 ) {

			// charged hadrons
			if(pfCandType[i] > 0.99 && pfCandType[i] < 1.01) {
				// get only candidates from PFNoPU
				if(pfCandPU[i] == 0) 
					ChNoPUSumPtDr03 += pfCandPt[i];
				else
					ChPUSumPtDr03 += pfCandPt[i];
			}


			// photons
			if(pfCandType[i] > 3.99 && pfCandType[i] < 4.01) {
				NhSumPtDr03 += pfCandPt[i];
				if(pfCandPt[i] > 1.) 
					NhSumPtDr03_Pt1GeV += pfCandPt[i];
			}

			// neutral hadrons
			if(pfCandType[i] > 4.99 && pfCandType[i] < 5.01) {
				NhSumPtDr03 +=  pfCandPt[i];
				if(pfCandPt[i] > 1.) 
					NhSumPtDr03_Pt1GeV += pfCandPt[i];
			}

			// to decide if adds electrons and muons... probably to add muons
			// some extra condition is needed.

		}
	}


	// corrections got from DATA

	// NOTE: The EA and DB paramenters are needed to be recomputed for Muons 


	float Rho03 = rho * M_PI  * 0.3 * 0.3;
	float EA_barrel =  (0.054667  / 0.212706) * Rho03;
	float EA_endcap =  (0.0558634 / 0.215981) * Rho03;

	// 	float EA_barrel =  (0.056234 / 0.751912) * RhoCorrectionFactor;
	// 	float EA_endcap =  (0.051635 / 0.758839) * RhoCorrectionFactor;

	float DB_barrel =  (0.054667 / 0.151236) * ChPUSumPtDr03;
	float DB_endcap =  (0.0558634 / 0.157923) * ChPUSumPtDr03;

	float NhSumPtDr03_EACorr = 0.;
	float NhSumPtDr03_DBCorr = 0.;

	if(fabs(muEta) < 1.4442) {	  
		NhSumPtDr03_EACorr = NhSumPtDr03 - EA_barrel;
		NhSumPtDr03_DBCorr = NhSumPtDr03 - DB_barrel;
	}
	else {
		NhSumPtDr03_EACorr = NhSumPtDr03 - EA_endcap;
		NhSumPtDr03_DBCorr = NhSumPtDr03 - DB_endcap;
	}

	// this is under investation 
	if(NhSumPtDr03_EACorr < 0.)
		NhSumPtDr03_EACorr = 0.;
	if(NhSumPtDr03_DBCorr < 0.)
		NhSumPtDr03_DBCorr = 0.;



	iso_nopu        = ChNoPUSumPtDr03 + NhSumPtDr03;
	iso_nopu_eacorr = ChNoPUSumPtDr03 + NhSumPtDr03_EACorr;
	iso_nopu_dbcorr = ChNoPUSumPtDr03 + NhSumPtDr03_DBCorr;
	iso_nopu_pt1GeV = ChNoPUSumPtDr03 + NhSumPtDr03_Pt1GeV;


	const vector<float> & trackIso = *ev.getVectorFloatAdr("Muons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Muons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Muons_HCALISO");
	iso_std = trackIso[idx] + ecalIso[idx] + hcalIso[idx]; 
	iso_std_rhocorr =  iso_std - Rho03;

	return;
}
