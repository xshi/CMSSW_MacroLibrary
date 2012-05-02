// ROOT Libraries
#include <TFile.h>
#include <TLorentzVector.h>
// Standard Libraries
//#include <cmath>
#include <iostream>
#include <iomanip>
//#include <vector>
//#include <string>
//#include <sstream>
// Other
#include "preselectionCMG.h"
#include "event.h"
#include "options.h"
//#include "triggerinfo.h"
#include "toolbox.h"
#include "muon.h"
#include "electron.h"
#include "toolsCMG.h"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;
using std::setw;

void LeptonPreselectionCMG( const Options & opt, PreselType type ) {
	if (type == ELE)
		cout << "Entering ElectronPreselection() ..." << endl;
	else if (type == MU)
		cout << "Entering MuonPreselection() ..." << endl;

	TFile * file = new TFile( opt.checkStringOption("inpath").c_str() );
	if (!file->IsOpen())
		throw string("ERROR: Can't open the file: " + opt.checkStringOption("inpath") + "!");
	TDirectory * dir = (TDirectory *) file->Get("evAnalyzer");
	TTree * tree = ( TTree * ) dir->Get( "data" );
	Event ev( tree );
	LeptonVariables leptonVars( ev );
	ElectronVariables electronVars( ev );
	MuonVariables muonVars( ev );
	JetVariables jetVars(ev);
	const int * runP = ev.getSVA<int>("run"); 
	const int * lumiP = ev.getSVA<int>("lumi"); 
	const int * eventP = ev.getSVA<int>("event"); 
	const float * metPtA = ev.getAVA<float>("met_pt");
	const float * metPhiA = ev.getAVA<float>("met_phi");
	const float * rhoP = ev.getSVA<float>("rho25");
	const int * nvtxP = ev.getSVA<int>("nvtx"); 
	const int * niP = ev.getSVA<int>("ngenITpu"); 

	string outputFile( opt.checkStringOption("outpath") );
	size_t pos = outputFile.find(".root");
	if ( pos != string::npos )
		outputFile.erase(pos);

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

	// bool isData = opt.checkBoolOption("isData");

	// unsigned long nentries = tree->GetEntries();
	unsigned long nentries = 100000;
	for ( unsigned long i = 0; i < nentries; i++ ) {
		if ( i % 10000 == 0) {
			cout << string(40, '\b');
			cout << setw(10) << i << " / " << setw(10) << nentries << " done ..." << std::flush;
		}

		tree->GetEntry( i );

		run = -999;
		lumi = -999;
		event = -999;
		pfmet = -999;
		nele = -999;
		nmu = -999;
		nsoftmu = -999;
		l1pt = -999;
		l1eta = -999;
		l1phi = -999;
		l2pt = -999;
		l2eta = -999;
		l2phi = -999;
		zmass = -999;
		zpt = -999;
		zeta = -999;
		mt = -999;
		njet = -999;
		nsoftjet = -999;
		maxJetBTag = -999;
		minDeltaPhiJetMet = -999;
		minDeltaPhiSoftJetMet = -999;
		nvtx = -999;
		ni = -999;

		run = *runP;
		lumi = *lumiP;
		event = *eventP;

//		cout << run << ":" << lumi << ":" << event << endl;
//		if (isData) {
//			if (!triggerAccept(ev, type))
//				continue;
//		}

		pfmet = metPtA[0];

		vector<Electron> electrons = buildElectronCollection(ev, leptonVars, electronVars);
		vector<Muon> muons = buildMuonCollection(ev, leptonVars, muonVars);

		float rho = *rhoP;

		vector<Electron> selectedElectrons;
		for (unsigned i = 0; i < electrons.size(); ++i) {
			TLorentzVector lv = electrons[i].lorentzVector();
			if ( lv.Pt() > 20 && fabs(lv.Eta()) < 2.5 && !electrons[i].isInCrack() && electrons[i].passesTightID()
					&& electrons[i].isPFIsolatedTight(rho) )
				selectedElectrons.push_back(electrons[i]);
		}

		vector<Electron> looseElectrons;
		for (unsigned i = 0; i < electrons.size(); ++i) {
			TLorentzVector lv = electrons[i].lorentzVector();
			if ( lv.Pt() > 10 && fabs(lv.Eta()) < 2.5 && !electrons[i].isInCrack() && electrons[i].passesTightID()
					&& electrons[i].isPFIsolatedTight(rho) )
				looseElectrons.push_back(electrons[i]);
		}

		vector<Muon> selectedMuons;
		for (unsigned i = 0; i < muons.size(); ++i) {
			TLorentzVector lv = muons[i].lorentzVector();
			if ( lv.Pt() > 20 && fabs(lv.Eta()) < 2.4 && muons[i].isTightMuon()
					&& muons[i].isPFIsolatedTight(Electron::effAreaMC(lv.Eta()), rho) ) {
				selectedMuons.push_back(muons[i]);
			}
		}

		vector<Muon> looseMuons;
		for (unsigned i = 0; i < muons.size(); ++i) {
			TLorentzVector lv = muons[i].lorentzVector();
			if ( lv.Pt() > 10 && fabs(lv.Eta()) < 2.4 && muons[i].isTightMuon()
					&& muons[i].isPFIsolatedTight(Electron::effAreaMC(lv.Eta()), rho) )
				looseMuons.push_back(muons[i]);
		}

		int nLeptons = looseElectrons.size() + looseMuons.size();
		if ( nLeptons > 2 )
			continue;

		string leptonsType;
		Lepton * selectedLeptons[2] = {0};
		if (type == ELE) {
			if (selectedElectrons.size() < 2)
				continue;
		   else {
			   selectedLeptons[0] = &selectedElectrons[0];
			   selectedLeptons[1] = &selectedElectrons[1];
		   }
		} else if (type == MU) {
			if (selectedMuons.size() < 2)
				continue;
			else {
			   selectedLeptons[0] = &selectedMuons[0];
			   selectedLeptons[1] = &selectedMuons[1];
			}
		}

//		if ( selectedLeptons[0]->id != -selectedLeptons[1]->id )
//			continue;

		nele = looseElectrons.size();
		nmu = looseMuons.size();
		nsoftmu = 0;

		TLorentzVector lep1 = selectedLeptons[0]->lorentzVector();
		TLorentzVector lep2 = selectedLeptons[1]->lorentzVector();

		if (lep1.Pt() < 20 || lep2.Pt() < 20)
			continue;

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

		double metphi = metPhiA[0];
		double metPx = pfmet * cos( metphi );
		double metPy = pfmet * sin( metphi );
		double px = metPx + Zcand.Px();
		double py = metPy + Zcand.Py();
		double pt2 = px * px + py * py;
		double e = sqrt(zpt * zpt + zmass * zmass) + sqrt(pfmet * pfmet + zmass * zmass);
		double mt2 = e * e - pt2;
		mt = (mt2 > 0) ? sqrt(mt2) : 0;

		vector<unsigned> jets;
		selectJetsCMG( ev, jetVars, jets );
		njet = jets.size();

		vector<unsigned> softjets;
		selectJetsCMG( ev, jetVars, softjets, 15 );
		nsoftjet = softjets.size();

		const ArrayVariableContainer<float> * j_px = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_px));
		const ArrayVariableContainer<float> * j_py = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_py));
		const ArrayVariableContainer<float> * j_pz = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_pz));
		const ArrayVariableContainer<float> * j_en = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_en));
		const ArrayVariableContainer<float> * j_btag = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_btag));

		maxJetBTag = -999;
		for ( int j = 0; j < njet; ++j ) {
			unsigned jetIdx = jets[j];
			TLorentzVector jet( j_px->getVal(jetIdx), j_py->getVal(jetIdx), j_pz->getVal(jetIdx), j_en->getVal(jetIdx) );
			if ( j_btag->getVal(jetIdx) > maxJetBTag && fabs(jet.Eta()) < 2.4 )
				maxJetBTag = j_btag->getVal(jetIdx);
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

		minDeltaPhiJetMet = 999;
		for ( int j = 0; j < njet; ++j ) {
			unsigned jetIdx = jets[j];
			TLorentzVector jet( j_px->getVal(jetIdx), j_py->getVal(jetIdx), j_pz->getVal(jetIdx), j_en->getVal(jetIdx) );
			double tempDelPhiJetMet = deltaPhi(metphi, jet.Phi());
			if ( tempDelPhiJetMet < minDeltaPhiJetMet )
				minDeltaPhiJetMet = tempDelPhiJetMet;
		}
		minDeltaPhiSoftJetMet = 999;
		for ( int j = 0; j < nsoftjet; ++j ) {
			unsigned jetIdx = softjets[j];
			TLorentzVector jet( j_px->getVal(jetIdx), j_py->getVal(jetIdx), j_pz->getVal(jetIdx), j_en->getVal(jetIdx) );
			double tempDelPhiSoftJetMet = deltaPhi(metphi, jet.Phi());
			if ( tempDelPhiSoftJetMet < minDeltaPhiSoftJetMet )
				minDeltaPhiSoftJetMet = tempDelPhiSoftJetMet;
		}

		nvtx = *nvtxP;

		ni = *niP;
		
		smallTree->Fill();
	}
	cout << endl;
	delete file;
	smallTree->Write("", TObject::kOverwrite);
	delete smallTree;
	delete out;
}

vector<Muon> buildMuonCollection( const Event & ev, const LeptonVariables & leptonVars, const MuonVariables & muonVars ) {
	vector<Muon> muons;
	const ArrayVariableContainer<int> * m_idbits = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(muonVars.m_idbits));
	const ArrayVariableContainer<float> * m_nMatches = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_nMatches));
	const ArrayVariableContainer<float> * m_validMuonHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_validMuonHits));
	
	const SingleVariableContainer<int> * l1_id = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l1_id));
	if (fabs(l1_id->getVal()) == 13) {
		const SingleVariableContainer<float> * px = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_px));
		const SingleVariableContainer<float> * py = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_py));
		const SingleVariableContainer<float> * pz = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_pz));
		const SingleVariableContainer<float> * en = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_en));
		const SingleVariableContainer<float> * ptErr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ptErr));
		const SingleVariableContainer<float> * ecalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ecalIso));
		const SingleVariableContainer<float> * hcalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_hcalIso));
		const SingleVariableContainer<float> * trkIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkIso));
		const SingleVariableContainer<float> * gIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_gIso));
		const SingleVariableContainer<float> * chIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_chIso));
		const SingleVariableContainer<float> * puchIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_puchIso));
		const SingleVariableContainer<float> * nhIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_nhIso));
		const SingleVariableContainer<int> * genid = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l1_genid));
		const SingleVariableContainer<float> * ensf = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ensf));
		const SingleVariableContainer<float> * ensferr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ensferr));
		const SingleVariableContainer<float> * d0 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_d0));
		const SingleVariableContainer<float> * dZ = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_dZ));
		const SingleVariableContainer<float> * trkpt = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkpt));
		const SingleVariableContainer<float> * trketa = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trketa));
		const SingleVariableContainer<float> * trkphi = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkphi));
		const SingleVariableContainer<float> * trkchi2 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkchi2));
		const SingleVariableContainer<float> * trkValidPixelHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkValidPixelHits));
		const SingleVariableContainer<float> * trkValidTrackerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkValidTrackerHits));
		const SingleVariableContainer<float> * trkLostInnerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkLostInnerHits));
		const SingleVariableContainer<int> * pidC = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l1_pid));
		int pid = pidC->getVal();

		Muon tmp(px->getVal(), py->getVal(), pz->getVal(), en->getVal(), ptErr->getVal(), ecalIso->getVal(), hcalIso->getVal(), trkIso->getVal(), gIso->getVal(),
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l1_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), m_idbits->getVal(pid), m_nMatches->getVal(pid), m_validMuonHits->getVal(pid));
		muons.push_back(tmp);
	}
	
	const SingleVariableContainer<int> * l2_id = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l2_id));
	if (fabs(l2_id->getVal()) == 13) {
		const SingleVariableContainer<float> * px = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_px));
		const SingleVariableContainer<float> * py = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_py));
		const SingleVariableContainer<float> * pz = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_pz));
		const SingleVariableContainer<float> * en = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_en));
		const SingleVariableContainer<float> * ptErr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ptErr));
		const SingleVariableContainer<float> * ecalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ecalIso));
		const SingleVariableContainer<float> * hcalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_hcalIso));
		const SingleVariableContainer<float> * trkIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkIso));
		const SingleVariableContainer<float> * gIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_gIso));
		const SingleVariableContainer<float> * chIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_chIso));
		const SingleVariableContainer<float> * puchIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_puchIso));
		const SingleVariableContainer<float> * nhIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_nhIso));
		const SingleVariableContainer<int> * genid = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l2_genid));
		const SingleVariableContainer<float> * ensf = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ensf));
		const SingleVariableContainer<float> * ensferr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ensferr));
		const SingleVariableContainer<float> * d0 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_d0));
		const SingleVariableContainer<float> * dZ = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_dZ));
		const SingleVariableContainer<float> * trkpt = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkpt));
		const SingleVariableContainer<float> * trketa = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trketa));
		const SingleVariableContainer<float> * trkphi = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkphi));
		const SingleVariableContainer<float> * trkchi2 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkchi2));
		const SingleVariableContainer<float> * trkValidPixelHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkValidPixelHits));
		const SingleVariableContainer<float> * trkValidTrackerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkValidTrackerHits));
		const SingleVariableContainer<float> * trkLostInnerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkLostInnerHits));
		const SingleVariableContainer<int> * pidC = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l2_pid));
		int pid = pidC->getVal();

		Muon tmp( px->getVal(), py->getVal(), pz->getVal(), en->getVal(), ptErr->getVal(), ecalIso->getVal(), hcalIso->getVal(), trkIso->getVal(), gIso->getVal(),
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l2_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), m_idbits->getVal(pid), m_nMatches->getVal(pid), m_validMuonHits->getVal(pid));
		muons.push_back(tmp);
	}

	const SingleVariableContainer<int> * lnC = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.ln));
	int ln = lnC->getVal();
	const ArrayVariableContainer<float> * px = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_px));
	const ArrayVariableContainer<float> * py = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_py));
	const ArrayVariableContainer<float> * pz = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_pz));
	const ArrayVariableContainer<float> * en = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_en));
	const ArrayVariableContainer<float> * ptErr = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ptErr));
	const ArrayVariableContainer<float> * ecalIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ecalIso));
	const ArrayVariableContainer<float> * hcalIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_hcalIso));
	const ArrayVariableContainer<float> * trkIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkIso));
	const ArrayVariableContainer<float> * gIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_gIso));
	const ArrayVariableContainer<float> * chIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_chIso));
	const ArrayVariableContainer<float> * puchIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_puchIso));
	const ArrayVariableContainer<float> * nhIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_nhIso));
	const ArrayVariableContainer<int> * id = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(leptonVars.ln_id));
	const ArrayVariableContainer<int> * genid = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(leptonVars.ln_genid));
	const ArrayVariableContainer<float> * ensf = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ensf));
	const ArrayVariableContainer<float> * ensferr = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ensferr));
	const ArrayVariableContainer<float> * d0 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_d0));
	const ArrayVariableContainer<float> * dZ = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_dZ));
	const ArrayVariableContainer<float> * trkpt = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkpt));
	const ArrayVariableContainer<float> * trketa = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trketa));
	const ArrayVariableContainer<float> * trkphi = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkphi));
	const ArrayVariableContainer<float> * trkchi2 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkchi2));
	const ArrayVariableContainer<float> * trkValidPixelHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkValidPixelHits));
	const ArrayVariableContainer<float> * trkValidTrackerHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkValidTrackerHits));
	const ArrayVariableContainer<float> * trkLostInnerHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkLostInnerHits));
	const ArrayVariableContainer<int> * pidC = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(leptonVars.ln_pid));

	for (int i = 0; i < ln; ++i) {
		if (fabs(id->getVal(i)) == 13) {
			int pid = pidC->getVal(i);
			Muon tmp(px->getVal(i), py->getVal(i), pz->getVal(i), en->getVal(i), ptErr->getVal(i), ecalIso->getVal(i), hcalIso->getVal(i), trkIso->getVal(i), gIso->getVal(i),
				chIso->getVal(i), puchIso->getVal(i), nhIso->getVal(i), id->getVal(i), genid->getVal(i), ensf->getVal(i), ensferr->getVal(i), d0->getVal(i), dZ->getVal(i),
				trkpt->getVal(i), trketa->getVal(i), trkphi->getVal(i), trkchi2->getVal(i), trkValidPixelHits->getVal(i), trkValidTrackerHits->getVal(i),
				trkLostInnerHits->getVal(i), m_idbits->getVal(pid), m_nMatches->getVal(pid), m_validMuonHits->getVal(pid));
			muons.push_back(tmp);
		}
	}
	return muons;
}

vector<Electron> buildElectronCollection(const Event & ev, const LeptonVariables & leptonVars, const ElectronVariables & electronVars) {
	vector<Electron> electrons;
	const ArrayVariableContainer<int> * e_idbits = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(electronVars.e_idbits));
	const ArrayVariableContainer<float> * e_hoe = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_hoe));
	const ArrayVariableContainer<float> * e_dphiin = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_dphiin));
	const ArrayVariableContainer<float> * e_detain = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_detain));
	const ArrayVariableContainer<float> * e_sihih = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_sihih));
	const ArrayVariableContainer<float> * e_sipip = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_sipip));
	const ArrayVariableContainer<float> * e_r9 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_r9));
	const ArrayVariableContainer<float> * e_sce = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_sce));
	const ArrayVariableContainer<float> * e_sceta = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_sceta));
	const ArrayVariableContainer<float> * e_scphi = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_scphi));
	const ArrayVariableContainer<float> * e_e2x5max = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_e2x5max));
	const ArrayVariableContainer<float> * e_e1x5 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_e1x5));
	const ArrayVariableContainer<float> * e_e5x5 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_e5x5));
	const ArrayVariableContainer<float> * e_h2te = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_h2te));
	const ArrayVariableContainer<float> * e_h2tebc = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_h2tebc));
	const ArrayVariableContainer<float> * e_ooemoop = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_ooemoop));
	const ArrayVariableContainer<float> * e_fbrem = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_fbrem));
	const ArrayVariableContainer<float> * e_eopin = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_eopin));
	
	const SingleVariableContainer<int> * l1_id = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l1_id));
	if (fabs(l1_id->getVal()) == 11) {
		const SingleVariableContainer<float> * px = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_px));
		const SingleVariableContainer<float> * py = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_py));
		const SingleVariableContainer<float> * pz = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_pz));
		const SingleVariableContainer<float> * en = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_en));
		const SingleVariableContainer<float> * ptErr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ptErr));
		const SingleVariableContainer<float> * ecalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ecalIso));
		const SingleVariableContainer<float> * hcalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_hcalIso));
		const SingleVariableContainer<float> * trkIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkIso));
		const SingleVariableContainer<float> * gIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_gIso));
		const SingleVariableContainer<float> * chIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_chIso));
		const SingleVariableContainer<float> * puchIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_puchIso));
		const SingleVariableContainer<float> * nhIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_nhIso));
		const SingleVariableContainer<int> * genid = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l1_genid));
		const SingleVariableContainer<float> * ensf = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ensf));
		const SingleVariableContainer<float> * ensferr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ensferr));
		const SingleVariableContainer<float> * d0 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_d0));
		const SingleVariableContainer<float> * dZ = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_dZ));
		const SingleVariableContainer<float> * trkpt = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkpt));
		const SingleVariableContainer<float> * trketa = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trketa));
		const SingleVariableContainer<float> * trkphi = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkphi));
		const SingleVariableContainer<float> * trkchi2 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkchi2));
		const SingleVariableContainer<float> * trkValidPixelHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkValidPixelHits));
		const SingleVariableContainer<float> * trkValidTrackerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkValidTrackerHits));
		const SingleVariableContainer<float> * trkLostInnerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_trkLostInnerHits));
		const SingleVariableContainer<int> * pidC = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l1_pid));
		int pid = pidC->getVal();

		Electron tmp(px->getVal(), py->getVal(), pz->getVal(), en->getVal(), ptErr->getVal(), ecalIso->getVal(), hcalIso->getVal(), trkIso->getVal(), gIso->getVal(),
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l1_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), e_idbits->getVal(pid), e_hoe->getVal(pid), e_dphiin->getVal(pid), e_detain->getVal(pid),
				e_sihih->getVal(pid), e_sipip->getVal(pid), e_r9->getVal(pid), e_sce->getVal(pid), e_sceta->getVal(pid), e_scphi->getVal(pid), e_e2x5max->getVal(pid),
				e_e1x5->getVal(pid), e_e5x5->getVal(pid), e_h2te->getVal(pid), e_h2tebc->getVal(pid), e_ooemoop->getVal(pid), e_fbrem->getVal(pid), e_eopin->getVal(pid));
		electrons.push_back(tmp);
	}
	const SingleVariableContainer<int> * l2_id = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l2_id));
	if (fabs(l2_id->getVal()) == 11) {
		const SingleVariableContainer<float> * px = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_px));
		const SingleVariableContainer<float> * py = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_py));
		const SingleVariableContainer<float> * pz = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_pz));
		const SingleVariableContainer<float> * en = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_en));
		const SingleVariableContainer<float> * ptErr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ptErr));
		const SingleVariableContainer<float> * ecalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ecalIso));
		const SingleVariableContainer<float> * hcalIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_hcalIso));
		const SingleVariableContainer<float> * trkIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkIso));
		const SingleVariableContainer<float> * gIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_gIso));
		const SingleVariableContainer<float> * chIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_chIso));
		const SingleVariableContainer<float> * puchIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_puchIso));
		const SingleVariableContainer<float> * nhIso = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_nhIso));
		const SingleVariableContainer<int> * genid = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l2_genid));
		const SingleVariableContainer<float> * ensf = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ensf));
		const SingleVariableContainer<float> * ensferr = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ensferr));
		const SingleVariableContainer<float> * d0 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_d0));
		const SingleVariableContainer<float> * dZ = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_dZ));
		const SingleVariableContainer<float> * trkpt = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkpt));
		const SingleVariableContainer<float> * trketa = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trketa));
		const SingleVariableContainer<float> * trkphi = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkphi));
		const SingleVariableContainer<float> * trkchi2 = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkchi2));
		const SingleVariableContainer<float> * trkValidPixelHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkValidPixelHits));
		const SingleVariableContainer<float> * trkValidTrackerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkValidTrackerHits));
		const SingleVariableContainer<float> * trkLostInnerHits = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_trkLostInnerHits));
		const SingleVariableContainer<int> * pidC = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.l2_pid));
		int pid = pidC->getVal();

		Electron tmp(px->getVal(), py->getVal(), pz->getVal(), en->getVal(), ptErr->getVal(), ecalIso->getVal(), hcalIso->getVal(), trkIso->getVal(), gIso->getVal(),
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l2_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), e_idbits->getVal(pid), e_hoe->getVal(pid), e_dphiin->getVal(pid), e_detain->getVal(pid),
				e_sihih->getVal(pid), e_sipip->getVal(pid), e_r9->getVal(pid), e_sce->getVal(pid), e_sceta->getVal(pid), e_scphi->getVal(pid), e_e2x5max->getVal(pid),
				e_e1x5->getVal(pid), e_e5x5->getVal(pid), e_h2te->getVal(pid), e_h2tebc->getVal(pid), e_ooemoop->getVal(pid), e_fbrem->getVal(pid), e_eopin->getVal(pid));
		electrons.push_back(tmp);
	}

	const SingleVariableContainer<int> * lnC = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(leptonVars.ln));
	int ln = lnC->getVal();
	const ArrayVariableContainer<float> * px = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_px));
	const ArrayVariableContainer<float> * py = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_py));
	const ArrayVariableContainer<float> * pz = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_pz));
	const ArrayVariableContainer<float> * en = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_en));
	const ArrayVariableContainer<float> * ptErr = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ptErr));
	const ArrayVariableContainer<float> * ecalIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ecalIso));
	const ArrayVariableContainer<float> * hcalIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_hcalIso));
	const ArrayVariableContainer<float> * trkIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkIso));
	const ArrayVariableContainer<float> * gIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_gIso));
	const ArrayVariableContainer<float> * chIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_chIso));
	const ArrayVariableContainer<float> * puchIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_puchIso));
	const ArrayVariableContainer<float> * nhIso = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_nhIso));
	const ArrayVariableContainer<int> * id = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(leptonVars.ln_id));
	const ArrayVariableContainer<int> * genid = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(leptonVars.ln_genid));
	const ArrayVariableContainer<float> * ensf = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ensf));
	const ArrayVariableContainer<float> * ensferr = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ensferr));
	const ArrayVariableContainer<float> * d0 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_d0));
	const ArrayVariableContainer<float> * dZ = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_dZ));
	const ArrayVariableContainer<float> * trkpt = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkpt));
	const ArrayVariableContainer<float> * trketa = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trketa));
	const ArrayVariableContainer<float> * trkphi = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkphi));
	const ArrayVariableContainer<float> * trkchi2 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkchi2));
	const ArrayVariableContainer<float> * trkValidPixelHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkValidPixelHits));
	const ArrayVariableContainer<float> * trkValidTrackerHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkValidTrackerHits));
	const ArrayVariableContainer<float> * trkLostInnerHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_trkLostInnerHits));
	const ArrayVariableContainer<int> * pidC = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(leptonVars.ln_pid));

	for (int i = 0; i < ln; ++i) {
		if (fabs(id->getVal(i)) == 11) {
			int pid = pidC->getVal(i);
			Electron tmp(px->getVal(i), py->getVal(i), pz->getVal(i), en->getVal(i), ptErr->getVal(i), ecalIso->getVal(i), hcalIso->getVal(i), trkIso->getVal(i), gIso->getVal(i),
				chIso->getVal(i), puchIso->getVal(i), nhIso->getVal(i), id->getVal(i), genid->getVal(i), ensf->getVal(i), ensferr->getVal(i), d0->getVal(i), dZ->getVal(i),
				trkpt->getVal(i), trketa->getVal(i), trkphi->getVal(i), trkchi2->getVal(i), trkValidPixelHits->getVal(i), trkValidTrackerHits->getVal(i),
				trkLostInnerHits->getVal(i), e_idbits->getVal(pid), e_hoe->getVal(pid),
					e_dphiin->getVal(pid), e_detain->getVal(pid),	e_sihih->getVal(pid), e_sipip->getVal(pid), e_r9->getVal(pid),
					e_sce->getVal(pid), e_sceta->getVal(pid), e_scphi->getVal(pid), e_e2x5max->getVal(pid), e_e1x5->getVal(pid),
					e_e5x5->getVal(pid), e_h2te->getVal(pid), e_h2tebc->getVal(pid), e_ooemoop->getVal(pid), e_fbrem->getVal(pid),
					e_eopin->getVal(pid));
			electrons.push_back(tmp);
		}
	}
	return electrons;
}

void selectJetsCMG(const Event & ev, const JetVariables & jetVars, vector<unsigned> & jets, double ptMin, double etaMax) {

	const SingleVariableContainer<int> * jn = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(jetVars.jn));
	const ArrayVariableContainer<float> * j_px = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_px));
	const ArrayVariableContainer<float> * j_py = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_py));
	const ArrayVariableContainer<float> * j_pz = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_pz));
	const ArrayVariableContainer<float> * j_en = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_en));

	jets.clear();
	for ( int i = 0; i < jn->getVal(); ++i ) {
		TLorentzVector jet(j_px->getVal(i), j_py->getVal(i), j_pz->getVal(i), j_en->getVal(i));
		if ( jet.Pt() > ptMin && fabs(jet.Eta()) < etaMax )
			jets.push_back( i );
	}
}

/*
void selectPhotons(const Event & ev, vector<unsigned> & photons) {
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
*/
