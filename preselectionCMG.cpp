#include <algorithm>
#include "electron.h"
#include "event.h"
#include "eventPrinter.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include "jet.h"
#include "jet.h"
#include "muon.h"
#include "options.h"
#include "photon.h"
#include "photonPrescale.h"
#include "preselectionCMG.h"
#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include "toolbox.h"
#include "toolsCMG.h"
#include <TRandom.h>

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;
using std::setw;
using std::ofstream;

void LeptonPreselectionCMG( const Options & opt, PreselType type, RooWorkspace * w ) {
	if (type == ELE)
		cout << "Running Electron Preselection :" << endl;
	else if (type == MU)
		cout << "Running Muon Preselection :" << endl;
	else if (type == EMU)
		cout << "Running Electron-Muon Preselection() ..." << endl;
	else if (type == PHOT)
		cout << "Running Photon Preselection :" << endl;

	string inputDir = opt.checkStringOption("INPUT_DIR");
	string outputDir = opt.checkStringOption("OUTPUT_DIR");
	string sampleName = opt.checkStringOption("SAMPLE_NAME");
	string inputFile = inputDir + '/' + sampleName + ".root";
	cout << "\tInput file: " << inputFile << endl;

	TFile * file = new TFile( inputFile.c_str() );
	if (!file->IsOpen())
		throw string("ERROR: Can't open the file: " + inputFile + "!");
	TDirectory * dir = (TDirectory *) file->Get("evAnalyzer");
	TH1D * histo = (TH1D *) ((TDirectory *) dir->Get("h2zz"))->Get("cutflow");
	TTree * tree = ( TTree * ) dir->Get( "data" );
	Event ev( tree );
	LeptonVariables leptonVars( ev );
	ElectronVariables electronVars( ev );
	MuonVariables muonVars( ev );
	JetVariables jetVars(ev);
	PhotonVariables photonVars(ev);
	const int * runP = ev.getSVA<int>("run"); 
	const int * lumiP = ev.getSVA<int>("lumi"); 
	const int * eventP = ev.getSVA<int>("event"); 
	const bool * trigP = ev.getSVA<bool>("hasTrigger");
	const float * metPtA = ev.getAVA<float>("met_pt");
	const float * metPhiA = ev.getAVA<float>("met_phi");
	//const float * rhoP = ev.getSVA<float>("rho25");
	const float * rhoP = ev.getSVA<float>("rho");
	const int * nvtxP = ev.getSVA<int>("nvtx"); 
	const int * niP = ev.getSVA<int>("ngenITpu"); 
	const int * cat = ev.getSVA<int>("cat"); 
	const int * phPrescale = 0;

	EventPrinter evPrint(ev, "events.txt");
	//evPrint.readInEvents("output1.txt");
	//evPrint.printElectrons();
	//evPrint.printMuons();
	//evPrint.printHeader();

	string outputFile = outputDir + '/' + sampleName;

	if (type == ELE)
		outputFile += "_elePresel.root";
	else if (type == MU)
		outputFile += "_muPresel.root";
	else if (type == EMU)
		outputFile += "_emuPresel.root";
	else if (type == PHOT)
		outputFile += "_phPresel.root";
	cout << "\tOutput file: " << outputFile << endl;

	TFile * out = new TFile( outputFile.c_str(), "recreate" );
	TH1D * outHisto = new TH1D("nevt", "nevt", 1, 0, 1);
	outHisto->SetBinContent(1, histo->GetBinContent(1));
	outHisto->Write("nevt");

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
	int nsoftjet;
	int nhardjet;
	double maxJetBTag;
	double minDeltaPhiJetMet;
	double detajj;
	double mjj;
	int nvtx;
	int ni;
	int category;
	double weight;

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
	smallTree->Branch( "NSOFTJET", &nsoftjet, "NSOFTJET/I" );
	smallTree->Branch( "NHARDJET", &nhardjet, "NHARDJET/I" );
	smallTree->Branch( "MAXJETBTAG", &maxJetBTag, "MAXJETBTAG/D" );
	smallTree->Branch( "MINDPJETMET", &minDeltaPhiJetMet, "MINDPJETMET/D" );
	smallTree->Branch( "DETAJJ", &detajj, "DETAJJ/D" );
	smallTree->Branch( "MJJ", &mjj, "MJJ/D" );
	smallTree->Branch( "NVTX", &nvtx, "NVTX/I" );
	smallTree->Branch( "nInter" , &ni, "nInter/I" );
	smallTree->Branch( "CATEGORY", &category, "CATEGORY/I" );
	smallTree->Branch( "Weight" , &weight, "Weight/D" );

	bool isData = opt.checkBoolOption("DATA");

	unsigned long nentries = tree->GetEntries();

	RooDataSet * events = nullptr;

	PhotonPrescale photonPrescales;

	if (type == PHOT) {
		if (w == nullptr)
			throw string("ERROR: No mass peak pdf!");
		RooRealVar * zmass = w->var("mass");
		zmass->setRange(76.2, 106.2);
		RooAbsPdf * pdf = w->pdf("massPDF");
		events = pdf->generate(*zmass, nentries);

		photonPrescales.addTrigger("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly", 22);
		photonPrescales.addTrigger("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly", 36);
		photonPrescales.addTrigger("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly", 50);
		photonPrescales.addTrigger("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly", 75);
		photonPrescales.addTrigger("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly", 90);
		photonPrescales.addTrigger("HLT_Photon135", 135);
		photonPrescales.addTrigger("HLT_Photon150", 150);
		photonPrescales.addTrigger("HLT_Photon160", 160);
		//photonPrescales.addTrigger( "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly", 22, opt.checkStringOption("PHOTON_PRESCALE_DIRECTORY") + "/HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PhotonPrescales.txt" );
		//photonPrescales.addTrigger( "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly", 36, opt.checkStringOption("PHOTON_PRESCALE_DIRECTORY") + "/HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PhotonPrescales.txt" );
		//photonPrescales.addTrigger( "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly", 50, opt.checkStringOption("PHOTON_PRESCALE_DIRECTORY") + "/HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PhotonPrescales.txt" );
		//photonPrescales.addTrigger( "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly", 75, opt.checkStringOption("PHOTON_PRESCALE_DIRECTORY") + "/HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PhotonPrescales.txt" );
		//photonPrescales.addTrigger( "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly", 90, opt.checkStringOption("PHOTON_PRESCALE_DIRECTORY") + "/HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PhotonPrescales.txt" );

		phPrescale = ev.getSVA<int>("gn_prescale");
	}

	const double TRIGGERTHR = 160;
	TH1D histoPT("histoPT", "histoPT", 100, 150, 200);
	TH1D eff("eff", "eff", 100, 150, 200);
	TH1D prescale("prescale", "prescale", 100, 0, 500);
	TH1D outRun("runs", "runs", 100, 194000, 204000);

	for ( unsigned long iEvent = 0; iEvent < nentries; iEvent++ ) {

		if ( iEvent % 10000 == 0) {
			cout << string(40, '\b');
			cout << setw(10) << iEvent << " / " << setw(10) << nentries << " done ..." << std::flush;
		}

		tree->GetEntry( iEvent );

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
		nsoftjet = -999;
		nhardjet = -999;
		maxJetBTag = -999;
		minDeltaPhiJetMet = -999;
		detajj = -999;
		mjj = -999;
		nvtx = -999;
		ni = -999;
		weight = -999;
		category = -1;

		run = *runP;
		lumi = *lumiP;
		event = *eventP;

		if (type == ELE && (*cat) != 2) {
			continue;
		}
		if (type == MU && (*cat) != 1) {
			continue;
		}
		if (type == EMU && (*cat) != 3) {
			continue;
		}

		if (! *trigP) {
			continue;
		}
		//		cout << run << ":" << lumi << ":" << event << endl;
		//		if (isData) {
		//			if (!triggerAccept(ev, type))
		//				continue;
		//		}

		vector<Electron> electrons = buildElectronCollection(ev, leptonVars, electronVars);
		vector<Muon> muons = buildMuonCollection(ev, leptonVars, muonVars);

		float rho = *rhoP;

		vector<Electron> looseElectrons;
		vector<Electron> selectedElectrons;
		for (unsigned j = 0; j < electrons.size(); ++j) {
			TLorentzVector lv = electrons[j].lorentzVector();
			if ( lv.Pt() > 10 && fabs(lv.Eta()) < 2.5 && !electrons[j].isInCrack()
					&& electrons[j].passesVetoID()
					&& electrons[j].pfIsolation(rho, isData) < 0.15
					) {
				looseElectrons.push_back(electrons[j]);
			}
			if ( lv.Pt() > 20 && fabs(lv.Eta()) < 2.5 && !electrons[j].isInCrack()
					&& electrons[j].passesMediumID()
					&& electrons[j].pfIsolation(rho, isData) < 0.15
					) {
				selectedElectrons.push_back(electrons[j]);
			}
		}

		vector<Muon> looseMuons;
		vector<Muon> softMuons;
		vector<Muon> selectedMuons;
		for (unsigned j = 0; j < muons.size(); ++j) {
			TLorentzVector lv = muons[j].lorentzVector();
			if ( lv.Pt() > 10 && fabs(lv.Eta()) < 2.4 && muons[j].isLooseMuon() && muons[j].isPFIsolatedLoose() ) {
				looseMuons.push_back(muons[j]);
			} else if ( lv.Pt() > 3 && fabs(lv.Eta()) < 2.4 && muons[j].isSoftMuon() )
				softMuons.push_back(muons[j]);
			if ( lv.Pt() > 20 && fabs(lv.Eta()) < 2.4 && muons[j].isTightMuon() && muons[j].isPFIsolatedTight() ) {
				selectedMuons.push_back(muons[j]);
			}
		}

		vector<Photon> photons = selectPhotonsCMG( ev, photonVars );
		vector<Photon> selectedPhotons;
		for (unsigned i = 0; i < photons.size(); ++i) {
			if (photons[i].isSelected(rho))
				selectedPhotons.push_back( photons[i] );
		}

		if (type == PHOT) {
			vector<Electron> tmpElectrons;
			for (unsigned i = 0; i < selectedPhotons.size(); ++i) {
				TLorentzVector phVec = selectedPhotons[i].lorentzVector();
				for (unsigned j = 0; j < looseElectrons.size(); ++j) {
					TLorentzVector elVec = looseElectrons[j].lorentzVector();
					double dR = deltaR(phVec.Eta(), phVec.Phi(), elVec.Eta(), elVec.Phi());
					if ( dR > 0.05 )
						tmpElectrons.push_back( looseElectrons[j] );
				}
			}
			looseElectrons = tmpElectrons;
		}

		string leptonsType;
		Lepton * selectedLeptons[2] = {0};
		if (type == ELE) {
			if (selectedElectrons.size() < 2) {
				continue;
			} else {
				selectedLeptons[0] = &selectedElectrons[0];
				selectedLeptons[1] = &selectedElectrons[1];
			}
		} else if (type == MU) {
			if (selectedMuons.size() < 2) {
				continue;
			} else {
				selectedLeptons[0] = &selectedMuons[0];
				selectedLeptons[1] = &selectedMuons[1];
			}
		} else if (type == EMU) {
			if (selectedElectrons.size() < 1 || selectedMuons.size() < 1) {
				continue;
			} else {
				selectedLeptons[0] = &selectedElectrons[0];
				selectedLeptons[1] = &selectedMuons[0];
			}
		} else if (type == PHOT) {
			if (selectedPhotons.size() < 1) {
				continue;
			}
		}

		nele = looseElectrons.size();
		nmu = looseMuons.size();
		nsoftmu = softMuons.size();

		TLorentzVector Zcand;

		if (type == ELE || type == MU || type == EMU) {
			TLorentzVector lep1 = selectedLeptons[0]->lorentzVector();
			TLorentzVector lep2 = selectedLeptons[1]->lorentzVector();

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

			Zcand = lep1 + lep2;
			zmass = Zcand.M();
		} else if (type == PHOT) {
			Zcand = selectedPhotons[0].lorentzVector();
			zmass = events->get(iEvent)->getRealValue("mass");
		}

		zpt = Zcand.Pt();
		zeta = Zcand.Eta();

		if (type == PHOT) {
			if ( (*cat) < 10) {
				continue;
			} else {
				int trgThreshold = ((*cat) - 22) / 1000;
				weight = *phPrescale;
				double nextT = photonPrescales.nextThreshold(trgThreshold);
				eff.Fill(zpt, weight);
				if (trgThreshold == TRIGGERTHR) {
					histoPT.Fill(zpt, weight);
					prescale.Fill(weight);
					outRun.Fill(run);
				}
				if (nextT > 0 && zpt < nextT)
					continue;
			}
		}

		vector<Jet> jets = selectJetsCMG( ev, jetVars );
		if (type == PHOT) {
			vector<Jet> tmpJets;
			for (unsigned i = 0; i < selectedPhotons.size(); ++i) {
				TLorentzVector phVec = selectedPhotons[i].lorentzVector();
				for (unsigned j = 0; j < jets.size(); ++j) {
					TLorentzVector jVec = jets[j].lorentzVector();
					double dR = deltaR(phVec.Eta(), phVec.Phi(), jVec.Eta(), jVec.Phi());
					if ( dR > 0.4 )
						tmpJets.push_back( jets[j] );
				}
			}
			jets = tmpJets;
		}

		TLorentzVector jetDiff = smearJets( jets );
		if (isData && jetDiff != TLorentzVector())
			throw std::string("Jet Corrections different from zero in DATA!");

		TLorentzVector met;
		met.SetPtEtaPhiM(metPtA[2], 0.0, metPhiA[2], 0.0);
		met -= jetDiff;
		pfmet = met.Pt();

		double px = met.Px() + Zcand.Px();
		double py = met.Py() + Zcand.Py();
		double pt2 = px * px + py * py;
		double e = sqrt(zpt * zpt + zmass * zmass) + sqrt(pfmet * pfmet + zmass * zmass);
		double mt2 = e * e - pt2;
		mt = (mt2 > 0) ? sqrt(mt2) : 0;

		vector<Jet> hardjets;
		vector<Jet> softjets;
		maxJetBTag = -999;
		minDeltaPhiJetMet = 999;
		for ( unsigned j = 0; j < jets.size(); ++j ) {
			TLorentzVector jet = jets[j].lorentzVector();
			if ( jet.Pt() > 30 ) {
				hardjets.push_back( jets[j] );
			}
			if ( jet.Pt() > 15 )
				softjets.push_back( jets[j] );
		}
		nhardjet = hardjets.size();
		nsoftjet = softjets.size();
		if ( type == PHOT && nsoftjet == 0 )
			continue;

		if (nhardjet > 1) {
			sort(hardjets.begin(), hardjets.end(), [](const Jet & a, const Jet & b) {
					return a.lorentzVector().Eta() < b.lorentzVector().Eta();
				});
			for (unsigned j = 0; j < hardjets.size() - 1; ++j) {
				TLorentzVector jet1 = hardjets[j].lorentzVector();
				TLorentzVector jet2 = hardjets[j + 1].lorentzVector();
				double tmpDelEta = jet2.Eta() - jet1.Eta();
				TLorentzVector diJetSystem = jet1 + jet2;
				double tmpMass = diJetSystem.M();
				if (tmpDelEta > 4.0 && tmpMass > 500 && zeta > jet1.Eta() && jet2.Eta() > zeta) {
					detajj = tmpDelEta;
					mjj = tmpMass;
				}
			}
		}

		category = evCategory(nhardjet, nsoftjet, detajj, mjj, type == PHOT);

		minDeltaPhiJetMet = 10;
		if (category == 1) {
			for ( unsigned j = 0; j < softjets.size(); ++j ) {
				TLorentzVector jet = softjets[j].lorentzVector();
				double tempDelPhiJetMet = deltaPhi(met.Phi(), jet.Phi());
				if ( tempDelPhiJetMet < minDeltaPhiJetMet )
					minDeltaPhiJetMet = tempDelPhiJetMet;
			}
		} else {
			for ( unsigned j = 0; j < hardjets.size(); ++j ) {
				TLorentzVector jet = hardjets[j].lorentzVector();
				if ( hardjets[j].btag > maxJetBTag && fabs(jet.Eta()) < 2.4 )
					maxJetBTag = hardjets[j].btag;
				double tempDelPhiJetMet = deltaPhi(met.Phi(), jet.Phi());
				if ( tempDelPhiJetMet < minDeltaPhiJetMet )
					minDeltaPhiJetMet = tempDelPhiJetMet;
			}
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

		nvtx = *nvtxP;

		ni = *niP;

		if ( opt.checkBoolOption("ADDITIONAL_LEPTON_VETO") && (type == ELE || type == MU || type == EMU) && ((nele + nmu + nsoftmu) > 2) )
			continue;
		if ( opt.checkBoolOption("ADDITIONAL_LEPTON_VETO") && (type == PHOT) && ((nele + nmu + nsoftmu) > 0) )
			continue;
		if ( opt.checkBoolOption("ZPT_CUT") && zpt < 30 )
			continue;
		// for different background estimation methods different window should be applied:
		// * sample for photons should have 76.2 < zmass < 106.2
		// * sample for non-resonant background should not have this cut applied
		if ( opt.checkBoolOption("TIGHT_ZMASS_CUT") && (type == ELE || type == MU) && (zmass < 76.2 || zmass > 106.2))
			continue;
		if ( opt.checkBoolOption("WIDE_ZMASS_CUT") && (type == ELE || type == MU) && (zmass < 76.2 || zmass > 106.2))
			continue;
		if ( opt.checkBoolOption("BTAG_CUT") && ( maxJetBTag > 0.275) )
			continue;
		
		smallTree->Fill();

		evPrint.setElectronCollection(electrons);
		evPrint.setMuonCollection(muons);
		evPrint.print();
	}
	cout << endl;
	
	TCanvas canv("canv", "canv", 800, 600);
	histoPT.Sumw2();
	eff.Sumw2();
	histoPT.Divide(&eff);
	histoPT.Draw();
	canv.SaveAs("threshold.ps");
	canv.Clear();
	prescale.Draw();
	canv.SaveAs("prescale.ps");
	canv.Clear();
	outRun.Draw();
	canv.SaveAs("runs.ps");

	delete file;
	smallTree->Write("", TObject::kOverwrite);
	delete smallTree;
	delete out;
}

vector<Muon> buildMuonCollection( const Event & ev, const LeptonVariables & leptonVars, const MuonVariables & muonVars ) {
	vector<Muon> muons;
	const ArrayVariableContainer<int> * m_idbits = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(muonVars.m_idbits));
	const ArrayVariableContainer<float> * m_nMatches = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_nMatches));
	const ArrayVariableContainer<float> * m_nMatchedStations = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_nMatchedStations));
	const ArrayVariableContainer<float> * m_validMuonHits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_validMuonHits));
	const ArrayVariableContainer<float> * m_innerTrackChi2 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_innerTrackChi2));
	const ArrayVariableContainer<float> * m_trkLayersWithMeasurement = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_trkLayersWithMeasurement));
	const ArrayVariableContainer<float> * m_pixelLayersWithMeasurement = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(muonVars.m_pixelLayersWithMeasurement));
	
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
		const SingleVariableContainer<float> * ip3d = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ip3d));
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
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l1_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(), ip3d->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), m_idbits->getVal(pid), m_nMatches->getVal(pid), m_nMatchedStations->getVal(pid), m_validMuonHits->getVal(pid), m_innerTrackChi2->getVal(pid),
				m_trkLayersWithMeasurement->getVal(pid), m_pixelLayersWithMeasurement->getVal(pid));
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
		const SingleVariableContainer<float> * ip3d = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ip3d));
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
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l2_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(), ip3d->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), m_idbits->getVal(pid), m_nMatches->getVal(pid), m_nMatchedStations->getVal(pid), m_validMuonHits->getVal(pid), m_innerTrackChi2->getVal(pid),
				m_trkLayersWithMeasurement->getVal(pid), m_pixelLayersWithMeasurement->getVal(pid));
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
	const ArrayVariableContainer<float> * ip3d = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ip3d));
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
				chIso->getVal(i), puchIso->getVal(i), nhIso->getVal(i), id->getVal(i), genid->getVal(i), ensf->getVal(i), ensferr->getVal(i), d0->getVal(i), dZ->getVal(i), ip3d->getVal(i),
				trkpt->getVal(i), trketa->getVal(i), trkphi->getVal(i), trkchi2->getVal(i), trkValidPixelHits->getVal(i), trkValidTrackerHits->getVal(i),
				trkLostInnerHits->getVal(i), m_idbits->getVal(pid), m_nMatches->getVal(pid), m_nMatchedStations->getVal(pid), m_validMuonHits->getVal(pid), m_innerTrackChi2->getVal(pid),
				m_trkLayersWithMeasurement->getVal(pid), m_pixelLayersWithMeasurement->getVal(pid));
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
	const ArrayVariableContainer<float> * e_dEtaCalo = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_dEtaCalo));
	const ArrayVariableContainer<float> * e_kfchi2 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_kfchi2));
	const ArrayVariableContainer<float> * e_kfhits = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_kfhits));
	const ArrayVariableContainer<float> * e_etawidth = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_etawidth));
	const ArrayVariableContainer<float> * e_phiwidth = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_phiwidth));
	const ArrayVariableContainer<float> * e_e1x5e5x5 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_e1x5e5x5));
	const ArrayVariableContainer<float> * e_preShowerOverRaw = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_preShowerOverRaw));
	const ArrayVariableContainer<float> * e_eopout = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(electronVars.e_eopout));
	
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
		const SingleVariableContainer<float> * ip3d = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l1_ip3d));
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
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l1_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(), ip3d->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), e_idbits->getVal(pid), e_hoe->getVal(pid), e_dphiin->getVal(pid), e_detain->getVal(pid),
				e_sihih->getVal(pid), e_sipip->getVal(pid), e_r9->getVal(pid), e_sce->getVal(pid), e_sceta->getVal(pid), e_scphi->getVal(pid), e_e2x5max->getVal(pid),
				e_e1x5->getVal(pid), e_e5x5->getVal(pid), e_h2te->getVal(pid), e_h2tebc->getVal(pid), e_ooemoop->getVal(pid), e_fbrem->getVal(pid),
				e_eopin->getVal(pid), e_dEtaCalo->getVal(pid), e_kfchi2->getVal(pid), e_kfhits->getVal(pid), e_etawidth->getVal(pid), e_phiwidth->getVal(pid),
				e_e1x5e5x5->getVal(pid), e_preShowerOverRaw->getVal(pid), e_eopout->getVal(pid) );
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
		const SingleVariableContainer<float> * ip3d = dynamic_cast<const SingleVariableContainer<float> *>(ev.getVariable(leptonVars.l2_ip3d));
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
				chIso->getVal(), puchIso->getVal(), nhIso->getVal(), l1_id->getVal(), genid->getVal(), ensf->getVal(), ensferr->getVal(), d0->getVal(), dZ->getVal(), ip3d->getVal(),
				trkpt->getVal(), trketa->getVal(), trkphi->getVal(), trkchi2->getVal(), trkValidPixelHits->getVal(), trkValidTrackerHits->getVal(),
				trkLostInnerHits->getVal(), e_idbits->getVal(pid), e_hoe->getVal(pid), e_dphiin->getVal(pid), e_detain->getVal(pid),
				e_sihih->getVal(pid), e_sipip->getVal(pid), e_r9->getVal(pid), e_sce->getVal(pid), e_sceta->getVal(pid), e_scphi->getVal(pid), e_e2x5max->getVal(pid),
				e_e1x5->getVal(pid), e_e5x5->getVal(pid), e_h2te->getVal(pid), e_h2tebc->getVal(pid), e_ooemoop->getVal(pid), e_fbrem->getVal(pid),
				e_eopin->getVal(pid), e_dEtaCalo->getVal(pid), e_kfchi2->getVal(pid), e_kfhits->getVal(pid), e_etawidth->getVal(pid), e_phiwidth->getVal(pid),
				e_e1x5e5x5->getVal(pid), e_preShowerOverRaw->getVal(pid), e_eopout->getVal(pid) );
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
	const ArrayVariableContainer<float> * ip3d = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(leptonVars.ln_ip3d));
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
				chIso->getVal(i), puchIso->getVal(i), nhIso->getVal(i), id->getVal(i), genid->getVal(i), ensf->getVal(i), ensferr->getVal(i), d0->getVal(i), dZ->getVal(i), ip3d->getVal(i),
				trkpt->getVal(i), trketa->getVal(i), trkphi->getVal(i), trkchi2->getVal(i), trkValidPixelHits->getVal(i), trkValidTrackerHits->getVal(i),
				trkLostInnerHits->getVal(i), e_idbits->getVal(pid), e_hoe->getVal(pid),
					e_dphiin->getVal(pid), e_detain->getVal(pid),	e_sihih->getVal(pid), e_sipip->getVal(pid), e_r9->getVal(pid),
					e_sce->getVal(pid), e_sceta->getVal(pid), e_scphi->getVal(pid), e_e2x5max->getVal(pid), e_e1x5->getVal(pid),
					e_e5x5->getVal(pid), e_h2te->getVal(pid), e_h2tebc->getVal(pid), e_ooemoop->getVal(pid), e_fbrem->getVal(pid),
					e_eopin->getVal(pid), e_dEtaCalo->getVal(pid), e_kfchi2->getVal(pid), e_kfhits->getVal(pid), e_etawidth->getVal(pid), e_phiwidth->getVal(pid),
					e_e1x5e5x5->getVal(pid), e_preShowerOverRaw->getVal(pid), e_eopout->getVal(pid));
			electrons.push_back(tmp);
		}
	}
	return electrons;
}

vector<Jet> selectJetsCMG(const Event & ev, const JetVariables & jetVars, double ptMin, double etaMax) {

	const SingleVariableContainer<int> * jn = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(jetVars.jn));
	const ArrayVariableContainer<float> * j_px = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_px));
	const ArrayVariableContainer<float> * j_py = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_py));
	const ArrayVariableContainer<float> * j_pz = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_pz));
	const ArrayVariableContainer<float> * j_en = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_en));
	const ArrayVariableContainer<float> * j_btag = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_btag));
	const ArrayVariableContainer<float> * j_genpt = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(jetVars.j_genpt));

	vector<Jet> jets;
	for ( int i = 0; i < jn->getVal(); ++i ) {
		TLorentzVector jet(j_px->getVal(i), j_py->getVal(i), j_pz->getVal(i), j_en->getVal(i));
		if ( jet.Pt() > ptMin && fabs(jet.Eta()) < etaMax )
			jets.push_back( Jet(j_px->getVal(i), j_py->getVal(i), j_pz->getVal(i), j_en->getVal(i),
						j_btag->getVal(i), j_genpt->getVal(i)));
	}
	return jets;
}

TLorentzVector smearJets(vector<Jet> & jets) {
	TLorentzVector jetDiff;
	for ( unsigned j = 0; j < jets.size(); ++j ) {
		Jet smJet = smearedJet( jets[j] );
		jetDiff += (smJet.lorentzVector() - jets[j].lorentzVector());
		jets[j] = smJet;
	}
	return jetDiff;
}

Jet smearedJet(const Jet & origJet) {
	if (origJet.genpt <= 0)
		return origJet;

	//smearing factors are described in https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
	double eta = fabs(origJet.lorentzVector().Eta());
	double pt = origJet.lorentzVector().Pt();
	double ptSF = 1.0;
	double ptSF_err = 0.06;
	if (eta < 0.5) {
		ptSF = 1.066;
		ptSF_err = sqrt(pow(0.007, 2) + pow(0.5 * (0.07 + 0.072), 2));
	} else if (eta >= 0.5 && eta < 1.7) {
		ptSF = 1.191;
		ptSF_err = sqrt(pow(0.019, 2) + pow(0.5 * (0.06 + 0.062), 2));
	} else if (eta >= 1.7 && eta < 2.3) {
		ptSF = 1.096;
		ptSF_err = sqrt(pow(0.030, 2) + pow(0.5 * (0.08 + 0.085),2));
	} else if (eta >= 2.3 && eta < 5.0) {
		ptSF = 1.166;
		ptSF_err = sqrt(pow(0.050, 2) + pow(0.5 * (0.19 + 0.199), 2));
	}
	
	ptSF = max(0., (origJet.genpt + gRandom->Gaus(ptSF, ptSF_err) * (pt - origJet.genpt))) / pt;  //deterministic version
	if (ptSF <= 0)
		return origJet;

	double px = origJet.px * ptSF;
	double py = origJet.py * ptSF;
	double pz = origJet.pz;
	double mass = origJet.lorentzVector().M();
	double en = sqrt(mass * mass + px * px + py * py + pz * pz);

	Jet smearedJet(px, py, pz, en, origJet.btag, origJet.genpt);
	return smearedJet;
}


vector<Photon> selectPhotonsCMG(const Event & ev, const PhotonVariables & photonVars) {

	const SingleVariableContainer<int> * gn = dynamic_cast<const SingleVariableContainer<int> *>(ev.getVariable(photonVars.gn));
	const ArrayVariableContainer<float> * g_px = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_px));
	const ArrayVariableContainer<float> * g_py = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_py));
	const ArrayVariableContainer<float> * g_pz = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_pz));
	const ArrayVariableContainer<float> * g_en = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_en));
	const ArrayVariableContainer<float> * g_iso1 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_iso1));
	const ArrayVariableContainer<float> * g_iso2 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_iso2));
	const ArrayVariableContainer<float> * g_iso3 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_iso3));
	const ArrayVariableContainer<float> * g_sihih = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_sihih));
	const ArrayVariableContainer<float> * g_sipip = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_sipip));
	const ArrayVariableContainer<float> * g_r9 = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_r9));
	const ArrayVariableContainer<float> * g_hoe = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_hoe));
	const ArrayVariableContainer<float> * g_htoe = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_htoe));
	const ArrayVariableContainer<float> * g_corren = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_corren));
	const ArrayVariableContainer<float> * g_correnerr = dynamic_cast<const ArrayVariableContainer<float> *>(ev.getVariable(photonVars.g_correnerr));
	const ArrayVariableContainer<int> * g_idbits = dynamic_cast<const ArrayVariableContainer<int> *>(ev.getVariable(photonVars.g_idbits));

	vector<Photon> photons;
	for ( int i = 0; i < gn->getVal(); ++i ) {
		Photon tmpPhoton(g_px->getVal(i), g_py->getVal(i), g_pz->getVal(i), g_en->getVal(i), g_iso1->getVal(i), g_iso2->getVal(i), g_iso3->getVal(i), g_sihih->getVal(i),
				g_sipip->getVal(i), g_r9->getVal(i), g_hoe->getVal(i), g_htoe->getVal(i), g_corren->getVal(i), g_correnerr->getVal(i), g_idbits->getVal(i));
		photons.push_back(tmpPhoton);
	}
	return photons;
}

/*
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
