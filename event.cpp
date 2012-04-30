// ROOT Libraries
#include <TLeaf.h>
// Standard Libraries
#include <cstdlib>
#include <iostream>
#include <sstream>
// Other
#include "event.h"
#include "triggerinfo.h"

using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;
using std::map;

Event::Event( TTree *tree ) {
	TObjArray * treeLeaves = tree->GetListOfLeaves();
	int nEntries = treeLeaves->GetEntries();
	TObjArray::Iterator_t iter( treeLeaves );
	for( int i = 0; i < nEntries; ++i ) {
		iter.Next();
		TLeaf * leafPointer = dynamic_cast<TLeaf *>( *iter );
		//cout << leafPointer->GetTypeName() << endl;
		//cout << leafPointer->GetLen() << endl;
		//cout << leafPointer->GetBranch()->GetName() << endl;
		//cout << leafPointer->GetBranch()->GetClassName() << endl;

		string varName = leafPointer->GetBranch()->GetName();
		if (findVariable(varName)) {
			cout << "ERROR: Variable \"" << varName << "\" has already been	found!" << endl;
			exit( EXIT_FAILURE );
		}

//		cout << varName << endl;
		int leafLength = leafPointer->GetLen();
//		leafPointer->Print();
//		cout << leafLength << endl;
//		if (leafPointer->GetMaximum() > leafLength)
//			leafLength = leafPointer->GetMaximum();
//		cout << leafLength << endl;
		TLeaf * leafCount = leafPointer->GetLeafCount();
		if (leafCount && leafCount->GetMaximum() > leafLength) {
			leafLength = leafCount->GetMaximum();
//			leafCount->Print();
		}
//		cout << leafLength << endl;
//		std::cin.get();

		if ( leafPointer->GetBranch()->GetClassName() != string("")
				&& string(leafPointer->GetTypeName()) != string(leafPointer->GetBranch()->GetClassName()) )
			continue;
		if( leafPointer ) {
			if ( leafPointer->GetTypeName() == string("UInt_t") && leafLength == 1 ) {
				SingleVariableContainer<unsigned> * tempWrapper = new SingleVariableContainer<unsigned>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Int_t") && leafLength == 1 ) {
				SingleVariableContainer<int> * tempWrapper = new SingleVariableContainer<int>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Float_t") && leafLength == 1 ) {
				SingleVariableContainer<float> * tempWrapper = new SingleVariableContainer<float>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Double_t") && leafLength == 1 ) {
				SingleVariableContainer<double> * tempWrapper = new SingleVariableContainer<double>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Bool_t") && leafLength == 1 ) {
				SingleVariableContainer<bool> * tempWrapper = new SingleVariableContainer<bool>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("UInt_t") && leafLength > 1 ) {
				ArrayVariableContainer<unsigned> * tempWrapper = new ArrayVariableContainer<unsigned>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Int_t") && leafLength > 1 ) {
				ArrayVariableContainer<int> * tempWrapper = new ArrayVariableContainer<int>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Float_t") && leafLength > 1 ) {
				ArrayVariableContainer<float> * tempWrapper = new ArrayVariableContainer<float>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Double_t") && leafLength > 1 ) {
				ArrayVariableContainer<double> * tempWrapper = new ArrayVariableContainer<double>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Bool_t") && leafLength > 1 ) {
				ArrayVariableContainer<bool> * tempWrapper = new ArrayVariableContainer<bool>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("vector<int>") && leafLength == 1 ) {
				vector<int> ** tmpPtr = new vector<int> *;
				*tmpPtr = 0;
				tree->SetBranchAddress( leafPointer->GetBranch()->GetName(), tmpPtr );
				pair<string, vector<int> **> temp(leafPointer->GetBranch()->GetName(), tmpPtr );
				VectorsInt.push_back(temp);
			} else if ( leafPointer->GetTypeName() == string("vector<float>") && leafLength == 1 ) {
				vector<float> ** tmpPtr = new vector<float> *;
				*tmpPtr = 0;
				tree->SetBranchAddress( leafPointer->GetBranch()->GetName(), tmpPtr );
				pair<string, vector<float> **> temp(leafPointer->GetBranch()->GetName(), tmpPtr );
				VectorsFloat.push_back(temp);
			} else if ( leafPointer->GetTypeName() == string("vector<double>") && leafLength == 1 ) {
				vector<double> ** tmpPtr = new vector<double> *;
				*tmpPtr = 0;
				tree->SetBranchAddress( leafPointer->GetBranch()->GetName(), tmpPtr );
				pair<string, vector<double> **> temp(leafPointer->GetBranch()->GetName(), tmpPtr );
				VectorsDouble.push_back(temp);
			} else if ( leafPointer->GetTypeName() == string("TriggerInfo") && leafLength == 1 ) {
				TriggerInfo ** tmpPtr = new TriggerInfo *;
				*tmpPtr = 0;
				tree->SetBranchAddress( leafPointer->GetBranch()->GetName(), tmpPtr );
				pair<string, TriggerInfo **> temp(leafPointer->GetBranch()->GetName(), tmpPtr );
				Triggers.push_back(temp);
			} else {
				cout << "ERROR: Don't know what to do with this type: " << leafPointer->GetTypeName() << '[' << leafLength << "] " << varName << " !" << endl;
				//exit (EXIT_FAILURE);
			}
		}
	}
}

Event::~Event() {
	for (auto iter = variables.begin(); iter != variables.end(); ++iter)
		delete *iter;

	for (unsigned i = 0; i < VectorsInt.size(); ++i)
		delete (VectorsInt[i].second);
	for (unsigned i = 0; i < VectorsFloat.size(); ++i)
		delete (VectorsFloat[i].second);
	for (unsigned i = 0; i < VectorsDouble.size(); ++i)
		delete (VectorsDouble[i].second);

	for (unsigned i = 0; i < Triggers.size(); ++i)
		delete (Triggers[i].second);
}

//==================================================================================

VariableContainer * Event::findVariable(const std::string & name) const {
	auto location = variables.begin();
	for (; location != variables.end(); ++location) {
		if ( name == (*location)->getName() )
			break;
	}
	if (location != variables.end())
		return *location;
	return 0;
}

unsigned Event::findVariableIndex(const std::string & name) const {
	unsigned location = 0;
	for (; location < variables.size(); ++location) {
		if ( name == variables[location]->getName() )
			break;
	}
	if (location < variables.size())
		return location;
	return -1;
}
		
const VariableContainer * Event::getVariable(unsigned i) const {
	if (i < variables.size())
		return variables[i];
	return 0;
}

//==================================================================================

const vector<int> * Event::getVectorIntAdr(const string & name) const {
	for (unsigned i = 0; i < VectorsInt.size(); ++i) {
		if ( VectorsInt[i].first == name )
			return *(VectorsInt[i].second);
	}
	cout << "ERROR: Couldn't find " << name << " variable!" << endl;
	exit( EXIT_FAILURE );
}

const vector<float> * Event::getVectorFloatAdr(const string & name) const {
	for (unsigned i = 0; i < VectorsFloat.size(); ++i) {
		if ( VectorsFloat[i].first == name )
			return *(VectorsFloat[i].second);
	}
	cout << "ERROR: Couldn't find " << name << " variable!" << endl;
	exit( EXIT_FAILURE );
}

const vector<double> * Event::getVectorDoubleAdr(const string & name) const {
	for (unsigned i = 0; i < VectorsDouble.size(); ++i) {
		if ( VectorsDouble[i].first == name )
			return *(VectorsDouble[i].second);
	}
	cout << "ERROR: Couldn't find " << name << " variable!" << endl;
	exit( EXIT_FAILURE );
}

//==================================================================================

const TriggerInfo * Event::getTriggerInfo(const string & name) const {
	for (unsigned i = 0; i < Triggers.size(); ++i) {
		if ( Triggers[i].first == name )
			return *(Triggers[i].second);
	}
	cout << "ERROR: Couldn't find " << name << " variable!" << endl;
	exit( EXIT_FAILURE );
}

