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

		int leafLength = leafPointer->GetLen();
		TLeaf * leafCount = leafPointer->GetLeafCount();
		if (leafCount && leafCount->GetMaximum() > leafLength) {
			leafLength = leafCount->GetMaximum();
		}
		bool array = (leafCount || leafLength > 1);

		if ( leafPointer->GetBranch()->GetClassName() != string("")
				&& string(leafPointer->GetTypeName()) != string(leafPointer->GetBranch()->GetClassName()) )
			continue;
		if( leafPointer ) {
			if ( leafPointer->GetTypeName() == string("UInt_t") && !array ) {
				SingleVariableContainer<unsigned> * tempWrapper = new SingleVariableContainer<unsigned>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Int_t") && !array ) {
				SingleVariableContainer<int> * tempWrapper = new SingleVariableContainer<int>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Float_t") && !array ) {
				SingleVariableContainer<float> * tempWrapper = new SingleVariableContainer<float>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Double_t") && !array ) {
				SingleVariableContainer<double> * tempWrapper = new SingleVariableContainer<double>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Bool_t") && !array ) {
				SingleVariableContainer<bool> * tempWrapper = new SingleVariableContainer<bool>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("UInt_t") && array ) {
				ArrayVariableContainer<unsigned> * tempWrapper = new ArrayVariableContainer<unsigned>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Int_t") && array ) {
				ArrayVariableContainer<int> * tempWrapper = new ArrayVariableContainer<int>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Float_t") && array ) {
				ArrayVariableContainer<float> * tempWrapper = new ArrayVariableContainer<float>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Double_t") && array ) {
				ArrayVariableContainer<double> * tempWrapper = new ArrayVariableContainer<double>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("Bool_t") && array ) {
				ArrayVariableContainer<bool> * tempWrapper = new ArrayVariableContainer<bool>(varName, leafLength);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("vector<int>") && !array ) {
				VectorVariableContainer<int> * tempWrapper = new VectorVariableContainer<int>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtrToPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("vector<float>") && !array ) {
				VectorVariableContainer<float> * tempWrapper = new VectorVariableContainer<float>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtrToPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("vector<double>") && !array ) {
				VectorVariableContainer<double> * tempWrapper = new VectorVariableContainer<double>(varName);
				tree->SetBranchAddress( varName.c_str(), tempWrapper->getPtrToPtr() );
				variables.push_back(tempWrapper);
			} else if ( leafPointer->GetTypeName() == string("TriggerInfo") && !array ) {
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
	else
		throw string("ERROR: Event::getVariable(unsigned i) : Index out of bounds!");
	return 0;
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

