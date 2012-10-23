// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
// ROOT Libraries
#include <TLorentzVector.h>
// Other
#include "cut.h"
#include "event.h"

using std::setw;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;

string DCut::GetName() const {
	stringstream temp;
	temp << varName_;
	if ( cT_ == GT ) {
		temp << "_GT_";
		temp << GetCutValue();
	} else if ( cT_ == LT ) {
		temp << "_LT_";
		temp << GetCutValue();
	} else if ( cT_ == EQ ) {
		temp << "_EQ_";
		temp << GetCutValue();
	} else if ( cT_ == WIN ) {
		temp << "_WIN_";
		temp << GetCutValue();
		temp << "_";
		temp << GetCutValue1();
	}
	return temp.str();
}

bool Cut::operator()() const {
	if ( desc_.GetType() == DCut::GT )
		return ( var_->getValue() > desc_.GetCutValue() );
	else if ( desc_.GetType() == DCut::LT )
		return ( var_->getValue() < desc_.GetCutValue() );
	else if ( desc_.GetType() == DCut::EQ )
		return ( var_->getValue() == desc_.GetCutValue() );
	else if ( desc_.GetType() == DCut::WIN )
		return ( var_->getValue() < desc_.GetCutValue() || var_->getValue() > desc_.GetCutValue1() );
	else
		return true;
}

string Cut::GetName() const {
	return desc_.GetName();
}
