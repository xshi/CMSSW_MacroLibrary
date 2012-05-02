// ROOT Libraries
#include <TLorentzVector.h>
// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
// Other
#include "event.h"
#include "variableGetter.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;

VariableGetter::VariableGetter(string hName, const Event & ev) {
	if ( (valueU_ = ev.getSingleVariableAddress<unsigned>(hName)) ) {
		type_ = UNSIGNED;	
	} else if ( (valueI_ = ev.getSingleVariableAddress<int>(hName)) ) {
		type_ = INT;
	} else if ( (valueD_ = ev.getSingleVariableAddress<double>(hName)) ) {
		type_ = DOUBLE;
	} else
		throw string("ERROR: VariableGetter::VariableGetter : Don't know what to do with this variable: " + hName + "!");
}

double VariableGetter::getValue() const {
	if (type_ == UNSIGNED)
		return (double) (*valueU_);
	else if (type_ == INT)
		return (double) (*valueI_);
	else if (type_ == DOUBLE)
		return (*valueD_);
	else
		return -999;
}
