// Other
#include "variableGetter.h"

#ifndef CUT_H
#define CUT_H

class Event;
class Samples;

class DCut {
	public :
		enum TYPE { GT, LT, EQ, WIN, NO };
	private:
		std::string varName_;
		DCut::TYPE cT_;
		double cutValue_;
		double cutValue1_;
	public:
		DCut( std::string variable = "", DCut::TYPE cType = NO, double val = 0, double val1 = 0 )
			: varName_(variable), cT_(cType), cutValue_(val), cutValue1_(val1) {};
		virtual ~DCut() {};
		virtual std::string GetName() const;
		virtual std::string GetVarName() const {
			return varName_;
		};
		virtual DCut::TYPE GetType() const {
			return cT_;
		}
		virtual double GetCutValue() const {
			return cutValue_;
		};
		virtual double GetCutValue1() const {
			return cutValue1_;
		};
		virtual void SetCutValue(double val) {
			cutValue_ = val;
		}
		virtual void SetCutValue1(double val) {
			cutValue1_ = val;
		}
};

class Cut {
	private:
		DCut desc_;
		VariableGetter * var_;
	public:
		Cut (const Event & ev, const DCut & cut) : desc_(cut) {
			if ( cut.GetVarName() == "" && cut.GetType() == DCut::NO )
				var_ = new VariableGetter("Event", ev);
			else
				var_ = new VariableGetter(cut.GetVarName(), ev);
		}
		Cut (const Cut & cut) : desc_(cut.desc_) {
			var_ = new VariableGetter( *cut.var_ );
		}
		virtual ~Cut() {
			delete var_;
		}
		virtual bool operator()() const;
		virtual std::string GetName() const;
		virtual std::string GetVarName() const {
			return desc_.GetVarName();
		};
		virtual double GetCutValue() const {
			return desc_.GetCutValue();
		};
		virtual void SetCutValue (double val) {
			desc_.SetCutValue(val);
		}
//		virtual double Optimize(const Samples & samples, OptConf & conf);
};

#endif // CUT_H
