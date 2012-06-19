#ifndef VARIABLEGETTER_H
#define VARIABLEGETTER_H

class Event;
class VariableContainer;

class VariableGetter {
	public:
		enum VGTYPE { UNS, INT, DBL, VEI, VEF };
		VariableGetter(std::string name, const Event & ev);
		~VariableGetter() {};
		double getValue() const;
		double getValue(unsigned i) const;
		std::string getName() const {
			return name;
		}
		VGTYPE getType() const {
			return type;
		}
		int size() const;
	private:
		std::string name;
		VariableContainer * var;
		VGTYPE type;
};

#endif // VARIABLEGETTER_H
