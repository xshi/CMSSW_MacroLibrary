#ifndef VARIABLEGETTER_H
#define VARIABLEGETTER_H

class Event;

class VariableGetter {
	public:
		VariableGetter(std::string name, const Event & ev);
		~VariableGetter() {};
		double getValue() const;
	private:
		enum TYPE { UNSIGNED, INT, DOUBLE };
		TYPE type_;
		const unsigned * valueU_;
		const int * valueI_;
		const double * valueD_;
};

#endif // VARIABLEGETTER_H
