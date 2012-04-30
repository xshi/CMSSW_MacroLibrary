#ifndef EVENT_H
#define EVENT_H

// ROOT Libraries
#include <TTree.h>
// Standard Libraries
#include <string>
#include <utility>
#include <vector>
#include <map>


class VariableContainer {
	private:
		std::string name;
		VariableContainer(const VariableContainer &);
		VariableContainer & operator=(const VariableContainer &);
	public:
		VariableContainer(const std::string & n) : name(n) {};
		virtual ~VariableContainer() {};
		const std::string & getName() const {
			return name;
		};
};

template <typename T> class SingleVariableContainer : public VariableContainer {
	private:
		T * varPtr;
		SingleVariableContainer(const SingleVariableContainer &);
		SingleVariableContainer & operator=(const SingleVariableContainer &);
	public:
		SingleVariableContainer(const std::string & n) : VariableContainer(n) {
			varPtr = new T;
		}
		virtual ~SingleVariableContainer() {
			delete varPtr;
		}
		T getVal() const {
			return *varPtr;
		}
		T * getPtr() const {
			return varPtr;
		}
};

template <typename T> class ArrayVariableContainer : public VariableContainer {
	private:
		unsigned size;
		T * varPtr;
		ArrayVariableContainer(const ArrayVariableContainer &);
		ArrayVariableContainer & operator=(const ArrayVariableContainer &);
	public:
		ArrayVariableContainer(const std::string & nm, unsigned n ) : VariableContainer(nm), size(n) {
			varPtr = new T[size];
		}
		virtual ~ArrayVariableContainer() {
			delete[] varPtr;
		}
		T getVal(unsigned i) const {
			if (i >= size)
				throw "ERROR: ArrayVariableContainer::getVal : Index out of bounds!";
			return varPtr[i];
		}
		T * getPtr() const {
			return varPtr;
		}
};

class TriggerInfo;

class Event {
	private:
		std::vector<VariableContainer *> variables;

		std::vector<std::pair<std::string, std::vector<int> **> > VectorsInt;
		std::vector<std::pair<std::string, std::vector<float> **> > VectorsFloat;
		std::vector<std::pair<std::string, std::vector<double> **> > VectorsDouble;
		
		std::vector<std::pair<std::string, TriggerInfo **> > Triggers;

		Event(const Event &);
		Event & operator=(const Event &);
	public :
		Event(TTree *tree);
		~Event();
		VariableContainer * findVariable(const std::string & name) const;
		unsigned findVariableIndex(const std::string & name) const;

		template <typename T> T getSingleVariableValue(const std::string & name) const {
			VariableContainer * tempPtr = findVariable(name);
			if (tempPtr) {
				SingleVariableContainer<T> * varPtr = dynamic_cast<SingleVariableContainer<T> *>(tempPtr);
				if (varPtr)
					return varPtr->getVal();
				else
					throw "ERROR: Event::getSingleVariableValue : Variable name (" + name + ") does not match expected type!";
			} else
				throw "ERROR: Event::getSingleVariableValue : Variable (" + name + ") can't be found!";
		}
		
		template <typename T> inline T getSVV(const std::string & name) const {
			return getSingleVariableValue<T>(name);
		}

		template <typename T> T * getSingleVariableAddress(const std::string & name) const {
			VariableContainer * tempPtr = findVariable(name);
			if (tempPtr) {
				SingleVariableContainer<T> * varPtr = dynamic_cast<SingleVariableContainer<T> *>(tempPtr);
				if (varPtr)
					return varPtr->getPtr();
			}
			return 0;
		}
		
		template <typename T> inline T * getSVA(const std::string & name) const {
			return getSingleVariableAddress<T>(name);
		}

		template <typename T> T getArrayVariableValue(const std::string & name, unsigned i) const {
			VariableContainer * tempPtr = findVariable(name);
			if (tempPtr) {
				ArrayVariableContainer<T> * varPtr = dynamic_cast<ArrayVariableContainer<T> *>(tempPtr);
				if (varPtr)
					return varPtr->getVal(i);
				else
					throw "ERROR: Event::getArrayVariableValue : Variable name (" + name + ") does not match expected type!";
			} else
				throw "ERROR: Event::getArrayVariableValue : Variable (" + name + ") can't be found!";
		}

		template <typename T> inline T getAVV(const std::string & name, unsigned i) const {
			return getArrayVariableValue<T>(name, i);
		}

		template <typename T> T * getArrayVariableAddress(const std::string & name) const {
			VariableContainer * tempPtr = findVariable(name);
			if (tempPtr) {
				ArrayVariableContainer<T> * varPtr = dynamic_cast<ArrayVariableContainer<T> *>(tempPtr);
				if (varPtr)
					return varPtr->getPtr();
			}
			return 0;
		}

		template <typename T> inline T * getAVA(const std::string & name) const {
			return getArrayVariableAddress<T>(name);
		}

		const std::vector<int> * getVectorIntAdr(const std::string & name) const;
		const std::vector<float> * getVectorFloatAdr(const std::string & name) const;
		const std::vector<double> * getVectorDoubleAdr(const std::string & name) const;
		
		const TriggerInfo * getTriggerInfo(const std::string & name) const;
};

#endif // EVENT_H
