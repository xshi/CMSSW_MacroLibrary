#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <map>

class Options {
	private:
		std::map<std::string, bool> boolOptions;
		std::map<std::string, std::string> stringOptions;
	public:
		Options();
		void readInOptions(const std::string & fileName);
		bool checkBoolOption(const std::string & name) const;
		const std::string & checkStringOption(const std::string & name) const;
		void addBoolOption(const std::string & name, bool value);
		void addStringOption( const std::string & name, const std::string & value );
};

#endif
