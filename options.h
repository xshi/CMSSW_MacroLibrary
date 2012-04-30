#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <vector>

class Options {
	private:
		std::vector<std::string> boolOptions;
		std::vector<std::pair<std::string, std::string> > stringOptions;
	public:
		Options();
		void readOptions(int argc, const char * argv[]);
		bool checkBoolOption(const std::string & name) const;
		const std::string & checkStringOption(const std::string & name) const;
};

#endif
