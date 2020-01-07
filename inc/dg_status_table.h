#ifndef DG_STATUS_TABLE
#define DG_STATUS_TABLE

#include <unordered_map>

class Status{

public:
	static std::unordered_map<char, char[4]> status_lookup;


	// constructor
	Status();

};

#endif
