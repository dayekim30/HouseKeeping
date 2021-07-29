#pragma once
#include <string>
#include <unordered_map>
#include "Sequence.h"

using namespace std;
typedef unordered_map<string, char> nuclo;
class frames {
public:
	frames() {}
	~frames() {}
	frames(nuclo& hmap);
	
	string setPro(string seq);
	string StrOrderchange(int a, string seq);
	void getAllframes(string seq);
	void setFrames(string g[6]);
private:
	nuclo map;
public:
	string result[6];

};