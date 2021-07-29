#pragma once
#include <string>

using namespace std;
class Sequence {

public:
	Sequence(const string input) {
		seq = input;
	}
	~Sequence() {}

public: 
	string forwardSeq();
	string backwardSeq();

public:
	string seq;
};