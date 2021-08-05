#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <unordered_set>

#include "Sequence.h"
#include "Frames.h"

using namespace std;
typedef vector<string> *heads;
typedef unordered_set<string>* sset;
typedef vector<sset> *mutation;

class Parsing {

public:
	Parsing(const int &any) {
		kmerFromgene = new unordered_map<string, heads>();
		mutFromkmer = new unordered_map<string, mutation>();
		k = any;
	}
	~Parsing() {
		delete kmerFromgene, mutFromkmer, nuciomap;
	}

public:
	void Aread(const string& filename);
	void Qread(const string& filename);
	void Nucio();

public:
	unordered_map<string, heads>* kmerFromgene;
	unordered_map<string, mutation>* mutFromkmer;
	unordered_map<string, char>* nuciomap;
	int k;
};
