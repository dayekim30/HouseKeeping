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
typedef vector<string> mutation;

class Parsing {

public:
	Parsing(int any) {
		kmerFromgene = new unordered_map<string, vector<string>>();
		mutFromkmer = new unordered_map<string, vector<mutation>>();
		k = any;
	}
	~Parsing() {}

public:
	void Aread(const string& filename);
	void Qread(const string& filename);
	void Nucio();

public:
	unordered_map<string, vector<string>>* kmerFromgene;
	unordered_map<string, vector<mutation>>* mutFromkmer;
	unordered_map<string, char>* nuciomap;
	int k;
};
