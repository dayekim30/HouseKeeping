#include "Parsing.h"
using namespace std;



void Parsing::Aread(const string &filename)
{
	
	fstream newfile;
	string id;
	string seq;

	newfile.open(filename, ios::in);
	cout << "Reading Gene database, creating k - mer mapping (k = " << k << ")" << endl;
	
	if (newfile.is_open()) {
		string tp;
		int i = 0;
		//int l = 0;
		getline(newfile, tp);

		while (!newfile.eof()) {
			
			if (tp[0] != '>') { cout << "wrong FASTA format" << endl; exit(0); }
			id = tp;
			seq = "";
			
			
			do {
				
				getline(newfile, tp);
				//l++;

				if (tp[0] == '>'|| newfile.eof())break;
				seq = seq + tp;
				
			} while (tp[0] != '>');
			
			i++;
			//cout << "l:" << l<< " ";
			//cout << "whats wrong with you length? " << seq.length() << "\n";
			//id and seq put in the hash
			if (seq.size() < k)continue;
			for (int j = 0; j < seq.size() - k + 1; j++) {

				string kmer = seq.substr(j, k);
				
				if (!(kmerFromgene->count(kmer))) {
					
					vector<string> ai;
					ai.push_back(id);
					
					kmerFromgene->insert(make_pair(kmer, ai));
				}
				else {
					vector<string> &ai = kmerFromgene->at(kmer);
					ai.push_back(id);
					
					kmerFromgene->insert_or_assign(kmer, ai);
				}

			}
			

			// triming -> id spliting |,; 
			size_t found = id.find('|');
			short first = found;
			string tamp = "";
			if (found != string::npos) {
				found = id.find('|', found + 1);
				if (found != string::npos) {
					tamp = id.substr(first+1, found - first-1);
				}
			}
			
			//cout << "tamp: " << tamp << endl;
			vector <string> tokens;
			
			while (tamp.find(';') != string::npos) {
				int com = tamp.find(';');
				tokens.push_back(tamp.substr(0, com));
				tamp = tamp.substr(com+1);
				//cout << "the editted tamp: " << tamp << endl;
			}
			if (tamp.find(';') == string::npos) tokens.push_back(tamp);

			//getting mutation hash table from tokens
			
			vector<mutation>* muts = new vector<mutation>();
			mutFromkmer->insert(make_pair(id, *muts));

			for (string a : tokens) {
				int n;
				int st;
				char c;
				for (int h = 0; h < a.length(); h++) {
					if (isdigit(a.at(h))) { st = h; break; }
				}


				size_t found = a.find("STOP");
				if (found != string::npos) {
					n = stoi(a.substr(st, found - 1));
					c = '*';
				}
				else {
					n = stoi(a.substr(st, a.length() - 2));
					c = a.at(a.length() - 1);
				}
				 
			/*	cout << "the target index: " << n << endl;
			
				cout << "the target character: " << c << endl;
				
				cout << "real index n-1 : " << n - 1 << endl;

				
				cout << "whats wrong with you length? " << seq.length() << "\n";
				cout << "seq.at(n-1): " << seq.at(n - 1) << "\n";*/
				
				// check there is targeted mutation in this sequence
				if (seq.at(n-1) == c) {
					//cout << "yes it is here" << endl;
					
					//make mutaion hashtable
					mutation mt;
					for (int x = 1; x < k + 1; x++) {
						int ft = (n - 1) - (k - x);
						if (ft < 0 || n+x-1>seq.length()) break;
						mt.push_back(seq.substr(ft, k));
					}
					muts->push_back(mt);
					mutFromkmer->insert_or_assign(id, *muts);				
				}
				

				
				

			}
			if (i % 1000 == 0) {
				cout << i << ".. ";
			}
			
		}
		cout << "\nthe number of genes: " << i <<"\n";
		newfile.close();
	}
		
}

void Parsing::Nucio()
{
	string filename = "neclo.txt";
	fstream newfile;
	newfile.open(filename, ios::in); //open a file to perform read operation using file object

	nuciomap = new unordered_map<string, char>();



	if (newfile.is_open()) {   //checking whether the file is open
		string tp;

		int i = 0;
		while ((getline(newfile, tp)) && i < 64) {
			vector<string> result;
			stringstream ss(tp);

			while (ss.good()) {
				string substr;
				getline(ss, substr, ',');
				result.push_back(substr);
			}

			nuciomap->insert(make_pair(result[0], result[1][0]));
			i++;
		}
		//read data from file object and put it into string. 
		//cout << "the number of set is " << i << "and" << nuciomap->size() << endl;
	}
	newfile.close(); //close the file object.
	

}

bool sortByVal(const pair<string, int>& a,
	const pair<string, int>& b)
{
	return (a.second < b.second);
}


void Parsing::Qread(const string &filename)
{
	fstream infile;
	string id;
	string seq;
	Nucio();
	infile.open(filename, ios::in);
	cout << "Translation starts (k = " << k << ")" << endl;

	if (infile.is_open()) {
		
		fstream fout;
		fout.open("AntibioticsResistence.csv", ios::out | ios::app);
		fout << "ids, matching_max_number ,Best_Matching_Seq" << "\n";
		
		string tp;
		int i = 0;
		
		while (getline(infile, tp)) {
			if (tp[0] != '@') {
				cout << "this is wrong FASTQ file format." << endl;
				exit(0);		
			}
			id = tp;
			getline(infile, tp);
			seq = tp;
			getline(infile, tp);
			getline(infile, tp);
			frames fr = frames(*nuciomap);
			fr.getAllframes(seq);
			i++;
			// max hit search
			int tmax = 0;
			string tmagen = "";
			for (int f = 0; f < 6; f++) {
				unordered_map<string, float> weight;
				unordered_set<string> fkmers;
				string fseq = fr.result[f];

				//cout << "fseq -frame(" << f << ") :" << fseq << endl;
				for (int i = 0; i < fseq.size() - k + 1; i++) {

					string kmer = fseq.substr(i, k);
					fkmers.insert(kmer);
					if (kmerFromgene->count(kmer)) {

						auto &flist = kmerFromgene->at(kmer);
						//cout << "kmerFromgene hashmap size: " << flist.size() << endl;
						for (string fid : flist) {
							float fwg = (float)(1 / (float)flist.size());

							if (!weight.count(fid)) {
								weight.insert(make_pair(fid, fwg));
								//cout << "float number 1: " << fwg << endl;
							}
							else {
								weight.insert_or_assign(fid, weight.at(fid) + fwg);
								//cout << "float number 2: " << fwg <<"float sum: "<< weight.at(fid) << endl;
							}
						
						}
												
					}

				}
				vector<pair<string, float>> vec;
				for (auto at : weight) {
					vec.push_back(make_pair(at.first, at.second));
				}
				sort(vec.begin(), vec.end(), sortByVal);
				//cout << "vec size(Hash set size): " << vec.size() << endl;

				float cmax = 0;
				string cmid = "";
				if (vec.empty()) {
					//cout << "nothing mathched" << endl;
					continue;
				
				}
				else {
					cmax = vec.back().second;
					cmid = vec.back().first;

					//cout << "cmax: " << (float)cmax << endl;;
					//cout << "cmid: " << cmid<<endl;

				}
				//check it has mutation in all lists.
				
				auto &mutlist = mutFromkmer->at(cmid);
				bool flag = false;
				for (auto muts : mutlist) {
					flag = false;
					for (auto st : muts) {
						if (fkmers.count(st)) flag = true;
					}

					if (!flag) {
						//cout << "there is no mut matching - ignore" << endl;
						break;
					}
				}
		 
				if (flag) {
					/*tmax = cmax;
					tmagen = cmid;*/
					//cout << "this is ture!" << endl;
					// file write
					if (tmax < cmax) {
						tmax = cmax;
						tmagen = cmid;
					}				
				}
			}

			//file write : tamx and cmax
			if (tmax>0) {
			//write ???- no matching
				//cout << "here is assigned one - tmagen: " << tmagen << endl;
				//cout << "tmax: " << tmax << endl;
				fout << id << "," << tmax << "," << tmagen << "\n";
			}
			else {
			
				fout << id << ",? , ?"<<"\n";
			}
			if ((i % 10000) == 0) cout << i << ".. ";
			
		}
		cout << "\nthe number of sequences: " << i << "\n";
		fout.close();
		infile.close();
	}
	
}

