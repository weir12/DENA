#ifndef __MOTIF_H_
#define __MOTIF_H_

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include "random.h"

using namespace std;

class motif
{
public:
	motif();
	motif(bool is_nondegenerate);
	motif(string input);
	virtual ~motif();
	string sequence;
	string complement;
	int count1;
	int count2;
	double pvalue;
	double score;	// aka fitness
	double degeneracy;
	double relative_ratio;
	unsigned int max_length;		// Maximum possible length of motif.

	double ExpVal;	// Expected number of offspring per generation

	void CalcDegeneracy();
	bool IsMotifInSequence(const char* inputseq, const char* wild, int allowed_mutations);
	vector<int> returnMotifLocations(const char* inputseq, const char* wild, int allowed_mutations);
	void get_max_length();
	inline bool basematch(const char base1, const char base2);
	char GetRandomBase(bool is_nondegenerate);
	void doComplement();
	char getComplement(const char base);
};

#endif

