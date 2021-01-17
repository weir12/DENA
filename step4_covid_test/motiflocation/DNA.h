// DNA.h: interface for the DNA class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DNA_H__99F6EF20_97C1_4A64_B705_D286EBA19B12__INCLUDED_)
#define AFX_DNA_H__99F6EF20_97C1_4A64_B705_D286EBA19B12__INCLUDED_

#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>      // old <ctype.h>
#include <string>
#include "motif.h"
#include "gzstream.h"

using namespace std ;

class DNA  
{
public:
	DNA();
	DNA(char *seq);
	virtual ~DNA();

	vector<string> sequence;
	vector<string> header;
	vector<unsigned int> SeqLength;
	int NSeqs;

	bool readFASTA(const char* filename, bool uppercase_only, bool lowercase_only, bool KG_ref);
	bool readFASTA(const char* filename, const string &header_name, bool uppercase_only, bool lowercase_only, bool KG_ref);

	int CountMotifOccurances(motif &Motif);

	struct ToLower
	{
		char operator() (char c) const  { return tolower(c); }
	};

	struct ToUpper
	{
		char operator() (char c) const  { return toupper(c); }
	};

	private:
		void Trim(string &str);
		void remove_upper_or_lowercase(string &string1, bool uppercase_only, bool lowercase_only);

};

#endif // !defined(AFX_DNA_H__99F6EF20_97C1_4A64_B705_D286EBA19B12__INCLUDED_)
