#include "motif.h"


motif::motif()
{
}

motif::motif(bool is_nondegenerate)
{
	// Generate a random motif.
	int i;
    int maxlength=6;
	int length;
	string out;
	ExpVal = 0.0;

	length = ran2(2, maxlength+1);
	sequence.resize(length);
	for (i=0; i<length; i++)
	{
		sequence[i] = GetRandomBase(is_nondegenerate);
	}
	doComplement();

	CalcDegeneracy();
	get_max_length();
}

motif::motif(string input)
{
	// Set motif to input sting.
	ExpVal = 0.0;
	sequence = input;
	CalcDegeneracy();
	get_max_length();
	doComplement();
}

motif::~motif()
{
}

char motif::getComplement(char base)
{
	switch (base)
	{
	case 'a': return 't';
	case 'c': return 'g';
	case 'g': return 'c';
	case 't': return 'a';
	case 'r': return 'y';
	case 'y': return 'r';
	case 'k': return 'm';
	case 'm': return 'k';
	case 's': return 's';
	case 'w': return 'w';
	case 'b': return 'v';
	case 'd': return 'h';
	case 'h': return 'd';
	case 'v': return 'b';
	case 'n': return 'n';
	case '1': return '1';
	case '2': return '2';
	case '3': return '3';
	case '4': return '4';
	case '5': return '5';
	case '6': return '6';
	case '7': return '7';
	case '8': return '8';
	case '9': return '9';
	default:
		cout << "Cannot find base complement of : " << base << endl;
		abort();
		return base;
	}
}


void motif::doComplement()
{
	unsigned int i;
	char base;
	unsigned int seqlen;
	seqlen = (unsigned int)sequence.length();
	complement.resize(seqlen);
	for (i=0; i<seqlen; i++)
	{
		base = sequence[i];
		complement[seqlen-i-1] = getComplement(base);
	}
}

inline bool motif::basematch(const char base1, const char base2)
{
	switch (base1)
	{
		case 'a':
			switch (base2)
			{
			/*case 'c': case 'g': case 't': case 'y': case 'k': case 's': case 'b':	// quicker to test the -ve case
				return false;
			default:
				return true;*/
			case 'a': case 'r': case 'm': case 'w': case 'd': case 'h': case 'v': case 'n':
				return true; 
			default:
			
				return false; 
			}
			break;
		case 'c': 
			switch (base2)
			{
			/*case 'a': case 'g': case 't': case 'r': case 'k': case 'w': case 'h':	// quicker to test the -ve case
				return false;
			default:
				return true;*/
			case 'c': case 'y': case 'm': case 's': case 'b': case 'h': case 'v': case 'n':
				return true; 
			default:
				return false; 
			}
			break;
		case 'g':
			switch (base2)
			{
			/*case 'a': case 'c': case 't': case 'y': case 'm': case 'w': case 'h':	// quicker to test the -ve case
				return false;
			default:
				return true;*/
			case 'g': case 'r': case 'k': case 's': case 'b': case 'd': case 'v': case 'n':
				return true; 
			default:
				return false; 
			}
			break;
		case 't':
			switch (base2)
			{
			/*case 'a': case 'c': case 'g': case 'r': case 'm': case 's': case 'v':	// quicker to test the -ve case
				return false;
			default:
				return true;*/
			case 't': case 'y': case 'k': case 'w': case 'b': case 'd': case 'h': case 'n':
				return true; 
			default:
				return false; 
			}
			break;
		case 'n': 
			return false; 
		default:
			cout << "\nNon-standard character in sequence file: " << base1 << " : " << base2 << endl;
			//abort();
			return false;
		}
	return false;
}


// The main search routine. This code is running 99.9% of the time so needs to be fast
// (hence the use of char*)
bool motif::IsMotifInSequence(const char* inputseq, const char* wild, int allowed_mutations)
{
	const char *cp = inputseq;
	const char *mp = wild;
	bool match;
	int max_gap_size;
	int max_gap_size_plus_one;
	char *buffer;
	size_t sizeof_char = sizeof( char ); 
	char null_char = (char)NULL;
	int mutations = 0;

	while ((*inputseq) && (*wild))
	{
		if (basematch(*inputseq, *wild))	// Compare degenerate bases
		{
			wild++;		// Move to next character
			inputseq++;
		} 
		else if ((*wild > '0') && (*wild <= '9'))	// This bit allows for gaps in the motif.
		{
			max_gap_size = (*wild - '0') + max_length;
			max_gap_size_plus_one = max_gap_size+1;
			while ((*wild) && (*wild > '0') && (*wild <= '9'))		// Remove all prefixed numbers
				wild++;
			buffer = new char[max_gap_size_plus_one];
			strncpy(buffer, inputseq, sizeof_char * max_gap_size_plus_one);
			buffer[max_gap_size] = null_char;
			match = IsMotifInSequence(buffer, wild, allowed_mutations - mutations);
			delete [] buffer;
			if (match)
				return true;	// The call to the recursion returned true, so we should return true.
			else
			{
				wild = mp;	// We do not have a match,
				cp++;		// so reset the wild string, and goto the next char on the inputsequence.
				inputseq = cp;
				mutations = 0;
			}
		}
		else if (mutations < allowed_mutations)
		{
			wild++;		// Move to next character
			inputseq++;
			mutations++;
		}
		else 
		{
			wild = mp;	// We do not have a match,
			cp++;		// so reset the wild string, and goto the next char on the inputsequence.
			inputseq = cp;
			mutations = 0;
		}
	}

	return (bool)!*wild;
}


vector<int> motif::returnMotifLocations(const char* inputseq, const char* wild, int allowed_mutations)
{
	const char *cp = inputseq;
	const char *mp = wild;
	bool match;
	int max_gap_size;
	vector<int> positions;
	int pos = 0;
	int max_gap_size_plus_one;
	char *buffer;
	size_t sizeof_char = sizeof( char ); 
	char null_char = (char)NULL;
	int mutations = 0;

	while (*inputseq)
	{
		if (basematch(*inputseq, *wild))	// Compare degenerate bases
		{
			wild++;		// Move to next character
			inputseq++;
			//pos++;
		}
		else if ((*wild > '0') && (*wild <= '9'))	// This bit allows for gaps in the motif.
		{
			max_gap_size = (*wild - '0') + max_length;
			max_gap_size_plus_one = max_gap_size+1;
			while ((*wild) && (*wild > '0') && (*wild <= '9'))		// Remove all prefixed numbers
				wild++;
			buffer = new char[max_gap_size_plus_one];
			strncpy(buffer, inputseq, sizeof_char * max_gap_size_plus_one);
			buffer[max_gap_size] = null_char;
			match = IsMotifInSequence(buffer, wild, allowed_mutations - mutations);
			delete [] buffer;
			if (match)
				positions.push_back(pos);	// The call to the recursion returned true, so we should save the position.
			else
			{
				wild = mp;	// We do not have a match,
				cp++;		// so reset the wild string, and goto the next char on the inputsequence.
				inputseq = cp;
				pos++;
				mutations = 0;
			}
		}
		else if (mutations < allowed_mutations)
		{
			wild++;		// Move to next character
			inputseq++;
			mutations++;
		}
		else 
		{
			wild = mp;	// We do not have a match,
			cp++;		// so reset the wild string, and goto the next char on the inputsequence.
			inputseq = cp;
			pos++;
			mutations = 0;
		}
		if ((bool)!*wild)
		{
			positions.push_back(pos);
			wild = mp;	
			cp++;		
			inputseq = cp;
			pos++;
			mutations = 0;
		}
	}

	return positions;
}

void motif::get_max_length()
{
	unsigned int out = 0;
	const char* wild = sequence.c_str();
	while (*wild)
	{
		if ((*wild > '1') && (*wild <= '9'))
		{
			out += (*wild - '0');
		}
		else
			out++;
		wild++;
	}
	//return out;
	max_length = out;
}

void motif::CalcDegeneracy()
{
	unsigned int i;
	char base;
	degeneracy= 0.0;
	for (i=0; i<sequence.length(); i++)
	{
		base = sequence[i];
		switch (base)
		{
			case 'a': case 't': case 'g': case 'c':
				degeneracy++;
				break;
			case 'r': case 'y': case 'k': case 'm': case 's': case 'w':
				degeneracy += 2;
				break;
			case 'b': case 'd': case 'h': case 'v':
				degeneracy += 3;
				break;
			case 'n':
				degeneracy += 4;
				break;
			case '1': case '2': case '3': case '4': case '5': case '6': case '7':case '8': case '9':
				degeneracy += (4 * (base - '0' + 1));
				break;
			default:
				printf("Non-standard base in motif\n");
				cout << base << endl;
				cout << sequence << endl;
				abort();
		}
	}
	degeneracy /= (double)sequence.length();
}


char motif::GetRandomBase(bool is_nondegenerate)
{
	int letter;
	if (is_nondegenerate == true)
	{
		letter = ran2(0,6);
		if (letter == 4) letter = 14;
		if (letter == 5) letter = 15;
	}
	else
		letter = ran2(0, 16);
	switch(letter)
	{
	case 0:	return 'a'; break;
	case 1:	return 'c'; break;
	case 2:	return 'g'; break;
	case 3:	return 't'; break;
	case 4:	return 'r'; break;
	case 5:	return 'y'; break;
	case 6:	return 'k'; break;
	case 7:	return 'm'; break;
	case 8:	return 's'; break;
	case 9:	return 'w'; break;
	case 10: return 'b'; break;
	case 11: return 'd'; break;
	case 12: return 'h'; break;
	case 13: return 'v'; break;
	case 14: return 'n'; break;
	case 15: 
		{	// insert a gap of size 1 to 9.
			int i; char out;
			i=ran2(1, 10);
			out = '0' + i;
			return out;
			break;
		}

	default:
		printf("error - illegal nucleotide\n");
		abort();
		return '-';
	}
}


