// DNA.cpp: implementation of the DNA class.
//
//////////////////////////////////////////////////////////////////////

#include "DNA.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DNA::DNA()
{
}

DNA::DNA(char *seq)
{
	string seqstring = seq;
	std::transform (seqstring.begin(), seqstring.end(), seqstring.begin(), ToLower());
	Trim(seqstring);
	sequence.resize(1);
	sequence[0].append(seqstring);
}

DNA::~DNA()
{

}

// Read in only a specific sequence.
bool DNA::readFASTA(const char* filename, const string &header_name, bool uppercase_only, bool lowercase_only, bool KG_ref=false)
{
	bool output = true;
	char buffer[5120];
	int bufsize = 5120;
	string junk, junk2;
	junk.reserve(bufsize);
	junk2.reserve(bufsize);
	igzstream in;

	NSeqs = 0;
	in.open(filename, ios::in);
	if (in.good())
	{
		while(!in.eof())
		{
			in.getline(buffer, bufsize);
			if (buffer[0] == '>')
			{
				NSeqs++;
				string s = buffer;
				s = s.substr(1);
				if (KG_ref == true)
				{
					string::size_type pos=s.find(' ',0);
					s = s.substr(0,pos);
				}
				cout << s << endl;
				if (s == header_name)
					break;
			}
		}
		header.resize(1);
		sequence.resize(1);

		junk = buffer;
		if (KG_ref == true)
		{
			string::size_type pos=junk.find(' ',0);
			junk = junk.substr(0,pos);
		}
		header[0].append(junk.substr(1));

		in.getline(buffer, bufsize);
		junk = buffer;
		remove_upper_or_lowercase(junk, uppercase_only, lowercase_only);
		std::transform (junk.begin(), junk.end(), junk.begin(), ToLower());
		Trim(junk);
		while(!in.eof())
		{
			if (junk[0] == '>')
			{
				break;
			}
			else
			{
				junk2 = "";
				while ((junk[0] != '>') && (!in.eof()))
				{
					junk2.append(junk);
					in.getline(buffer, bufsize);
					junk = buffer;
					remove_upper_or_lowercase(junk, uppercase_only, lowercase_only);
					std::transform (junk.begin(), junk.end(), junk.begin(), ToLower());
					Trim(junk);
				}
				sequence[0].append(junk2);
			}
		}
		in.close();

		SeqLength.resize(1);
		SeqLength[0] = (unsigned int)(sequence[0].length());
	}
	else
	{
		cout << "Unable to read file: " << filename << endl;
		abort();
	}

	return output;
}


// Read all sequences
bool DNA::readFASTA(const char* filename, bool uppercase_only, bool lowercase_only, bool KG_ref=false)
{
	bool output = true;
	int count;
	int i;
	char buffer[5120];
	int bufsize = 5120;
	int nonstandard;
	string junk, junk2;
	junk.reserve(bufsize);
	junk2.reserve(bufsize);
	igzstream in, in2;

	NSeqs = 0;
	in2.open(filename, ios::in);
	if (in2.good())
	{
		while(!in2.eof())
		{
			in2.getline(buffer, bufsize);
			if (buffer[0] == '>')
			{
				string s = buffer;
				s = s.substr(1);
				if (KG_ref == true)
				{
					string::size_type pos=s.find(' ',0);
					s = s.substr(0,pos);
				}
				cout << s << endl;
				NSeqs++;
			}
		}
		header.resize(NSeqs);
		sequence.resize(NSeqs);
		cout << "\t" << NSeqs << " sequences" << endl;

		in2.close();

		in.open(filename, ios::in);
		in.getline(buffer, bufsize);
		junk = buffer;
		if (junk[0] != '>')
		{
			cout << "Error first line is not a headerline" << endl;
			exit(0);
		}
		Trim(junk);
		i=0;
		while(!in.eof())
		{
			if (junk[0] == '>')
			{
				if (KG_ref == false)
				{
					int pos = junk.find(" ");
					while (pos != -1)
					{
						junk.replace(pos,1,"_");
						pos = junk.find(" ");
					}
				}
				else
				{
					string::size_type pos=junk.find(' ',0);
					junk = junk.substr(0,pos);
				}
				header[i].append(junk.substr(1));
				i++;
				count = 0;
				in.getline(buffer, bufsize);
				junk = buffer;
				remove_upper_or_lowercase(junk, uppercase_only, lowercase_only);
				std::transform (junk.begin(), junk.end(), junk.begin(), ToLower());
				Trim(junk);
			}
			else
			{
				junk2 = "";
				while ((junk[0] != '>') && (!in.eof()))
				{
					junk2.append(junk);
					in.getline(buffer, bufsize);
					junk = buffer;
					Trim(junk);
					if (junk[0] != '>')
					{
						remove_upper_or_lowercase(junk, uppercase_only, lowercase_only);
						std::transform (junk.begin(), junk.end(), junk.begin(), ToLower());
						nonstandard = 1;
						while (nonstandard != -1)
						{
							nonstandard = (int)junk.find_first_not_of(" /t/nacgtn/0");
							if (nonstandard != -1)
							{
								cout << "Non-standard character found on line " << i << endl;
								junk[nonstandard] = 'n';
							}
						}
					}
				}
				sequence[i-1].append(junk2);
			}
		}
		in.close();

		SeqLength.resize(NSeqs);
		for (i=0; i < NSeqs; i++)
		{
			SeqLength[i] = (unsigned int)(sequence[i].length());
		}
	}
	else
	{
		cout << "Unable to open file: " << filename << endl;
		abort();
	}

	return output;
}

void DNA::Trim(string &str)
{
		// trim leading whitespace
		string::size_type  notwhite = str.find_first_not_of(" \t\n\v\b\r\f\a\\\?\'\"\0");
		str.erase(0,notwhite);

		// trim trailing whitespace
		notwhite = str.find_last_not_of(" \t\n\v\b\r\f\a\\\?\'\"\0"); 
		str.erase(notwhite+1); 
}

int DNA::CountMotifOccurances(motif &Motif)
{
	int i;
	int count = 0;
	bool match;

	for (i=NSeqs; i--;)
	{
		match = Motif.IsMotifInSequence(sequence[i].c_str(), Motif.sequence.c_str(), 0);
		if (match)
		{
			// have found motif so we're done
			count++;
		}
		else
		{
			// we didn't find the motif, so we need to test the complement
			match = Motif.IsMotifInSequence(sequence[i].c_str(), Motif.complement.c_str() ,0 );
			if (match)
			{
				count++;
			}
		}
	}
	return count;
}

void DNA::remove_upper_or_lowercase(string &string1, bool uppercase_only, bool lowercase_only)
{
	if (uppercase_only == true)
	{
		size_t position = string1.find_first_of("acgt");
		while ( position != string::npos )
		{
			string1.replace( position, 1, "N" );
			position = string1.find_first_of("acgt", position + 1 );
		}
	}
	else if (lowercase_only == true)
	{
		size_t position = string1.find_first_of("ACGT");
		while ( position != string::npos )
		{
			string1.replace( position, 1, "N" );
			position = string1.find_first_of("ACGT", position + 1 );
		}
	}
}
