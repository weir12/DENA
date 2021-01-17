
#include "motiflocation.h"

long *idum;

struct ToLower
{
	char operator() (char c) const  { return tolower(c); }
};

int main(int argc, char* argv[])
{
	int i;
	unsigned int ui;
	char filename1[255], outfilename[255];
	string startstr = "";
	string chrom_str = "";
	DNA SequenceSet1;
	bool f1 = false, f2 = false;
	int allowed_mutations = 0;
	bool chosen_chr = false;
	bool output_seqlen = false;
	bool uppercase_only = false, lowercase_only = false;
	bool KG_ref = false;

	for(i = 0; i < argc; i++)
	{
		if(*argv[i] == '-')
		{ 
			if((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-?") == 0) || (strcmp(argv[i], "-help") == 0)) 
			{
				print_help();
				exit (0);
			}
			if(strcmp(argv[i], "-file1") == 0) { strcpy(filename1 , argv[i+1]); f1 = true; }
			if(strcmp(argv[i], "-out") == 0) { strcpy(outfilename , argv[i+1]); f2 = true; }
			if(strcmp(argv[i], "-startstr") == 0) startstr = argv[i+1];
			if(strcmp(argv[i], "-chrom") == 0) {chrom_str = argv[i+1]; chosen_chr = true; }
			if(strcmp(argv[i], "-mutations") == 0) allowed_mutations = atoi(argv[i+1]);
			if(strcmp(argv[i], "-outputseqlen") == 0) output_seqlen = true;
			if(strcmp(argv[i], "-uppercase_only") == 0) uppercase_only = true;
			if(strcmp(argv[i], "-lowercase_only") == 0) lowercase_only = true;
			if(strcmp(argv[i], "-1000G_ref") == 0) KG_ref = true;
		}
	}

	if (chosen_chr == false)
	{
		chrom_str = "n.a";
	}

	if (!((f1 == true) && (f2 == true)))
	{
		cout << "Missing input parameter" << endl;
		print_help();
		exit(0);
	}

	
	std::transform(startstr.begin(), startstr.end(), startstr.begin(), ToLower());

	motif my_motif = motif(startstr);

	cout << "\nReading in Files" << endl;
	if ((strstr( filename1, ".fa" ) != NULL) || (strstr( filename1, ".fa.gz" ) != NULL))
	{
		cout << "File 1 Fasta file" << endl;
		if (chosen_chr == false)
			SequenceSet1.readFASTA(filename1, uppercase_only, lowercase_only, KG_ref);
		else
			SequenceSet1.readFASTA(filename1, chrom_str, uppercase_only, lowercase_only, KG_ref);
	}
	else
	{
		cout << "Appears to not be a fasta file? Should end in .fa" << endl;
		exit(0);
	}

	cout << "Finished reading files\n" << endl;
	for (unsigned int ui=0; ui<SequenceSet1.sequence.size(); ui++)
		cout << "Read " << SequenceSet1.sequence[ui].size() << " bases in " << SequenceSet1.header[ui] << endl;

	cout << "Searching for motif" << endl;
	vector<string> chroms;
	vector<int> positions;
	vector<int> seqnum;
	vector<string> revchroms;
	vector<int> reverse_positions;
	vector<int> revseqnum;

	for (unsigned int ui=0; ui<SequenceSet1.sequence.size(); ui++)
	{
		vector<int> tmp_positions;
		tmp_positions = my_motif.returnMotifLocations(SequenceSet1.sequence[ui].c_str(), my_motif.sequence.c_str(), allowed_mutations);
		vector<string> tmp_chroms(tmp_positions.size(), SequenceSet1.header[ui]);
		positions.insert( positions.end(), tmp_positions.begin(), tmp_positions.end() );
		chroms.insert(chroms.end(), tmp_chroms.begin(), tmp_chroms.end() );
		vector<int> tmp_seqnum(tmp_positions.size(), ui);
		seqnum.insert( seqnum.end(), tmp_seqnum.begin(), tmp_seqnum.end());
	}

	cout << "Searching for reverse complement" << endl;
	for (unsigned int ui=0; ui<SequenceSet1.sequence.size(); ui++)
	{
		vector<int> tmp_positions;
		tmp_positions = my_motif.returnMotifLocations(SequenceSet1.sequence[ui].c_str(), my_motif.complement.c_str(), allowed_mutations);
		vector<string> tmp_chroms(tmp_positions.size(), SequenceSet1.header[ui]);
		reverse_positions.insert( reverse_positions.end(), tmp_positions.begin(), tmp_positions.end() );
		revchroms.insert(revchroms.end(), tmp_chroms.begin(), tmp_chroms.end() );
		vector<int> tmp_seqnum(tmp_positions.size(), ui);
		revseqnum.insert( revseqnum.end(), tmp_seqnum.begin(), tmp_seqnum.end());
	}


	ofstream out_file;
	out_file.open(outfilename);
	//out_file << "chrom\tchromStart\tchromEnd\tstrand\tsequence" << endl;
	int plus = startstr.length();
	for (ui=0; ui < positions.size(); ui++)
	{
		out_file << chroms[ui] << "\t" << positions[ui] << "\t" << positions[ui] + plus << "\t+\t" << SequenceSet1.sequence[seqnum[ui]].substr(positions[ui], plus);
		if (output_seqlen == true)
			out_file << "\t" << SequenceSet1.sequence[seqnum[ui]].size();
		out_file << endl;
	}

	for (ui=0; ui < reverse_positions.size(); ui++)
	{
		out_file << revchroms[ui] << "\t" << reverse_positions[ui] << "\t" << reverse_positions[ui] + plus << "\t-\t" << SequenceSet1.sequence[revseqnum[ui]].substr(reverse_positions[ui], plus);
		if (output_seqlen == true)
			out_file << "\t" << SequenceSet1.sequence[revseqnum[ui]].size();
		out_file << endl;
	}
	out_file.close();
	
	return 1;
}


void print_help() {

	printf("\nmotifhunter %0.1f\n", version);
	printf("Returns every location of a motif in a given file.\n\n");
	printf("*** Positions are zero based! ***.\n\n");
	
	printf("Options :\n");
	printf("-file1 <file>                Sequence data file 1 (required)\n");
	printf("-out <file>                  Output file (required)\n");
	printf("-startstr <string>           Starting Motif\n");
	printf("-chrom <string>              Output string used for chromosome\n");
	printf("-mutations <int>             Number of allowed mutations\n");
	printf("-outputseqlen                Output the chromosome length as well\n");
	printf("-uppercase_only              Only output results in uppercase sequence\n");
	printf("-lowercase_only              Only output results in lowercase sequence\n");
	printf("-1000G_ref                   FASTA headers are in 1000G format\n");
	printf("\n");
}


