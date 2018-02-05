
#include <string>

using std::string;

// Fasta_Read
////////////////////////////////////////////////////////

bool Fasta_Read (FILE * fp, string & s, string & hdr)
{
	int  ch;

	s . erase ();
	hdr . erase ();

	// skip till next '>' if necessary
	while  ((ch = fgetc (fp)) != EOF && ch != '>')
		;
	
	if  (ch == EOF)
		return  false;
	
	// put rest of line into  hdr
	while  ((ch = fgetc (fp)) != EOF && ch != '\n')
		hdr.push_back (char (ch));
	
	// put everything up till next '>' into  s
	while  ((ch = fgetc (fp)) != EOF && ch != '>')
	{
		if  (! isspace (ch))
			s.push_back (char (ch));
	}
	
	if  (ch == '>')
		ungetc (ch, fp);

	return  true;
}


// Fastq_Read
/////////////////////////////////////////////////////////////////////////


bool Fastq_Read (FILE * fp, string & s, string & hdr, string & qual)
{
        int  ch;

        s . erase ();
        hdr . erase ();
	qual.erase ();

	while  ((ch = fgetc (fp)) != EOF && ch != '@')
                ;

        if  (ch == EOF)
                return  false;

	while  ((ch = fgetc (fp)) != EOF && ch != '\n')
                hdr.push_back (char (ch));

	while  ((ch = fgetc (fp)) != EOF && ch != '+')
        {
                if  (! isspace (ch))
                	s.push_back (char (ch));
        }
	while  ((ch = fgetc (fp)) != EOF && ch == '\n')
		;

	ungetc (ch, fp);

	while  ((ch = fgetc (fp)) != EOF && ch != '\n')
        {
		 if  (! isspace (ch))
		        qual.push_back (char (ch));
	}

	//while  ((ch = fgetc (fp)) != EOF && ch != '@')
        //        ;

        if  (ch == '@')
                ungetc (ch, fp);

        return  true;
}



