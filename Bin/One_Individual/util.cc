
#include <string>

using std::string;

/// Split
//////////////////////////////////////////////////////////

void split(string &str, char delim, std::vector<string> &splitstr)
{
    string buf = "";
    int i = 0;
    while (i < str.length()) {
        if (str[i] != delim)
            buf += str[i];
        else if (buf.length() > 0) {
            splitstr.push_back(buf);
            buf = "";
        }
        i++;
    }
    if (!buf.empty())
        splitstr.push_back(buf);
}


// To_UPPER_case
// /////////////////////////////////////////////////////

string to_upper_seq(string &seq)
{
	string up;
        for ( string::iterator it=seq.begin(); it != seq.end(); ++it)
        {
        	char letter = *it;
        	switch (letter)
        	{
            	case 'a':
	        	up.append("A"); break;
               	case 't':
                       	up.append("T"); break;
                case 'g':
                       	up.append("G"); break;
                case 'c':
                       	up.append("C"); break;
	        case 'n':
                       	up.append("N"); break;
              	case 'A':
	        	up.append("A"); break;
               	case 'T':
                       	up.append("T"); break;
                case 'G':
                       	up.append("G"); break;
                case 'C':
                       	up.append("C"); break;
                default:
		       	up.append("N"); break;
        	}
        }
	return(up);
}



// Reverse_Complement
// /////////////////////////////////////////////////////

string reverse_complement(string &seq)
{
        string rc;
        for ( string::iterator it=seq.end()-1; it != seq.begin()-1; --it)
        {
        	char letter = *it;
        	switch (letter)
        	{
        // try to find the complement of every char of the source string
                	case 'A':
	                	rc.append("T"); break;
               	 	case 'T':
                        	rc.append("A"); break;
                	case 'G':
                        	rc.append("C"); break;
                	case 'C':
                        	rc.append("G"); break;
	        	case 'N':
                        	rc.append("N"); break;
              	  	default:
                        	rc.append("N"); break;
        	}
        }
	return(rc);
}


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


bool Fastq_Read (FILE * fp, string & s, string & hdr)
{
        int  ch;

        s . erase ();
        hdr . erase ();

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

	while  ((ch = fgetc (fp)) != EOF && ch != '@')
                ;

        if  (ch == '@')
                ungetc (ch, fp);

        return  true;
}



