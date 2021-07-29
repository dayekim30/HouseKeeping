#include "Sequence.h"

string Sequence::forwardSeq()
{	
    string str(seq.length(), 0);
    for (int c = 0; c < str.length(); c++) {
        char base = seq[c];
        switch (base) {
        case 'A':
        case 'a':
            str[c] = 'A';
            break;
        case 'T':
        case 't':
            str[c] = 'T';
            break;
        case 'C':
        case 'c':
            str[c] = 'C';
            break;
        case 'G':
        case 'g':
            str[c] = 'G';
            break;
        default:
            str[c] = 'N';
        }
    }
    return str;
	
}

string Sequence::backwardSeq()
{
    string str(seq.length(), 0);
    for (int c = 0; c < seq.length(); c++) {
        char base = seq[c];
        switch (base) {
        case 'A':
        case 'a':
            str[seq.length() - c - 1] = 'T';
            break;
        case 'T':
        case 't':
            str[seq.length() - c - 1] = 'A';
            break;
        case 'C':
        case 'c':
            str[seq.length() - c - 1] = 'G';
            break;
        case 'G':
        case 'g':
            str[seq.length() - c - 1] = 'C';
            break;
        default:
            str[seq.length() - c - 1] = 'N';
        }
    }
    return str;
}
