/**
 * @file    ExtraTools.h
 * @author  Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                                  <a.afshinfard (at) gmail.com>
 *
 * @section LICENCE
 *
 * Copyright (C) 2017-2020
 *   Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                           <afshinfard (at) gmail.com>
 *	 Seyed Abolfazl Motahari <motahari (at) sharif.edu
 *
 **/


#ifndef EXTRATOOLS_H
#define EXTRATOOLS_H


#include <string>
using namespace std;

inline string revComplemACGT(string readSegment){
    // [Completed]
    for (char& c : readSegment) {
        switch(c) {
        case 'a':
            c = 't';
            break;
        case 'g':
            c = 'c';
            break;
        case 'c':
            c = 'g';
            break;
        case 't':
            c = 'a';
            break;
        case 'A':
            c = 'T';
            break;
        case 'G':
            c = 'C';
            break;
        case 'C':
            c = 'G';
            break;
        case 'T':
            c = 'A';
            break;
        default:
            ;
        }
    }
    reverse( readSegment.begin(), readSegment.end() );
    //readSegment = string ( readSegment.rbegin(), readSegment.rend() );
    return readSegment;
}


string convertNumToStr(long long number) {
    ostringstream ostr;
    ostr<<number;
    return ostr.str();
}


inline string getGenome(long long loci, long long length){
    // [not yet]
    if(loci + length > Ref.sz ){
        cerr<<"\n=======\n=======\n=======\n exceeds genome length!\n=======\n=======\n=======\n";
        return 0;
    }
    string genomeString;
    for(long long i = 0 ; i < length ; i++ ){
        genomeString += Ref.charAt(loci+i);
    }
    return genomeString;
}

int findCharacterPerRow(string genomeName) {
    int characterPerRow;
    ifstream ifstr(genomeName.c_str());
    string line;
    getline(ifstr, line);
    getline(ifstr, line);
    characterPerRow = (long)line.size();
    ifstr.close();
    return characterPerRow;
}
long long getGenomeLength(){

    ifstream file((MAINADDRESS+genomeName).c_str(), ifstream::in | ifstream::binary);
    if(!file.is_open())
        return -1;

    string line1;
    getline(file,line1);
    int line1Size = line1.size();

    file.seekg(0, ios::end);
    long long fileSize = file.tellg();
    file.close();

    fileSize = fileSize - line1Size - 1 - 1 ;/*the last line has no \n*/
    long remainedLineCount = (fileSize)/81 ;
    long long GenomeSize = fileSize - remainedLineCount;
    return GenomeSize;
}
string getGenomeHDD(long long first, int length) {
    // [Complete]
    struct stat sb;
    off_t len;
    char *p;
    int fd;
    string name = genomeName;// + ".fasta";

    getGenomeHDDMutex.lock();

    ifstream ifstr(name.c_str());
    string firstLine;
    getline(ifstr, firstLine);
    int firstLineCharCount = (int)firstLine.size();
    ifstr.close();
    first += firstLineCharCount+1 + first/characterPerRow;
    fd = open(name.c_str(), O_RDONLY);
    if (fd == -1)
        perror ("open");

    if (fstat (fd, &sb) == -1)
        perror ("fstat");

    p = (char *)mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

    if (close (fd) == -1)
        perror ("close");

    string temp = "";
    for (len = first; temp.size() < length; ){
        if (p[len] != '\n')
            temp += p[len];

        ++len;
    }
    if (munmap (p, sb.st_size) == -1)
        perror ("munmap");
    //close(fd);
    getGenomeHDDMutex.unlock();
    return temp;
}



#endif // EXTRATOOLS_H
