/**
 * @file    ExtraTools.h
 * @author  Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                                  <afshinfard (at) gmail.com>
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




#endif // EXTRATOOLS_H
