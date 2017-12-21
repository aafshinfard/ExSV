#ifndef WRITERESULTS_H
#define WRITERESULTS_H

#include <unistd.h>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <fstream>


void write2Informatives( vector<informativeReadsData> *informativeReads/*, ofstream *ofstr1*/ ){
    // [Complete]

    //cerr<<"\n\n 2inf++ "<<informativeReads->size()<<endl;
    cnt_2informatives += informativeReads->size();
    for ( vector<informativeReadsData>::iterator it = informativeReads->begin() ; it < informativeReads->end(); it++){

        for ( vector<solvedFragInfo>::iterator it2 = it->clusters.begin() ; it2 < it->clusters.end(); it2++){
            (ofstr_2Informatives)<< (it->readName) <<"\t|\t"
                     << (it2->startChunk)            <<  "\t"  <<  (it2->endChunk) <<  "\t"
                     << (it2->isReverse ? "-1":"+1") <<  "\t"  <<  (it2->shift)    <<  "\t"
                     << (it2->startPos)              <<  "\t"  <<  (it2->endPos)   <<  endl;
        }
    }

}

void writeInformatives( vector<informativeReadsData> *informativeReads/*, ofstream *ofstr1*/ ){
    // [Complete]

    //cerr<<"\n\n inf++ "<<informativeReads->size()<<endl;
    cnt_informatives += informativeReads->size();
    for ( vector<informativeReadsData>::iterator it = informativeReads->begin() ; it < informativeReads->end(); it++){

        //ofstr1<< (it->readName) << endl;
        //ofstr1<< (it->readSeq) << endl;

        //style 1:
        //        ofstr1<< (it->index) <<"\t|\t";
        //        for ( vector<solvedFragInfo>::iterator it2 = it->clusters.begin() ; it2 != it->clusters.end(); ++it2)
        //            ofstr1<< (it2->startChunk)            <<  "\t"  <<  (it2->endChunk) <<  "\t"
        //                  << (it2->isReverse ? "-1","+1") <<  "\t"  <<  (it2->shift)    <<  "\t"
        //                  << (it2->startPos)              <<  "\t"  <<  (it2->endPos)   <<  "\t|\t";
        //        ofstr1<<endl;

        //style 2:
        for ( vector<solvedFragInfo>::iterator it2 = it->clusters.begin() ; it2 < it->clusters.end(); it2++){
            (ofstr_Informatives)<< (it->readName) <<"\t|\t"
                     << (it2->startChunk)            <<  "\t"  <<  (it2->endChunk) <<  "\t"
                     << (it2->isReverse ? "-1":"+1") <<  "\t"  <<  (it2->shift)    <<  "\t"
                     << (it2->startPos)              <<  "\t"  <<  (it2->endPos)   <<  endl;
        }
    }
}

void writeUnMapped( vector<unMappedReads>* unMappeds/*, ofstream *ofstr1*/ ){
    // [Complete]
    writeUnMappedMutex.lock();
    cnt_unmapped += unMappeds->size();
    for ( vector<unMappedReads>::iterator it = unMappeds->begin() ; it != unMappeds->end(); ++it)
        (ofstr_unMapped)<< (it->name) <<"\n"<< (it->read) <<"\n"<< (it->quality) <<endl;
    writeUnMappedMutex.unlock();
}


#endif // WRITERESULTS_H
