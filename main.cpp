/**
 * @file    main.cpp
 * @author  Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                                  <afshinfard (at) gmail.com>
 *          Bioinformatics Research Lab - Sharif Uiversity of Technology
 *
 * @cite
 *
 * @copyright (c) 2017
 *
 *      Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                              <afshinfard (at) gmail.com>
 *      Damoon Nashta-ali       <damoun_dna (at) yahoo.com>
 *      Seyed Abolfazl Motahari <motahari (at) sharif.edu
 *
 **/

#include "Header.h"
#include "ExtraTools.h"
#include "DataStructures.h"
#include "SureMap.h"
#include "LocalAligner.h"

// ================================|
//> Global Vars:                   |
// ================================|
//                                 |
int characterPerRow;
long long rdcntr = 0;
long long kcntr = 0;
int fileCounter = 1; // to remember current active file
long long cnt_2informatives = 0;
long long cnt_informatives = 0;
long long cnt_unmapped = 0;
std::mutex rdcntrMutex;
std::mutex readFileMutex;
std::mutex getGenomeHDDMutex;
std::mutex writeInformativesMutex;
std::mutex write2InformativesMutex;
std::mutex writeUnMappedMutex;
long long totalReads = 0;
long long readCounter = 0;
long long genomeLength = 0;
readMappingInfo globalRMI[MAXTHREADS];

ofstream ofstr_Informatives(informativeFileName.c_str());
ofstream ofstr_2Informatives(doubInformativeFileName.c_str());
ofstream ofstr_unMapped(unMappedFileName.c_str());

int shiftSteps = chunkSize / shiftIterations ;
//                                 |
// ================================'

// ================================|
//> Functions:                     |
// ================================|
//                                 |

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
mapResult uniqueMap( string readChunk , string quality, int idx ){
    // [Completed]
    // Func aim:
    // given a read chunk, it's quality an thread index
    // return the unique mapping result for this chunk

    globalRMI[idx].read = readChunk;
    globalRMI[idx].qual = quality;
    // globalRMI[idx].readName = "?";

    vector< mappingInfo > mps = getAllMapping( idx, globalRMI[idx] );

    string mos = readChunk;
    string rev = reverseComplement( readChunk );
    if( mps.size() > 1 ){
        cerr<<"not a unique search";
        return mapResult();
    }
    mappingInfo ix = mps[0];
    if( ix.flag != 4 ){
        ix.cigar = getCigar( ix.refPos, ( ix.flag == 0 ) ? mos : rev, ix.acceptedValue, idx  );
        ix.refPos = ix.refPos - ix.acceptedValue + ix.cigar.second;
    }else
        return mapResult( -1 , false  );
    return mapResult(ix.refPos,( ix.flag == 0 ? false : true ));

}

string getGenome(){
    // [not yet]
}

int findCharacterPerRow(string genomeName) {
    int characterPerRow;
    ifstream ifstr(genomeName.c_str());
    string line;
    getline(ifstr, line);
    getline(ifstr, line);
    characterPerRow = (int)line.size();
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


int ifAlignable( long long startLoci, string readSegment, bool isReverse, int *alignShiftPos
                  /*,bool isToRight, jumpsCount, (counter==0 && onlyBG ? true:false)*/ ){
    if( isReverse )
        readSegment = revComplemACGT( readSegment );
    double indelShiftLocal = indelShift /*+ (numGap*jumpsCount)/chunkSize*/;

//    bool binGenome;
//    if( useBinGenome ){
//        binGenome = checkBinGenome(startLoci, startLoci + this->chunkSize, false, (BGLR ? isToRight : true));
//        if(binGenome){
//            cntBG++;
//            return 1;
//        }
//        if( onlyBGlocal )
//            return 0;
//    }
    cntCOM++;

    if(startLoci - indelShiftLocal*readSegment.size() < 0 || startLoci + indelShiftLocal*readSegment.size()+chunkSize  > genomeLength-1 )
        return 0;

    string genomeString = getGenomeHDD(startLoci - indelShiftLocal*readSegment.size(), chunkSize + 2*readSegment.size()*indelShiftLocal );
    int readlen = readSegment.size();
    string seqConstantB = "";

    for(int i=0;i < (int)(indelShiftLocal*readlen);i++)
        seqConstantB += "A";
    readSegment += seqConstantB;

    LocalAligner *local;
    local = new LocalAligner(readSegment, genomeString, gapPen, misPen, matchPen, (int)(indelShiftLocal*readlen), (int)(editDistance*readlen) );
    local->process();
    local->backtrack();
    //cout//<<local->realCigar << endl;
    //cout<<local->cigar << endl;
    int score = local->mScore;
    char isAlignable = (score > localAlThreshold ? 2 : 0 );
    if( isAlignable ){
        *alignShiftPos = local->totalGapr - local->totalGapfa;
        //previous shift if needed.
    }
    delete local;
    return isAlignable;
}

long extendManager(long from, long until, string* read, solvedFragInfo* anchorExtendResult){
    // [Completed] FUNC aim:
    // to extend from from's neighbour (not from from itself) to until
    // from < until not guaranteed! based n diretion it's feasible and return number of last extended chunk
    if(from == until)
        return from;
    long long startLoci = (from > until ? anchorExtendResult->startPos :anchorExtendResult->endPos);
    int direction = (from > until ? -1 : +1);    
    bool isReverse = anchorExtendResult->isReverse;

    // ====================================
    // frags of read from from's neighour until until:
    vector<string> fragsOfRead;

    if(direction == 1){
        for( int i = from; i < (until-1) ; i++ )
            fragsOfRead.push_back( read->substr(  (i*chunkSize)  , chunkSize ));

        if( until == (int) ((read->size() + chunkSize - 20) / chunkSize) )
            fragsOfRead.push_back( read->substr( (until-1)*chunkSize /*, readLen-((readLen-1)/chunkSize)*chunkSize)*/ ));
        else
            fragsOfRead.push_back( read->substr( (until-1)*chunkSize , chunkSize ));
    }else
        for( int i = from-1; i >= until ; i-- )
            fragsOfRead.push_back( read->substr( ((i-1)*chunkSize)  , chunkSize ));

    // ===============================================
    // let's extend from from until until if possible:
    int counter = 0;
    int alignShiftPos = 0;
    int successFlag = 0;
    do{
        startLoci += (direction*(chunkSize + alignShiftPos))*( isReverse ? -1 : 1 );
        successFlag =  ifAlignable( startLoci, fragsOfRead.at( counter ), isReverse, &alignShiftPos
                                    /*,isToRight, jumpsCount, (counter==0 && onlyBG ? true:false)*/ );
        counter += ( successFlag ? 1:0 );
        //jumpsCount += ( successFlag == 1 ? 1:0);
    }while( successFlag && counter < abs(from - until) );

    // ========================
    // data changes and return:
    if(from > until){
        anchorExtendResult->startChunk = from - counter;
        anchorExtendResult->startPos -= (isReverse?-1:+1)*(counter*chunkSize+alignShiftPos);
        return anchorExtendResult->startChunk;
    }else{
        anchorExtendResult->endChunk = from + counter;
        anchorExtendResult->endPos += (isReverse?-1:+1)*(counter*chunkSize+alignShiftPos);
        return anchorExtendResult->endChunk;
    }
}

void fullExtension( string* read, vector<solvedFragInfo>* clusters, int idx ){
    int numFrags =  ( read->size() + chunkSize - 20 ) / chunkSize;
    sort( clusters->begin(), clusters->end() );

    if( clusters->front().startChunk > 1 ){
        extendManager( clusters->front().startChunk ,1 ,read ,&(clusters->front()) );
    }
    for( int i = 1 ; i < clusters->size() ; i++ ){
        if(clusters->at(i-1).endChunk > clusters->at(i).startChunk){
            cerr<<"\n baha'e"<<endl;
            continue;
        }
        extendManager( clusters->at(i-1).endChunk ,clusters->at(i).startChunk-1 ,read ,clusters->data()+i-1 );
        extendManager( clusters->at(i).startChunk ,clusters->at(i-1).endChunk+1 ,read ,&(clusters->at(i)) );
    }
    if( clusters->back().endChunk < numFrags ){
        extendManager( clusters->back().endChunk ,numFrags ,read ,&(clusters->back()) );
    }

}

interval anchorAndExtend(string* read ,string* qc ,vector<interval>* checkStack ,int it
                         ,array<vector<mapResult>, MAX_SHIFT_ITER>* mapResults
                         ,vector<solvedFragInfo>* clusters
                         ,int idx ){
    // Remained: inExtend
    // FUNC aim:
    // Given [a single read], its [quality], its [incomplete mapping records] for different shifting
    // and [asked Intervals] try to anchor and exend exactly for one cluster and push new intervals to checkStack


    solvedFragInfo anchorExtendResult;
    interval toCheck = checkStack->back();
    checkStack->pop_back();

    bool isFinalChunk = ( ((int)(read->length() + chunkSize - 20) / chunkSize) == toCheck.end ? true : false );

    // =====================
    // Check First and Last:

    if( mapResults->at(it)[toCheck.start-1].position == 0 )
        mapResults->at(it)[toCheck.start-1] = uniqueMap( read->substr((toCheck.start-1)*chunkSize, chunkSize)
                                                       ,qc->substr((toCheck.start-1)*chunkSize, chunkSize)
                                                       ,idx );
    if( mapResults->at(it)[toCheck.end-1].position == 0 )
        mapResults->at(it)[toCheck.end-1] = uniqueMap( ( isFinalChunk ? read->substr((toCheck.end-1)*chunkSize)
                                                                  : read->substr((toCheck.end-1)*chunkSize, chunkSize))
                                                     ,( isFinalChunk ? qc->substr((toCheck.end-1)*chunkSize)
                                                                     : qc->substr((toCheck.end-1)*chunkSize, chunkSize))
                                                     ,idx );
    if( mapResults->at(it)[toCheck.start-1].isReverse == mapResults->at(it)[toCheck.end-1].isReverse &&
            abs(mapResults->at(it)[toCheck.start-1].position - mapResults->at(it)[toCheck.end-1].position
            - pow(-1,(mapResults->at(it)[toCheck.start-1].isReverse == true ? 1:0 ))*(toCheck.start-toCheck.end)*chunkSize )
            <
            gapRate*abs(toCheck.start-toCheck.end) ){
            // nothing to push to stack, this interval has been completely solved.
            anchorExtendResult = solvedFragInfo(toCheck.start,mapResults->at(it)[toCheck.start-1].position,
                                                   toCheck.end, mapResults->at(it)[toCheck.end-1].position,
                                                   mapResults->at(it)[toCheck.start-1].isReverse, it * shiftSteps );
            clusters->push_back(anchorExtendResult);
            return toCheck;
    }

    // =====================
    // Check others:

    interval anchoredInterval;
    //interval extendedInterval;
    bool isAnchored = false;
    int half = (int) ceil((double) (toCheck.end-toCheck.start)/2);
//    if ((toCheck.end-toCheck.start) %2 != 0) {
//        half--;
//    }

    for (int i = 1; i < half; i++) {
        // one from left with all rights
        if( mapResults->at(it)[toCheck.start+i-1].position == 0 )
            mapResults->at(it)[toCheck.start+i-1] = uniqueMap( read->substr((toCheck.start+i-1)*chunkSize, chunkSize)
                                                             ,qc->substr((toCheck.start+i-1)*chunkSize, chunkSize)
                                                             ,idx );
        if( mapResults->at(it)[toCheck.start+i-1].position > 0 ){
            for(int j = 0 ; j < i ; j++){
                if( mapResults->at(it)[toCheck.start+i-1].isReverse == mapResults->at(it)[toCheck.end-j-1].isReverse &&
                        abs(mapResults->at(it)[toCheck.start+i-1].position - mapResults->at(it)[toCheck.end-j-1].position
                            - pow(-1,(mapResults->at(it)[toCheck.start+i-1].isReverse==true?1:0))*(toCheck.start+i-(toCheck.end-j))*chunkSize )
                        <
                        gapRate*abs(toCheck.start+i-toCheck.end-j) ){
                    anchoredInterval = interval( toCheck.start+i, toCheck.end-j);
                    isAnchored = true;
                    anchorExtendResult = solvedFragInfo(toCheck.start+i,mapResults->at(it)[toCheck.start+i-1].position,
                                                           toCheck.end-j, mapResults->at(it)[toCheck.end-j-1].position,
                                                           mapResults->at(it)[toCheck.start+i-1].isReverse, it * shiftSteps );
                    // add to readAnchorInfo DS
                    break;
                }
            }
        }
        if(isAnchored)
            break;
        // one from right with all lefts
        if( mapResults->at(it)[toCheck.end-i-1].position == 0 )
            mapResults->at(it)[toCheck.end-i-1] = uniqueMap( read->substr((toCheck.end-i-1)*chunkSize, chunkSize)
                                                           ,qc->substr((toCheck.end-i-1)*chunkSize, chunkSize)
                                                           ,idx );
        if( mapResults->at(it)[toCheck.end-i-1].position > 0 ){
            for(int j = 0 ; j <= i ; j++){
                if( mapResults->at(it)[toCheck.start+j-1].isReverse == mapResults->at(it)[toCheck.end-i-1].isReverse &&
                        abs(mapResults->at(it)[toCheck.start+j-1].position - mapResults->at(it)[toCheck.end-i-1].position
                            - pow(-1,(mapResults->at(it)[toCheck.start+j-1].isReverse==true?1:0))*(toCheck.start+j-(toCheck.end-i))*chunkSize )
                        <
                        gapRate*abs(toCheck.start+j-toCheck.end-i) ){
                    anchoredInterval = interval( toCheck.start+j, toCheck.end-i);
                    isAnchored = true;
                    anchorExtendResult = solvedFragInfo(toCheck.start+j,mapResults->at(it)[toCheck.start+j-1].position,
                                                           toCheck.end-i, mapResults->at(it)[toCheck.end-i-1].position,
                                                           mapResults->at(it)[toCheck.start+j-1].isReverse, it * shiftSteps );
                    break;
                }
            }
        }
        if(isAnchored)
            break;
    }
    // =======================
    // Check Middle if needed:
    if( (toCheck.end-toCheck.start)%2 == 0 && !isAnchored ){
        // check Middle with all
        long middle = toCheck.start+half-1;
        if( mapResults->at(it)[middle].position == 0 )
            mapResults->at(it)[middle] = uniqueMap( read->substr((toCheck.start+half-1)*chunkSize, chunkSize)
                                                    ,qc->substr((toCheck.start+half-1)*chunkSize, chunkSize)
                                                    ,idx );
        if( mapResults->at(it)[middle].position > 0 ){
            for (int j = 0; j < half; j++) {
                // with left
                if( mapResults->at(it)[toCheck.start+j-1].isReverse == mapResults->at(it)[middle].isReverse &&
                        abs(mapResults->at(it)[toCheck.start+j-1].position - mapResults->at(it)[middle].position
                            - pow(-1,(mapResults->at(it)[toCheck.start+j-1].isReverse==true?1:0))*(toCheck.start+j-middle)*chunkSize )
                        <
                        gapRate*abs(j-middle) ){
                    anchoredInterval = interval( toCheck.start+j, middle);
                    isAnchored = true;
                    anchorExtendResult = solvedFragInfo(toCheck.start+j,mapResults->at(it)[toCheck.start+j-1].position,
                                                           middle, mapResults->at(it)[middle].position,
                                                           mapResults->at(it)[toCheck.start+j-1].isReverse, it * shiftSteps );
                    break;
                }
                //with right
                if( mapResults->at(it)[toCheck.end-j-1].isReverse == mapResults->at(it)[middle].isReverse &&
                        abs(mapResults->at(it)[toCheck.end-j-1].position - mapResults->at(it)[middle].position
                            - pow(-1,(mapResults->at(it)[toCheck.end-j-1].isReverse==true?1:0))*(toCheck.end-j-middle)*chunkSize )
                        <
                        gapRate*abs(j-middle) ){
                    anchoredInterval = interval( middle, toCheck.end-j);
                    isAnchored = true;
                    anchorExtendResult = solvedFragInfo(toCheck.end-j,mapResults->at(it)[toCheck.end-j-1].position,
                                                           middle, mapResults->at(it)[middle].position,
                                                           mapResults->at(it)[toCheck.end-j-1].isReverse, it * shiftSteps );
                    break;
                }

            }
        }
    }
    // ========================
    // Recursivley Check sides:
    if(!isAnchored){
        interval tempa = interval(toCheck.start,toCheck.start+half-(1) );
        if( abs(tempa.start - tempa.end) > 0 )
            checkStack->push_back( tempa  );
        interval tempb = interval(toCheck.end-half+(1) , toCheck.end );
        if( abs(tempb.start - tempb.end) > 0)
            checkStack->push_back( tempb );
        return interval(0,0);
    }
    // ONLY Extend-needed ANCHORED Intervals:

    // =====================
    // Extend if Possible: (only to check intervals)
    //1: inExtend
    // not yet

    //2.a: outExtend left
    //    if(anchoredInterval.start != toCheck.start ){
    //        extendedInterval.start = extendManager( anchoredInterval.start, toCheck.start, read, &anchorExtendResult);
    //        // push unsolved regions:
    //        if((extendedInterval.start - toCheck.start) > 1)
    //            checkStack->push_back( interval(toCheck.start, extendedInterval.start-1) );
    //    }else
    //        extendedInterval.start = toCheck.start;
    //    //2.b: outExtend right
    //    if(anchoredInterval.end != toCheck.end ){
    //        extendedInterval.end = extendManager( anchoredInterval.end, toCheck.end, read, &anchorExtendResult);
    //        // push unsolved regions:
    //        if((toCheck.end - extendedInterval.end) > 1)
    //            checkStack->push_back( interval(extendedInterval.end+1, toCheck.end) );
    //    }else
    //        extendedInterval.end = toCheck.end;

    // =====================
    // return solved interval
    clusters->push_back(anchorExtendResult);
    //return extendedInterval;
    return anchoredInterval;
}

void inline determineCheckStack( array<vector<bool>, MAX_SHIFT_ITER>* isSolvedChunk , int it,
                                 array<vector<mapResult>, MAX_SHIFT_ITER>* mapResults){
    int shiftLen = it*shiftSteps;
    // previous chunkings: 0 to it in isSolvedChunk
    // nothing yet
}

void processing( long long* readIndex ,string* readsName ,string* reads ,string *qc ,long long cntReads ,int idx ){
    // [Complete itself]

    //ofstream ofstr_Informatives(informativeFileName.c_str());
    //ofstream ofstr_2Informatives(doubInformativeFileName.c_str());
    //ofstream ofstr_unMapped(unMappedFileName.c_str());
    //vector<readInfo> informativeReads;
    string currentRead, shiftedRead, currentQC;
    int shiftLen = 0;
    vector<informativeReadsData> informativeReads;
    vector<informativeReadsData> doubInformatives;
    vector<unMappedReads> unMappeds;
    // each read data strucures
    vector<solvedFragInfo> clusters;
    vector<interval> checkStack;
    array<vector<bool>, MAX_SHIFT_ITER> isSolvedChunk; // ALTERNATIVECODING: static data structure:
    array<vector<mapResult>, MAX_SHIFT_ITER> mapResults;// 0:not searched yet - -1: can't find uniqely - -2: can't fnd with least size
    long cntChunks = 0;

    for(long i = 0 ; i < cntReads ; i++){ // loop on reads of current thread
                                         // =============================== {
        // fullyAnchored = false;
        if( !clusters.empty() )
            clusters.clear();
        for(int j = 0 ; j < shiftIterations ; j++ ){ // loop on shift anchoring
                                                   // =============================== {
            // [[[Initialization: ]]]
            shiftLen = j*shiftSteps;
            currentRead = reads[i].substr( shiftLen );
            currentQC = qc[i].substr( shiftLen );
            cntChunks = (currentRead.length() + chunkSize - 20) / chunkSize ;
            isSolvedChunk[j] = vector<bool>(cntChunks,false);
            mapResults[j] = vector<mapResult>(cntChunks);

            // [[[Where to check: ]]]
            checkStack.clear();
            ( j == 0 ? checkStack.push_back( interval(1,cntChunks) ) :
                determineCheckStack( &isSolvedChunk, j, &mapResults ) );
            if( checkStack.empty() )
                break;
            // [[[Full Anchoring + Extension: ]]]
            //cerr<<"read: "<<readsName[i]<<endl;
            rdcntrMutex.lock();
            rdcntr++;
            if(readIndex[i] == 171914)
                string a = readsName[i];
            if(rdcntr % 100000 == 0){
                cerr<<readIndex[i];
                cerr<<"__________\n"<<++kcntr<<"00 k reads done\n";
                cerr<<" unMapped:       "<<((double) cnt_unmapped*100/rdcntr)<<"\%\t ("<<cnt_unmapped<<")\n";
                cerr<<" informative:    "<<((double) cnt_informatives*100/rdcntr)<<"\%\t ("<<cnt_informatives<<")\n";
                cerr<<" dbl informative:"<<((double) cnt_2informatives*100/rdcntr)<<"\%\t ("<<cnt_2informatives<<")\n";
            }
            rdcntrMutex.unlock();
            while( !checkStack.empty() ){

                interval solvedRegion =
                        anchorAndExtend( &currentRead, &currentQC, &checkStack, j,
                                         &mapResults, &clusters ,idx );
                if( solvedRegion.start != 0 )
                    for(int k = solvedRegion.start ; k <= solvedRegion.end ; k++)
                        // only needed for determining checkStack in future anchoring iterations
                        isSolvedChunk[j][k] = true;
            }
            if(clusters.size() > 0)
                fullExtension( &currentRead, &clusters, idx );
                                                   // =============================== }
        }
        if( clusters.size() > 1){
            doubInformatives.push_back(informativeReadsData( readIndex[i], readsName[i], reads[i], clusters));
            if( doubInformatives.size() >= WRITEWHEN){
                write2InformativesMutex.lock();
                write2Informatives(&doubInformatives/*, &ofstr_2Informatives*/);
                write2InformativesMutex.unlock();
                doubInformatives.clear();
            }
        }
        if( clusters.size() > 1
                ||
                (clusters.size() == 1
                 && (clusters[0].startChunk != 1
                     || clusters[0].endChunk != ( (int) (reads[i].size() + chunkSize - 20) / chunkSize) ) ) ){
            //isInformative
            informativeReads.push_back(informativeReadsData( readIndex[i], readsName[i], reads[i], clusters));
            if( informativeReads.size() >= WRITEWHEN){
                writeInformativesMutex.lock();
                writeInformatives(&informativeReads/*, &ofstr_Informatives*/);
                writeInformativesMutex.unlock();
                informativeReads.clear();
            }
            //informativeReads
            //add informative reads to a vector to be written into file or for further analyses
            // we can use depth data too
        }else if( clusters.size() == 0 ){
            unMappeds.push_back( unMappedReads(readsName[i],reads[i], qc[i]) );
            if( unMappeds.size() >= WRITEWHEN){
                writeUnMapped(&unMappeds/*, &ofstr_unMapped*/);
                unMappeds.clear();
            }
        }
                                         // =============================== }
    }
    if( doubInformatives.size() > 0){
        write2InformativesMutex.lock();
        write2Informatives(&doubInformatives/*, &ofstr_2Informatives*/);
        write2InformativesMutex.unlock();
        doubInformatives.clear();
    }
    if( informativeReads.size() > 0){
        writeInformativesMutex.lock();
        writeInformatives(&informativeReads/*, &ofstr_Informatives*/);
        writeInformativesMutex.unlock();
        informativeReads.clear();
    }
    if( unMappeds.size() > 0){
        writeUnMapped(&unMappeds/*, &ofstr_unMapped*/);
        unMappeds.clear();
    }
}

void mt_ReadAndProcess(ifstream* fin,int* idt){

    // ////////////////// Variables: (per thread)
    string readsName[MAXREADSIZE + 1000];
    string reads[MAXREADSIZE + 1000];
    long long readsIndex[MAXREADSIZE + 1000];
    string qc[MAXREADSIZE + 1000];
    //resetControllers();
    string line;
    long long cntReads = 0;
    int tt = 0, tmpName = 0, lpNum = 0;

    while( true ){

        // ////////////////// **********************************|
        // ////////////////// STEP 1 ---------------------------|
        // ////////////////// Threads compete to read from file |
        // ////////////////// **********************************|
        readFileMutex.lock();

        // Part 1: reading to this thread buffer here ( min(buffSize,RemainedSize) )

        //cerr<<"\n Thread "<<*idt<<" started to read from file";
        //for( int i = 0; i < seed2; i++ )
        //writing[i] = reading[i] = false;
        cntReads = 0;
        auto t1 = std::chrono::high_resolution_clock::now();
        int mxLen = 0;
        long long totbp = 0;
        //getline( *fin, line );
        //getline( *fin, line );
        while( !(cntReads >= MAXREADSIZE || totbp > 100000000) ){
            // if file ended: open a new file and make the buffer a full buffer
            if( !(getline( *fin, line )) ){
                if(++fileCounter <= readsFileCount){
                    fin->close();   fin->clear();
                    fin->open(readsFileName+"_"+convertNumToStr(fileCounter)+".fq");
                    cerr<<"\n  -- Started a new File to Read..\n";
                }
                else{
                    //readFileMutex.unlock();
                    break;
                }
            }
            //cerr<<"\n thread "<<(*idt)<<" Started to read file";
            if( line.length() == 0 ){
                //readFileMutex.unlock();
                break;
            }
            if( line[0] == '@' && line.length() > 1 ){
                readsName[cntReads] = line.substr( 1 );
                istringstream sin( readsName[cntReads] );
                string tnp;
                sin >> tnp;
                readsName[cntReads] = tnp;
            }
            else{
                ostringstream sout;
                sout << tmpName++;
                readsName[cntReads] = sout.str();
            }
            getline( *fin, line );
            for( int i = 0; i < line.length(); i++ ){
                if( line[i] >= 'A' && line[i] <= 'Z' )
                    line[i] = 'a' + ( line[i] - 'A' );
                if( line[i] != 'a' && line[i] != 'c' && line[i] != 'g' && line[i] != 't' )
                    line[i] = 'a';
            }
            totbp += line.length();
            mxLen = max( mxLen, (int)line.length());
            reads[cntReads] = line;
            readsIndex[cntReads] = readCounter++;
            getline( *fin, line );
            getline( *fin, line );
            //if( totalReads + cntReads >= 750000 )
            //	cerr << reads[cntReads] << endl;
            qc[cntReads++] = line;
            if( cntReads >= MAXREADSIZE || totbp > 100000000){
                //readFileMutex.unlock();
                break;
            }
        }

        if( cntReads == 0 ){
            readFileMutex.unlock();
            break;
            cerr<<"\n no more to read";
        }
        totalReads += cntReads;
        //cerr<<"\n thread "<<(*idt)<<" Finished read file";
        //cerr << "\nProcessing " << cntReads << " reads  (" << totbp << "bp)" << "...\n";
        // Part 2: Last read size and fseek to handle (only when using fread() instead of getline())

        readFileMutex.unlock();

        // ////////////////// **********************************|
        // ////////////////// STEP 2 ---------------------------|
        // ////////////////// Threads seperated processing      |
        // ////////////////// **********************************|

        // Part 1: Parsing buffer into reads (only when using fread() instead of getline())
        {

        }
        // Part 2: Process reads one-by-one
        processing( readsIndex ,readsName ,reads ,qc ,cntReads ,*idt);
        //cerr << "\n Alaki Masalan Processing Reads...\n";
        //        long long matched = 0;
        //        for( int i = 0; i < core; i++ )
        //            matched += mch[i];
        //        tt++;
        auto t2 = std::chrono::high_resolution_clock::now();
        auto dd = std::chrono::duration_cast< std::chrono::seconds>(t2 - t1);
        //        cerr << "###### " << lpNum++ << ": mapping time = " << dd.count() << " seconds ######" << endl;
        //        double nesbat = (double)matched * 100. / (double)totalReads;
        //        cerr << "aligned reads : " << matched << "/" << totalReads << " -- " << nesbat << "%\n\n";
        //        prepareForWriting();
        //break;
    }
}

int main(int argc, char *argv[]){
    QCoreApplication a(argc, argv);
    chdir(MAINADDRESS);
    // ===========================================================================================
    // SureMap main:
    for( int i = 0 ; i < MAXTHREADS ; i++ ){
        minVal[i] = 1000;
        mxFailed[i] = 1000;

        globalRMI[i].uniqeOption = globalUniqeOption;
        globalRMI[i].maxReport = globalMaxReport;
        globalRMI[i].bestOne = globalBestOne;
        globalRMI[i].maxDiffMismatch = globalMaxDiffMismatch;
        globalRMI[i].maxDiffEdit = globalMaxDiffEdit;
        globalRMI[i].gap = globalGap;
        globalRMI[i].noisePercent = globalNoisePercent;
    }
    if( globalMode == "fast" )
        Mode = 4, mc =  4;
    else if( globalMode == "normal" )
        Mode = 3, mc = 10;

    std::srand ( unsigned ( std::time(0) ) );
    vector< unsigned short > test;
    int32_t opt,nqrys,maxqry,i;
    char* idxname;
    char* qryname;
    ///////////////////////////
    // Setting:
    idxname = "/home/ameer/ExactSV/Chr19";
    qryname = "";
    outputAdr = "/outDir";

    FILE* f;
    uint8_t** queries;
    char buf[4096];
    uint32_t start,stop,cnt;
    string args[100];

    // p2

    ofstream fout( outputAdr );
    fout.close();
    string refPrefix = idxname;
    refPrefix = getValidPrefix( refPrefix );
    //cerr << refPrefix << endl;
    string fastqAdr = qryname;
    string fwAdr = refPrefix + ".fm";
    string rvAdr = refPrefix + ".rev.fm";
    string rfInf = refPrefix;
    genomeLength = loadRef( rfInf );
    //cerr<<"\n((("<<genomeLength<<endl;
    cerr << "--Ref loading is done!\n";
    /* load index */
    fmIdx.load(fwAdr);
    cerr << "--forward BWT loading is done!\n";
    revIdx.load( rvAdr );
    cerr << "--Reverse BWT loading is done!\n";
    totChar = fmIdx.sigma;
    //Aligner( fastqAdr );

    // ===========================================================================================


    characterPerRow = findCharacterPerRow(genomeName);
    //genomeLength = getGenomeLength();
    string currentFileName;
    if(readsFileCount == 1)
        currentFileName = readsFileName+".fq";
    else if(readsFileCount > 1)
        currentFileName = readsFileName+"_1.fq";
    else
    {
        cerr<<"\n Wrong number of files... (<1).\n Extited";
        return 0;
    }
    //cerr<<currentFileName<<endl;
    ifstream ifstr(currentFileName.c_str());

    int *ids;
    ids = new int[numberOfThreads];
    vector<thread> threads(numberOfThreads);
    auto T1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numberOfThreads; i++) {
        ids[i] = i;
        //cerr<<"\n before thread "<<i<<endl;
        threads[i] = thread(mt_ReadAndProcess, &ifstr ,  &ids[i] );
    }
    // wait until all threads have finished
    for (int i = 0; i < numberOfThreads; i++){
        threads[i].join();
    }
    auto T2 = std::chrono::high_resolution_clock::now();
    auto DD = std::chrono::duration_cast< std::chrono::seconds>(T2 - T1);
    cerr << "\n________________________________\n Total Run Time = " << DD.count() << "s" << endl;

    cerr<< "\n________________________________\n "<<" Reads count: "<<rdcntr<<endl;
    cerr<<" unMapped:       "<<((double) cnt_unmapped*100/rdcntr)<<"\%\t ("<<cnt_unmapped<<")\n";
    cerr<<" informative:    "<<((double) cnt_informatives*100/rdcntr)<<"\%\t ("<<cnt_informatives<<")\n";
    cerr<<" dbl informative:"<<((double) cnt_2informatives*100/rdcntr)<<"\%\t ("<<cnt_2informatives<<")\n";
    cerr<< "\n________________________________\n";
    cerr<< "\n Successfully Finished.";

    delete ids;
    ifstr.close();
    return a.exec();
}
