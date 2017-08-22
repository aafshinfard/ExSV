#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H


struct mapResult{
    long long position = 0; // 0:not searched yet - -1: can't find uniqely - -2: can't fnd with least size
    bool isReverse = false;
    //int flag;
    //string cigar;
    mapResult(){}
    mapResult(long long p, bool isReverse):position(p),isReverse(isReverse){}
};

struct solvedFragInfo{

    long startChunk;
    long long startPos = 0;

    long endChunk;
    long long endPos = 0;

    bool isReverse = false;
    int shift = -1;

    //int flag;
    //string cigar;
    solvedFragInfo(){}
    solvedFragInfo(long start_c, long long start_p, long end_c, long long end_p, bool rev, int shiftt):
        startChunk(start_c),startPos(start_p),endChunk(end_c),endPos(end_p),isReverse(rev),shift(shiftt){}
};

struct informativeReadsData{
    long long index;
    string readName;
    string readSeq;
    vector<solvedFragInfo> clusters;
    informativeReadsData(){}
    informativeReadsData( long long indexx, string readname, string readseq, vector<solvedFragInfo> clusterss):
        index(indexx),readName(readname),readSeq(readseq){
        clusters = clusterss;
    }

};

struct interval{
    long start;
    long end;
    interval(){}
    interval(long p, long q):start(p),end(q){}

};

struct singleBowtied{
    long long readIndex;
    int fragNum;
    long long snap;
    int flag;
    singleBowtied(){}
    singleBowtied(long long rI, int fN, long long s, int f):readIndex(rI),fragNum(fN),snap(s),flag(f){}
};
// usage : singleUniqes[ stoll(tokens[0])-1 ].push_back( singleBowtied(stoll(tokens[0]) , stoi(tokens[1]) , stoll(tokens[2]) , stoi(tokens[3]) ) );

#endif // DATASTRUCTURES_H
