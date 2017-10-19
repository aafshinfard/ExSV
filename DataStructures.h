/**
 * @file    DataStructures.h
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


#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H
#include <cstdint>

struct informativeChunk{
    string readName = "";

    uint32_t loci = 0;
    bool bpIsRight = 0;

    bool isRelated = 0;

    uint32_t related_loci = 0;
    bool related_bpIsRight = 0;

    informativeChunk(){}
    informativeChunk( const informativeChunk &a ):
    readName(a.readName), loci(a.loci), bpIsRight(a.bpIsRight), isRelated(a.isRelated), related_loci(a.related_loci), related_bpIsRight(a.related_bpIsRight){}
    informativeChunk( string readNamee, uint32_t locii, bool bpIsRightt, bool isRelatedd, uint32_t related_locii, bool related_bpIsRightt):
    readName(readNamee), loci(locii), bpIsRight(bpIsRightt), isRelated(isRelatedd), related_loci(related_locii), related_bpIsRight(related_bpIsRightt){}
    //    bool operator<(const informativeChunk &a, const informativeChunk &b)
    //    {
    //        return a.loci < b.loci;
    //    }
};

//std::ostream& dump(std::ostream &o, const informativeChunk& p)
//{
//    return o << p.readName << endl << p.loci << endl << p.bpIsRight << endl << p.isRelated << endl << p.related_loci << endl << p.related_bpIsRight << endl;
//}
//std::ostream& operator << (std::ostream &o,const informativeChunk &a){
//  return dump(o,a);
//}
std::ostream& operator<<(std::ostream& os, const informativeChunk& pr)
{
    return os << pr.readName << " " << pr.loci << " " << pr.bpIsRight << " " << pr.isRelated << " " << pr.related_loci << " " << pr.related_bpIsRight   << endl;
}
std::istream& operator>>(std::istream& is, informativeChunk& pr)
{
    string newName;
    informativeChunk new_pr;
    if(    is >> std::ws
        && is >> newName /*std::getline(is, newName, 'q')*/
        && is >> new_pr.loci >> new_pr.bpIsRight >> new_pr.isRelated >> new_pr.related_loci >> new_pr.related_bpIsRight )
    {
        new_pr.readName = newName;
        pr = new_pr; // could do more validation here
    }else cout<<"Reading format problem!";
    return is;
}

bool operator<(const informativeChunk &a, const informativeChunk &b)
{
    return a.loci < b.loci;
}
bool compLoci(const informativeChunk &a, const informativeChunk &b){
    return a.loci < a.loci;
}
bool compConnected(const informativeChunk &a, const informativeChunk &b){
    return a.related_loci < a.related_loci;
}

struct connectedEvent{
    long long index = 0;
    bool isRightBP = false;
    int weight = 0;
    connectedEvent(){}
    connectedEvent( long long indexx, bool isRightBPp, int weightt ):
    index(indexx), isRightBP(isRightBPp), weight(weightt){}
};
bool operator<(const connectedEvent &a, const connectedEvent &b)
{
    return a.weight < b.weight;
}

struct feasibleEvents{
    vector<informativeChunk> informatives;
    int maxDepth = 0;
    uint32_t start = 0;
    uint32_t end = 0;
    vector<connectedEvent> connectedEvents;
};

struct mapResult{
    long long position = 0; // 0:not searched yet - -1: can't find uniqely - -2: can't fnd with least size
    bool isReverse = false;
    //int flag;
    //string cigar;
    mapResult(){}
    mapResult(long long p, bool isReverse):position(p),isReverse(isReverse){}
};

struct solvedFragInfo{

    long startChunk = 0;
    long long startPos = 0;

    long endChunk = 0;
    long long endPos = 0;

    bool isReverse = false;
    int shift = -1;

    //int flag;
    //string cigar;
    solvedFragInfo(){}
    solvedFragInfo(long start_c, long long start_p, long end_c, long long end_p, bool rev, int shiftt):
        startChunk(start_c),startPos(start_p),endChunk(end_c),endPos(end_p),isReverse(rev),shift(shiftt){}
//    bool operator<(const solvedFragInfo &a, const solvedFragInfo &b)
//    {
//        return a.startChunk < b.startChunk;
//    }
};
bool operator<(const solvedFragInfo &a, const solvedFragInfo &b)
{
    return a.startChunk < b.startChunk;
}

struct unMappedReads{
    string name = "";
    string read = "";
    string quality = "";
    unMappedReads(){}
    unMappedReads(string namee, string readd, string qualityy ):
        name(namee),read(readd),quality(qualityy){}
};


//
//
//
//
//
struct informativeReadsData{
    long long index = 0;
    string readName = "EMP";
    string readSeq = "EMP";
    vector<solvedFragInfo> clusters;
    informativeReadsData(){}
    informativeReadsData( long long indexx, string readname, string readseq, vector<solvedFragInfo> clusterss):
        index(indexx),readName(readname),readSeq(readseq){
        clusters = clusterss;
    }

};
//std::istream& operator>>(std::istream& is, informativeReadsData& pr)
//{
//    informativeReadsData new_pr;
//    if(    is >> std::ws
//        && std::getline(is, new_pr.Name, '~')
//        && is >> new_pr.StartSalary >> new_pr.Age )
//    {
//        pr = new_pr; // could do more validation here
//    }
//    return is;
//}

//std::ostream& operator<<(std::ostream& os, const informativeReadsData& pr)
//{
//    return os << pr.index << " " << pr.readName << " " << pr.readSeq << " " << pr.clusters << '\n';
//}

struct interval{
    long start = 0;
    long end = 0;
    interval(){}
    interval(long p, long q):start(p),end(q){}

};

struct singleBowtied{
    long long readIndex = 0;
    int fragNum = 0;
    long long snap = 0;
    int flag = 0;
    singleBowtied(){}
    singleBowtied(long long rI, int fN, long long s, int f):readIndex(rI),fragNum(fN),snap(s),flag(f){}
};
// usage : singleUniqes[ stoll(tokens[0])-1 ].push_back( singleBowtied(stoll(tokens[0]) , stoi(tokens[1]) , stoll(tokens[2]) , stoi(tokens[3]) ) );

#endif // DATASTRUCTURES_H
