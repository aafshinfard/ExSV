#ifndef LocalAligner_h
#define LocalAligner_h

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

class LocalAligner {

public:
    std::string mAlignmentSeqA;
    std::string mAlignmentSeqB;
    LocalAligner(std::string, std::string, int, int, int, int, int);
    void process();
    void backtrack();
    void Print();
    ~LocalAligner();
    int mScore, totalMismatch = 0, totalGap = 0, totalGapr = 0, totalGapfa = 0; //totalMatch = 0
    std::string cigar;
    std::string realCigar;
    void produceCigar();
    int MATCH;
    int MISMATCH;
    int GAP;
    int counter;
    int shift;
    int NegativeInf;
    int k;

private:
    std::string mSeqA;
    std::string mSeqB;
    int** mD;
    int weight(size_t, size_t);
    std::string reverseCigar();
};


#endif // LOCALALIGNER_H

//#ifndef LocalAligner_h
//#define LocalAligner_h

//#include <iostream>
//#include <string>
//#include <vector>
//#include <sstream>
//#include <algorithm>

//class LocalAligner {

//public:
//    std::string mAlignmentSeqA;
//    std::string mAlignmentSeqB;
//    LocalAligner(std::string, std::string, int, int, int, int, int);
//    void process();
//    void backtrack();
//    void Print();
//    ~LocalAligner();
//    int mScore, totalMismatch = 0, totalGap = 0; //totalMatch = 0
//    std::string cigar;
//    std::string realCigar;
//    void produceCigar();
//    int MATCH;
//    int MISMATCH;
//    int GAP;
//    int counter;
//    int shift;
//    int NegativeInf;
//    int k;

//private:
//    std::string mSeqA;
//    std::string mSeqB;
//    int** mD;
//    int weight(size_t, size_t);
//    std::string reverseCigar();
//};


//#endif // LOCALALIGNER_H

//// V2 :
////#ifndef LocalAligner_h
////#define LocalAligner_h

////#include <iostream>
////#include <string>
////#include <vector>
////#include <sstream>
////#include <algorithm>

////class LocalAligner {

////public:
////    std::string mAlignmentSeqA;
////    std::string mAlignmentSeqB;
////    LocalAligner(std::string, std::string, int, int, int, int, int);
////    void process();
////    void backtrack();
////    void Print();
////    ~LocalAligner();
////    int mScore, totalMismatch = 0, totalGap = 0, totalGapr = 0, totalGapfa = 0; //totalMatch = 0
////    std::string cigar;
////    std::string realCigar;
////    void produceCigar();
////    int MATCH;
////    int MISMATCH;
////    int GAP;
////    int counter;
////    int shift;
////    int NegativeInf;
////    int k;

////private:
////    std::string mSeqA;
////    std::string mSeqB;
////    int** mD;
////    int weight(size_t, size_t);
////    std::string reverseCigar();
////};


////#endif // LOCALALIGNER_H
