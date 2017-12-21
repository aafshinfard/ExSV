/**
 * @file    DataStructures.h
 * @author  Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                                  <afshinfard (at) gmail.com>
 *
 * @section LICENCE
 *
 * Copyright (C) 2017-2020
 *   Amirhosein Afshinfard   <afshinfard (at) ce.sharif.edu>
 *                           <a.afshinfard (at) gmail.com>
 *   Damoon Nashta-ali       <damoun_dna (at) yahoo.com>
 *	 Seyed Abolfazl Motahari <motahari (at) sharif.edu
 *
 **/

#include "LocalAligner.h"
using namespace std;

LocalAligner::LocalAligner(std::string seqA, std::string seqB, int gapPen, int misPen, int matchPen, int shift, int k)
{
    this->MATCH = matchPen;
    this->GAP = gapPen;
    this->counter = 0;
    this->MISMATCH = misPen;
    this->shift = shift;
    this->NegativeInf=-1000 * MATCH;
    this->k=k+2*shift;
    for (int i=0;i < seqA.size(); i++)
        seqA[i] = toupper(seqA[i]);
    for (int i=0;i < seqB.size(); i++)
        seqB[i] = toupper(seqB[i]);
    this->mSeqA = seqA;
    this->mSeqB = seqB;
    this->mD = new int*[this->mSeqA.size() + 1];
    for (int i = 0; i <= this->mSeqA.size(); i++) {
        this->mD[i] = new int[this->mSeqB.size() + 1];
    }
    for (int i = 0; i <= this->mSeqA.size(); i++) {
        for (int j = 0; j <= this->mSeqB.size(); j++) {
            if (!i) this->mD[i][j] = 0;
            else if (j == 0) {
                mD[i][j] = 0;
            } else {
                mD[i][j] = 0;
            }
        }
    }
}

LocalAligner::~LocalAligner()
{
    for (int i = 0; i <= this->mSeqA.size(); i++)
        delete[] this->mD[i];
    delete this->mD;
}

void LocalAligner::process()
{
    int scoreDiag = 0;
    int scoreLeft = 0;
    int scoreUp = 0;
    int end=this->mSeqB.size();
    for (int i = 1; i <= this->mSeqA.size(); i++)
    {
        if((((this->k/2)+i)<(this->mSeqB.size())))
            end = (this->k/2)+i;
        for (int j = std::max(1, i-(this->k/2)); j <= end; j++)
        {
            if(i < this->mSeqA.size()-shift+1)
            {
                scoreDiag = mD[i-1][j-1] + weight(i, j);
                scoreLeft = mD[i][j-1] + GAP;
                scoreUp = mD[i-1][j] + GAP;
                mD[i][j] = std::max(std::max(scoreDiag, scoreLeft), scoreUp);
            }
            else
            {
                scoreDiag = mD[i-1][j-1] + MISMATCH;
                scoreLeft = mD[i][j-1] + GAP - 10* MATCH;
                scoreUp = mD[i-1][j] + GAP + 10* MATCH;
                mD[i][j] = std::max(std::max(scoreDiag, scoreLeft), scoreUp);
            }
        }
        if(1<(i-(this->k/2)))
        {
            mD[i][i-(this->k/2)-1]=this->NegativeInf;
        }
        if(((this->k/2)+i)<(this->mSeqB.size())){
            mD[i][(this->k/2)+i+1]=this->NegativeInf;
        }
    }

}
void LocalAligner::Print(){
    for(int i=0;i<this->mSeqA.size();i++){
        for(int j=0;j<this->mSeqB.size();j++){
            cout << mD[i][j]<<"\t";
        }
        cout << endl;
    }
}

void LocalAligner::backtrack(){
    //size_t i = this->mSeqA.size() - shift;
    size_t i = this->mSeqA.size() - shift;
    int max = -1000 * MATCH , index = 0, flag = 0;
    for (int k = shift; k <= this->mSeqB.size(); k++)
    {
        if(mD[i][k] > max)
        {
            max = mD[i][k];
            index = k;
        }
    }
    size_t j = index;
    this->mScore = mD[i][j];
    this->cigar = "";
    this->realCigar = "";
    std::ostringstream os;

    while (i > 0 && j > 0) {

        if (mD[i][j] == mD[i-1][j] + GAP) {
            mAlignmentSeqA += this->mSeqA[i-1];
            mAlignmentSeqB += "-";
            totalGapfa ++;
            i--;
            continue;
        }
        else if (mD[i][j] == mD[i-1][j-1] + weight(i, j)) { // MATCH-MISMATCH
            mAlignmentSeqA += this->mSeqA[i-1];
            mAlignmentSeqB += this->mSeqB[j-1];
            if(weight(i,j)==MISMATCH)
                totalMismatch ++;
            i--;
            j--;
            continue;
        }
        else
        {
            mAlignmentSeqA += "-";
            totalGapr ++;
            mAlignmentSeqB += this->mSeqB[j-1];
            j--;
            continue;
        }
    }
    std::reverse(this->mAlignmentSeqA.begin(), mAlignmentSeqA.end());
    std::reverse(this->mAlignmentSeqB.begin(), mAlignmentSeqB.end());
    //std::cout << mAlignmentSeqA << "\t" << mAlignmentSeqB << std::endl;
    produceCigar();
}

void LocalAligner::produceCigar(){
    std::string temp;
    temp = this->mAlignmentSeqA;
    this->mAlignmentSeqA = this->mAlignmentSeqB;
    this->mAlignmentSeqB = temp;
    int counter = 0, realCounter = 0;
    cigar = "";
    realCigar = "";
    int gapCount= 0, Icounter = 0, Dcounter = 0;
    std::ostringstream os;
    for (int i = 0 ; i < this->mAlignmentSeqA.size(); i++) {
        if (mAlignmentSeqA[i] == mAlignmentSeqB[i]){
            gapCount = 0;
            counter++;
            if (Icounter) {
                os.str("");
                os<<Icounter;
                realCigar += os.str();
                realCigar += "I";
            }
            if (Dcounter) {
                os.str("");
                os<<Dcounter;
                realCigar += os.str();
                realCigar += "D";
            }
            Icounter = 0;
            Dcounter = 0;
            realCounter++;
        }
        else if (mAlignmentSeqA[i] != mAlignmentSeqB[i] && mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] != '-'){
            gapCount = 0;
            if (Icounter) {
                os.str("");
                os<<Icounter;
                realCigar += os.str();
                realCigar += "I";
            }
            if (Dcounter) {
                os.str("");
                os<<Dcounter;
                realCigar += os.str();
                realCigar += "D";
            }
            if (counter) {
                os.str("");
                os<<counter;
                cigar += os.str();
            }
            Icounter = 0;
            Dcounter = 0;
            counter = 0;
            cigar += mAlignmentSeqA[i];
            realCounter++;
        }
        else if (mAlignmentSeqA[i] == '-' && mAlignmentSeqB[i] != '-') {
            gapCount ++;
            Icounter ++;
            if (counter) {
                os.str("");
                os<<counter;
                cigar += os.str();
            }
            if (gapCount == 1) cigar += "^";
            cigar += mAlignmentSeqB[i];
            counter = 0;

            if (realCounter) {
                os.str("");
                os<<realCounter;
                realCigar += os.str();
                realCigar += "M";
            }
            realCounter = 0;

            if (Dcounter) {
                os.str("");
                os<<Dcounter;
                realCigar += os.str();
                realCigar += "D";
            }
            Dcounter = 0;
        }
        else if (mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] == '-'){
            gapCount = 0;
            Dcounter ++;
            if (realCounter) {
                os.str("");
                os<<realCounter;
                realCigar += os.str();
                realCigar += "M";
            }
            realCounter = 0;
            if (Icounter) {
                os.str("");
                os<<Icounter;
                realCigar += os.str();
                realCigar += "I";
            }
            Icounter = 0;
        }

    }
    if (counter) {
        os.str("");
        os<<counter;
        cigar += os.str();
    }

    if (realCounter) {
        os.str("");
        os<<realCounter;
        realCigar += os.str();
        realCigar += "M";
    }
    if (Dcounter) {
                os.str("");
                os<<Dcounter;
                realCigar += os.str();
                realCigar += "D";
    }
    if (Icounter) {
                os.str("");
                os<<Icounter;
                realCigar += os.str();
                realCigar += "I";
    }
}


int LocalAligner::weight(size_t i, size_t j) {
    ++counter;
    if (this->mSeqA[i - 1] == this->mSeqB[j - 1])
        return MATCH;
    else
        return MISMATCH;
}
////
////  LocalAligner.cpp
////  MetaAligner
////
////  Created by ahmad on 12/19/93.
////  Copyright (c) 1393 ahmad. All rights reserved.
////

//#include "LocalAligner.h"

//using namespace std;

//LocalAligner::LocalAligner(std::string seqA, std::string seqB, int gapPen, int misPen, int matchPen, int shift, int k) {
//    this->MATCH = matchPen;
//    this->GAP = gapPen;
//    this->counter = 0;
//    this->MISMATCH = misPen;
//    this->shift = shift;
//    this->NegativeInf=-1000;
//    this->k=k+2*shift;
//    for (int i=0;i < seqA.size(); i++){
//        seqA[i] = toupper(seqA[i]);
//        seqB[i] = toupper(seqB[i]);
//    }
//    this->mSeqA = seqA;
//    this->mSeqB = seqB;
//    this->mD = new int*[this->mSeqA.size() + 1];
//    for (int i = 0; i <= this->mSeqA.size(); i++) {
//        this->mD[i] = new int[this->mSeqB.size() + 1];
//    }
//    for (int i = 0; i <= this->mSeqA.size(); i++) {
//        for (int j = 0; j <= this->mSeqB.size(); j++) {
//            if (!i) this->mD[i][j] = 0;
//            else if (j == 0) {
//                mD[i][j] = 0;
//            } else {
//                mD[i][j] = 0;
//            }
//        }
//    }
//}

//LocalAligner::~LocalAligner() {
//    for (int i = 0; i <= this->mSeqA.size(); i++)
//        delete[] this->mD[i];
//    delete this->mD;
//}

//void LocalAligner::process() {
//    int scoreDiag = 0;
//    int scoreLeft = 0;
//    int scoreUp = 0;
//    int end=this->mSeqB.size();
//    for (int i = 1; i <= this->mSeqA.size(); i++) {
//        if((((this->k/2)+i)<(this->mSeqB.size())))
//            end = (this->k/2)+i;
//        for (int j = std::max(1, i-(this->k/2)); j <= end; j++) {
//            if(j < this->mSeqB.size()-shift+1){
//                scoreDiag = mD[i-1][j-1] + weight(i, j);
//                scoreLeft = mD[i][j-1] + GAP;
//                scoreUp = mD[i-1][j] + GAP;
//                mD[i][j] = std::max(std::max(scoreDiag, scoreLeft), scoreUp);
//            }
//            else{
//                scoreDiag = mD[i-1][j-1] + weight(i, j);
//                scoreLeft = mD[i][j-1] + GAP - 10;
//                scoreUp = mD[i-1][j] + GAP + 10;
//                mD[i][j] = std::max(std::max(scoreDiag, scoreLeft), scoreUp);
//            }
//        }
//        if(1<(i-(this->k/2))){
//            mD[i][i-(this->k/2)-1]=this->NegativeInf;
//        }
//        if(((this->k/2)+i)<(this->mSeqB.size())){
//            mD[i][(this->k/2)+i+1]=this->NegativeInf;
//        }
//    }
//}

//void LocalAligner::Print() {
//    for(int i=0;i<this->mSeqA.size();i++){
//        for(int j=0;j<this->mSeqB.size();j++){
//            cout << mD[i][j]<<"\t";
//        }
//        cout << endl;
//    }
//}

//void LocalAligner::backtrack() {
//    size_t i = this->mSeqA.size()-shift;
//    size_t j = this->mSeqB.size()-shift;
//    this->mScore = mD[i][j];
//    this->cigar = "";
//    this->realCigar = "";
//    std::ostringstream os;

//    while (i > 0 && j > 0) {

//        if (mD[i][j] == mD[i-1][j-1] + weight(i, j)) { // MATCH-MISMATCH
//            mAlignmentSeqA += mSeqA[i-1];
//            mAlignmentSeqB += mSeqB[j-1];
//            i--;
//            j--;
//            continue;
//        } else if (mD[i][j] == mD[i][j-1] + GAP) {
//            mAlignmentSeqA += "-";
//            mAlignmentSeqB += mSeqB[j-1];
//            j--;
//            continue;
//        } else {
//            mAlignmentSeqA += mSeqA[i-1];
//            mAlignmentSeqB += "-";
//            i--;
//            continue;
//        }
//    }
//    std::reverse(this->mAlignmentSeqA.begin(), mAlignmentSeqA.end());
//    std::reverse(this->mAlignmentSeqB.begin(), mAlignmentSeqB.end());
//    //std::cout << mAlignmentSeqA << "\t" << mAlignmentSeqB << std::endl;
//    produceCigar();
//}

//void LocalAligner::produceCigar() {
//    std::string temp;
//    temp = this->mAlignmentSeqA;
//    this->mAlignmentSeqA = this->mAlignmentSeqB;
//    this->mAlignmentSeqB = temp;
//    int counter = 0, realCounter = 0;
//    cigar = "";
//    realCigar = "";
//    int gapCount= 0, Icounter = 0, Dcounter = 0;
//    std::ostringstream os;
//    for (int i = 0 ; i < this->mAlignmentSeqA.size(); i++) {
//        if (mAlignmentSeqA[i] == mAlignmentSeqB[i]){
//            gapCount = 0;
//            counter++;
//            if (Icounter) {
//                os.str("");
//                os<<Icounter;
//                realCigar += os.str();
//                realCigar += "I";
//            }
//            if (Dcounter) {
//                os.str("");
//                os<<Dcounter;
//                realCigar += os.str();
//                realCigar += "D";
//            }
//            Icounter = 0;
//            Dcounter = 0;
//            realCounter++;
//        }
//        else if (mAlignmentSeqA[i] != mAlignmentSeqB[i] && mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] != '-'){
//            totalMismatch ++;
//            gapCount = 0;
//            if (Icounter) {
//                os.str("");
//                os<<Icounter;
//                realCigar += os.str();
//                realCigar += "I";
//            }
//            if (Dcounter) {
//                os.str("");
//                os<<Dcounter;
//                realCigar += os.str();
//                realCigar += "D";
//            }
//            if (counter) {
//                os.str("");
//                os<<counter;
//                cigar += os.str();
//            }
//            Icounter = 0;
//            Dcounter = 0;
//            counter = 0;
//            cigar += mAlignmentSeqA[i];
//            realCounter++;
//        }
//        else if (mAlignmentSeqA[i] == '-' && mAlignmentSeqB[i] != '-') {
//            totalGap ++;
//            gapCount ++;
//            Icounter ++;
//            if (counter) {
//                os.str("");
//                os<<counter;
//                cigar += os.str();
//            }
//            if (gapCount == 1) cigar += "^";
//            cigar += mAlignmentSeqB[i];
//            counter = 0;

//            if (realCounter) {
//                os.str("");
//                os<<realCounter;
//                realCigar += os.str();
//                realCigar += "M";
//            }
//            realCounter = 0;

//            if (Dcounter) {
//                os.str("");
//                os<<Dcounter;
//                realCigar += os.str();
//                realCigar += "D";
//            }
//            Dcounter = 0;
//        }
//        else if (mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] == '-'){
//            totalGap ++;
//            gapCount = 0;
//            Dcounter ++;
//            if (realCounter) {
//                os.str("");
//                os<<realCounter;
//                realCigar += os.str();
//                realCigar += "M";
//            }
//            realCounter = 0;
//            if (Icounter) {
//                os.str("");
//                os<<Icounter;
//                realCigar += os.str();
//                realCigar += "I";
//            }
//            Icounter = 0;
//        }

//    }
//    if (counter) {
//        os.str("");
//        os<<counter;
//        cigar += os.str();
//    }

//    if (realCounter) {
//        os.str("");
//        os<<realCounter;
//        realCigar += os.str();
//        realCigar += "M";
//    }
//    if (Dcounter) {
//        os.str("");
//        os<<Dcounter;
//        realCigar += os.str();
//        realCigar += "D";
//    }
//    if (Icounter) {
//        os.str("");
//        os<<Icounter;
//        realCigar += os.str();
//        realCigar += "I";
//    }
//}

//int LocalAligner::weight(size_t i, size_t j) {
//    ++counter;
//    if (this->mSeqA[i - 1] == this->mSeqB[j - 1])
//        return MATCH;
//    else
//        return MISMATCH;
//}
//// V2 :
////#include "LocalAligner.h"
////using namespace std;
////LocalAligner::LocalAligner(std::string seqA, std::string seqB, int gapPen, int misPen, int matchPen, int shift, int k){
////    this->MATCH = matchPen;
////    this->GAP = gapPen;
////    this->counter = 0;
////    this->MISMATCH = misPen;
////    this->shift = shift;
////    this->NegativeInf=-1000;
////    this->k=k+2*shift;
////    for (int i=0;i < seqA.size(); i++){
////        seqA[i] = toupper(seqA[i]);
////        seqB[i] = toupper(seqB[i]);
////    }
////    this->mSeqA = seqA;
////    this->mSeqB = seqB;
////    this->mD = new int*[this->mSeqA.size() + 1];
////    for (int i = 0; i <= this->mSeqA.size(); i++) {
////        this->mD[i] = new int[this->mSeqB.size() + 1];
////    }
////    for (int i = 0; i <= this->mSeqA.size(); i++) {
////        for (int j = 0; j <= this->mSeqB.size(); j++) {
////            if (!i) this->mD[i][j] = 0;
////            else if (j == 0) {
////                mD[i][j] = 0;
////            } else {
////                mD[i][j] = 0;
////            }
////        }
////    }
////}

/////*LocalAligner::~LocalAligner() {
////    for (int i = 0; i <= this->mSeqA.size(); i++)
////        delete[] this->mD[i];
////    delete this->mD;
////}*/

////void LocalAligner::process() {
////    int scoreDiag = 0;
////    int scoreLeft = 0;
////    int scoreUp = 0;
////    int end=this->mSeqB.size();
////    for (int i = 1; i <= this->mSeqA.size(); i++) {
////        if((((this->k/2)+i)<(this->mSeqB.size())))
////            end = (this->k/2)+i;
////        for (int j = std::max(1, i-(this->k/2)); j <= end; j++) {
////            if(j < this->mSeqB.size()-shift+1)
////            {
////                scoreDiag = mD[i-1][j-1] + weight(i, j);
////                scoreLeft = mD[i][j-1] + GAP;
////                scoreUp = mD[i-1][j] + GAP;
////                mD[i][j] = std::max(std::max(scoreDiag, scoreLeft), scoreUp);
////            }
////            else
////            {
////                scoreDiag = mD[i-1][j-1] + weight(i, j);
////                scoreLeft = mD[i][j-1] + GAP - 10 * MATCH;
////                scoreUp = mD[i-1][j] + GAP - 10 * MATCH;
////                mD[i][j] = std::max(std::max(scoreDiag, scoreLeft), scoreUp);
////            }
////        }
////        if(1<(i-(this->k/2))){
////            mD[i][i-(this->k/2)-1]=this->NegativeInf;
////        }
////        if(((this->k/2)+i)<(this->mSeqB.size())){
////            mD[i][(this->k/2)+i+1]=this->NegativeInf;
////        }
////    }
////}
////void LocalAligner::Print(){
////    for(int i=0;i<this->mSeqA.size();i++){
////        for(int j=0;j<this->mSeqB.size();j++){
////            cout << mD[i][j]<<"\t";
////        }
////        cout << endl;
////    }
////}

////void LocalAligner::backtrack(){
////    //size_t i = this->mSeqA.size() - shift;
////    size_t j = this->mSeqB.size() - shift;
////    int max = 0 , index = 0, flag = 0;
////    for (int k = this->mSeqA.size() - shift; k < this->mSeqA.size(); k++)
////    {
////        if(mD[k][j] > max)
////        {
////            max = mD[k][j];
////            index = k;
////        }
////    }
////    size_t i = index;
////    this->cigar = "";
////    this->realCigar = "";
////    std::ostringstream os;


////    while (i > 0 && j > 0) {

////        if (mD[i][j] == mD[i][j-1] + GAP) {
////                    mAlignmentSeqA += "-";
////                    mAlignmentSeqB += mSeqB[j-1];
////                    j--;
////                    continue;
////        }
////        else if (mD[i][j] == mD[i-1][j-1] + weight(i, j)) { // MATCH-MISMATCH
////            mAlignmentSeqA += mSeqA[i-1];
////            mAlignmentSeqB += mSeqB[j-1];
////            i--;
////            j--;
////            continue;
////        }
////        else {
////            mAlignmentSeqA += mSeqA[i-1];
////            mAlignmentSeqB += "-";
////            i--;
////            continue;
////        }
////    }
////    std::reverse(this->mAlignmentSeqA.begin(), mAlignmentSeqA.end());
////    std::reverse(this->mAlignmentSeqB.begin(), mAlignmentSeqB.end());
////    //std::cout << mAlignmentSeqA << "\t" << mAlignmentSeqB << std::endl;
////    produceCigar();
////}


////void LocalAligner::produceCigar(){
////    std::string temp;
////    temp = this->mAlignmentSeqA;
////    this->mAlignmentSeqA = this->mAlignmentSeqB;
////    this->mAlignmentSeqB = temp;
////    int counter = 0, realCounter = 0;
////    cigar = "";
////    realCigar = "";
////    int gapCount= 0, Icounter = 0, Dcounter = 0;
////    std::ostringstream os;
////    for (int i = 0 ; i < this->mAlignmentSeqA.size(); i++) {
////        if (mAlignmentSeqA[i] == mAlignmentSeqB[i]){
////            gapCount = 0;
////            counter++;
////            if (Icounter) {
////                os.str("");
////                os<<Icounter;
////                totalGapr += Icounter;
////                realCigar += os.str();
////                realCigar += "I";
////            }
////            if (Dcounter) {
////                os.str("");
////                os<<Dcounter;
////                totalGapfa += Dcounter;
////                realCigar += os.str();
////                realCigar += "D";
////            }
////            Icounter = 0;
////            Dcounter = 0;
////            realCounter++;
////        }
////        else if (mAlignmentSeqA[i] != mAlignmentSeqB[i] && mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] != '-'){
////            gapCount = 0;
////            totalMismatch ++;
////            if (Icounter) {
////                os.str("");
////                os<<Icounter;
////                totalGapr += Icounter;
////                realCigar += os.str();
////                realCigar += "I";
////            }
////            if (Dcounter) {
////                os.str("");
////                os<<Dcounter;
////                totalGapfa += Dcounter;
////                realCigar += os.str();
////                realCigar += "D";
////            }
////            if (counter) {
////                os.str("");
////                os<<counter;
////                cigar += os.str();
////            }
////            Icounter = 0;
////            Dcounter = 0;
////            counter = 0;
////            cigar += mAlignmentSeqA[i];
////            realCounter++;
////        }
////        else if (mAlignmentSeqA[i] == '-' && mAlignmentSeqB[i] != '-') {
////            gapCount ++;
////            Icounter ++;
////            if (counter) {
////                os.str("");
////                os<<counter;
////                cigar += os.str();
////            }
////            if (gapCount == 1) cigar += "^";
////            cigar += mAlignmentSeqB[i];
////            counter = 0;

////            if (realCounter) {
////                os.str("");
////                os<<realCounter;
////                realCigar += os.str();
////                realCigar += "M";
////            }
////            realCounter = 0;

////            if (Dcounter) {
////                os.str("");
////                os<<Dcounter;
////                totalGapfa += Dcounter;
////                realCigar += os.str();
////                realCigar += "D";
////            }
////            Dcounter = 0;
////        }
////        else if (mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] == '-'){
////            gapCount = 0;
////            Dcounter ++;
////            if (realCounter) {
////                os.str("");
////                os<<realCounter;
////                realCigar += os.str();
////                realCigar += "M";
////            }
////            realCounter = 0;
////            if (Icounter) {
////                os.str("");
////                os<<Icounter;
////                totalGapr += Icounter;
////                realCigar += os.str();
////                realCigar += "I";
////            }
////            Icounter = 0;
////        }

////    }
////    if (counter) {
////        os.str("");
////        os<<counter;
////        cigar += os.str();
////    }

////    if (realCounter) {
////        os.str("");
////        os<<realCounter;
////        realCigar += os.str();
////        realCigar += "M";
////    }
////    if (Dcounter) {
////        os.str("");
////        os<<Dcounter;
////        totalGapfa += Dcounter;
////        realCigar += os.str();
////        realCigar += "D";
////    }
////    if (Icounter) {
////        os.str("");
////        os<<Icounter;
////        totalGapr += Icounter;
////        realCigar += os.str();
////        realCigar += "I";
////    }
////}


////int LocalAligner::weight(size_t i, size_t j) {
////    ++counter;
////    if (this->mSeqA[i - 1] == this->mSeqB[j - 1])
////        return MATCH;
////    else
////        return MISMATCH;
////}
