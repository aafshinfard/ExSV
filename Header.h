#ifndef HEADER_H
#define HEADER_H

#include <QCoreApplication>
#include <stdio.h>
#include <stdlib.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

// #include <boost/atomic.hpp> // for sureMap
#include <unistd.h>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <stack>
#include <cmath>
#include <iomanip>

using namespace std;

// ================================|
// ================================|
//> Settings:                      |
// ================================|
//                                 |
// ==== Normal Settings:
#define seed1 1000000007
#define seed2  1000019
#define seed3 100000007
#define MAXREADSIZE 10000
#define MAX_SHIFT_ITER 5
#define WRITEWHEN 100
#define MAINADDRESS "/home/ameer/ExactSV/"
//#define BUFFERSIZE 1024*1024    // Size of fread Buffer
//#define MAX_BUFFER_SIZE 100     // each thread will read MAX_BUFFER_SIZE reads in each trial
// ====

int numberOfThreads = 8;
// ====
int chunkSize = 25;
// ==== Files:
string readsFileName = "reads3";    //string readName = "reads_k12.fq";
int readsFileCount = 1;
string genomeName = "Ref_SV.fa"; //string genomeName = "E_coli_K12_DH10B.fa";    //string genomeName = "chr19.fa";
string indexName = "Ref_SV"; //string indexName = "E_coli_K12_DH10B";        //string indexName = "chr19";
string outputDir = "/home/ameer/SVAnchoring/SV_out2/";

string informativeFileName = "informatives";
//string outputDir = "/home/ubuntu/SVAnchoring/SV_out2/";

// ==== Alignment Parameters:
int dashV = 1;
double gapRate = 5; // rate (#) of acceptable gap caused by mutation or sequencing error per chunk (for correspondance checking)
double indelShift = 0.4;
int editDistance = 2;
int gapPen = -8, misPen = -5, matchPen = 10;
//int numGap = 5;     //numGap in Anchoring | error in Assignment
int LocAlth = 130;
//int anchoringShift = 5;
//int assignmentShift = 10; //d = chunkSize - assignmentShift
bool runAnchoring = true;
bool BGLR = true;               // [true]: Left and Right Binary Genomes
                                // [false]: only 1 BinaryGenome (less accurate and slower But less Ram usage)
int shiftIterations = 2;        // max = MAX_SHIFT_ITER
bool logTxt = false;
bool pairedEndedReads = false;  // if true, then reads must be in two files "readName"_1.fq and "readName"_2.fq
bool runPhase2 = false;





long long cntCOM = 0, cntBG=0;
//                                 |
// ================================|
//> End of Settings:               |
// ================================|
// ================================'




#endif // HEADER_H
