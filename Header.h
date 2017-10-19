/**
 * @file    Header.h
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

#ifndef HEADER_H
#define HEADER_H


#include <QCoreApplication>
#include <stdio.h>
#include <stdlib.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <boost/atomic.hpp>
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
#define MAXTHREADS 40
#define MAXREADSIZE 5000
#define MAX_SHIFT_ITER 5
#define WRITEWHEN 100
#define MINACCEVENTDEPTH 0
#define MAINADDRESS "/home/ameer/ExactSV/"
//#define BUFFERSIZE 1024*1024    // Size of fread Buffer
//#define MAX_BUFFER_SIZE 100     // each thread will read MAX_BUFFER_SIZE reads in each trial
// ====

int numberOfThreads = 8;
// ====
int chunkSize = 25;
// ==== Files:
string readsFileName = "readsv2"; //string readsFileName = "reads_t1";    //string readName = "reads_k12.fq";
int readsFileCount = 2;
string genomeName = "Chr19.fa"; //string genomeName = "E_coli_K12_DH10B.fa";    //string genomeName = "chr19.fa";
//string indexName = "Chr19"; //string indexName = "E_coli_K12_DH10B";        //string indexName = "chr19";
string outputDir = "/home/ameer/ExactSV/SV_out/";

string rightEventsFileName = outputDir+"rightEvents";
string leftEventsFileName = outputDir+"leftEvents";
string informativeFileName = outputDir+"informatives";
string informativeFullyFileName = outputDir+"informativesFully";
string doubInformativeFileName = outputDir+"doubInformatives";
string unMappedFileName = outputDir+"unMapped";
//string outputDir = "/home/ubuntu/SVAnchoring/SV_out2/";

// ==== Alignment Parameters:
int dashV = 1;
double gapRate = 5; // rate (#) of acceptable gap caused by mutation or sequencing error per chunk (for correspondance checking)
double indelShift = 0.4;
int editDistance = 2;
int gapPen = -8, misPen = -5, matchPen = 10;
//int numGap = 5;     //numGap in Anchoring | error in Assignment
int LocAlth = 130;
int localAlThreshold = LocAlth;
//int anchoringShift = 5;
//int assignmentShift = 10; //d = chunkSize - assignmentShift
bool runAnchoring = true;
bool inClustersExtension = false;
long accMinDist = 200;
bool BGChanges = false;
bool BGChangesOnAnchor = false;
bool BGLR = true;               // [true]: Left and Right Binary Genomes
                                // [false]: only 1 BinaryGenome (less accurate and slower But less Ram usage)
int shiftIterations = 1;        // max = MAX_SHIFT_ITER
bool logTxt = false;
bool pairedEndedReads = false;  // if true, then reads must be in two files "readName"_1.fq and "readName"_2.fq
bool runPhase1 = false;
bool runPhase2 = true;
bool writeInfo = true;

long long cntCOM = 0, cntBG=0;


// SureMap running options:

/* runing options */
//int mc = 10;
//int core = 1;
//int globalMaxReport = 1;
//int maxReport[MAXTHREADS];
//int globalBestOne = 0;
//string globalMode = "normal";
//int bestOne[MAXTHREADS];
//double globalNoisePercent = -0.1;
//double noisePercent[MAXTHREADS];
//int globalMaxDiffMismatch = 1;
//int maxDiffMismatch[MAXTHREADS];
//int globalMaxDiffEdit = 1;
//int maxDiffEdit[MAXTHREADS];
//int globalGap = 0;
//int gap[MAXTHREADS];
//int globalUniqeOption = 1;
//int globalLongReadNoice = 30;
//int uniqeOption[MAXTHREADS];
//string outputAdr = "report.sam";
//bool longRead = false;
//int mxFailed[MAXTHREADS];
//int minVal[MAXTHREADS];
//int hardToMap = 0;
//int Mode = 20;
//int fragLen = 500;


//                                 |
// ================================|
//> End of Settings:               |
// ================================|
// ================================'


#endif // HEADER_H
