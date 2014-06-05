#ifndef USERDEFINEDDATA_H
#define USERDEFINEDDATA_H

const int NGenPtBins=25;
const int NETA = 82;
const double GenPt[NGenPtBins+1] = {12,15,20,27,35,45,57,72,90,120,150,200,300,400,550,750,1000,1500,2000,2500,3000,3500,4000,4500,5000,6000};
const char JetResponseFilename[1024]     = "L2L3Production_iterativeCone5.root"; 
const char FitterOutputFilename[1024]    = "FitterResults_iterativeCone5.root";
const char L3OutputROOTFilename[1024]    = "L3Graphs_iterativeCone5.root";
const char L3ResponseFilename[1024]      = "L3Response_iterativeCone5.txt";
const char L3CorrectionFilename[1024]    = "CMSSW_152_L3Absolute_iterativeCone5.txt";
const char L2OutputROOTFilename[1024]    = "L2Graphs_iterativeCone5.root";
const char L2CorrectionFilename[1024]    = "CMSSW_152_L2Relative_iterativeCone5.txt";
const double eta_boundaries[NETA+1]={-5.191,-4.889,-4.716,-4.538,-4.363,-4.191,-4.013,-3.839,-3.664, 
			             -3.489,-3.314,-3.139,-2.964,-2.853,-2.650,-2.500,-2.322,-2.172,
                                     -2.043,-1.930,-1.830,-1.740,-1.653,-1.566,-1.479,-1.392,-1.305,
			             -1.218,-1.131,-1.044,-0.957,-0.879,-0.783,-0.696,-0.609,-0.522,
			             -0.435,-0.348,-0.261,-0.174,-0.087,0.000,0.087,0.174,0.261,0.348,
                                      0.435,0.522,0.609,0.696,0.783,0.879,0.957,1.044,1.131,1.218,1.305,
                                      1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.500,
                                      2.650,2.853,2.964,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,
                                      4.716,4.889,5.191};

#endif
