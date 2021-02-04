//
//  apobec.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 14/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#include "apobec.hpp"
#include "options.h"
#include <iostream>
#include "biogenome.hpp"
#include "ghuman.hpp"
#include "signanalysis.hpp"
#include "mutation.hpp"
#include "dna.hpp"
#include <fstream>
#include <ctime>
#include "signanalysis.hpp"

void CAPOBEC::RunAnalysis(int argc, const char * argv[])
{


time_t t;
tm* now;

// Start date/time
t = time(0);   // get time now
now = localtime(&t);
cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

//

CHumanGenome human;
human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");

CSignatureAnalysis a;
    
a.signatures.push_back(CMutationSignature("TCA",2,"T"));
a.signatures.push_back(CMutationSignature("TCT",2,"T"));
a.signatures.push_back(CMutationSignature("TCA",2,"G"));
a.signatures.push_back(CMutationSignature("TCT",2,"G"));
a.signatures.push_back(CMutationSignature("TCC",2,"T"));
a.signatures.push_back(CMutationSignature("TCG",2,"T"));
a.signatures.push_back(CMutationSignature("TCC",2,"G"));
a.signatures.push_back(CMutationSignature("TCG",2,"G"));

vector<string> cancers;
vector<string> samples;
CMutations m;
m.LoadMutations(CANCER_MUTATIONS_FILEFORMAT, CANCER_MUTATIONS_PATH, cancers, samples, 1, &human);

a.ClassifyMutations(m,cancers,samples,&human);

/*a.CalculateAPOBECEnrichment(&human);*/
 
 // Replication timing
 
 //a.AnalyzeReplicationTiming(a.signatureMuts,"results_RT_APOBEC.txt");
 //a.AnalyzeReplicationTiming(a.otherMuts,"results_RT_OTHER.txt");
 a.CalculateTargetsinRTBins("TCX",&human,1);
 //a.CalculateTargetsinRTBins("ALL",&human,0);


// Expression

//a.AnalyzeExpression(a.apobecMuts, "results_exp_APOBEC.txt");
//a.AnalyzeExpression(a.otherMuts, "results_exp_OTHER.txt");

//if(argc > 1)
//    a.CalculateTargetsinExpressionBins("ALL_in_EXPbins",&human,argv[1],argv[2],0);
//else
//   a.CalculateTargetsinExpressionBins("ALL_in_EXPbins",&human);


// RT + expression

//a.AnalyzeRTExpression(a.apobecMuts,"results_RTexp_APOBEC.txt");
//a.AnalyzeRTExpression(a.otherMuts,"results_RTexp_OTHER.txt");

//if(argc > 1)
// a.CalculateTargetsinRTexpressionBins("TCW_in_RTEXPbins",&human,argv[1],argv[2],1);
//else
// a.CalculateTargetsinRTexpressionBins("TCW_in_RTEXPbins",&human);


// End date/time
t = time(0);
now = localtime(&t);
cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

}
