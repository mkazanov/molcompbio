//
//  main.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 07/03/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "apobec.hpp"
#include "allmotifs.hpp"
#include "options.h"
#include <fstream>
#include <iostream>
#include <set>
#include "mutation.hpp"
#include "signanalysis.hpp"

int main(int argc, const char * argv[]) {
    
    
    time_t t;
    tm* now;
    
    // Start date/time
    t = time(0);   // get time now
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";
    
    /*
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    CMutations m;
    int isHeader = 1;
    m.LoadMutations(CANCER_MUTATIONS, isHeader);
    
    CSignatureAnalysis a;
    a.CalculateExpressionAllMotifs(m, "MOTIFS_EXP/ALLMUTS_in_EXPbins", &human, argv[1], argv[2]);
    a.CalculateTargetsinExpBinAllMotifs("MOTIFS_EXP/ALLMOTIFS_in_EXPbins", &human, argv[1], argv[2]);
    */
 
    
    CAllMotifs c;
    c.AnalysisRTExp("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL",argv[1],argv[2]);
    
    // End date/time
    t = time(0);
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

    return 0;
    
    
    CAPOBEC a;
    a.RunAnalysis(argc, argv);
    
    return 0;
}
