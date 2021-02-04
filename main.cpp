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
#include "CASEclusters.hpp"

int main(int argc, const char * argv[]) {
    
    //ProcessClusters();
    MakeApobecMutationsFile();
    return 0;

    
    
    vector<string> cancers;
    vector<string> samples;
    
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    CMutations m;
    cancers.push_back("CESC");
    m.LoadMutations(CANCER_MUTATIONS_FILEFORMAT, CANCER_MUTATIONS_PATH, cancers, samples, 1, &human);
    m.RenameSamples("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_samples_list.csv",2,3,16);

    CSignatureAnalysis a;
    
    a.signatures.push_back(CMutationSignature("TCA",2,"T"));
    a.signatures.push_back(CMutationSignature("TCT",2,"T"));
    a.signatures.push_back(CMutationSignature("TCA",2,"G"));
    a.signatures.push_back(CMutationSignature("TCT",2,"G"));
    a.signatures.push_back(CMutationSignature("TCC",2,"T"));
    a.signatures.push_back(CMutationSignature("TCG",2,"T"));
    a.signatures.push_back(CMutationSignature("TCC",2,"G"));
    a.signatures.push_back(CMutationSignature("TCG",2,"G"));
    
    a.ClassifyMutations(m,cancers,samples,&human);

    a.CalculateAPOBECEnrichment(&human);
    
    return 0;
    
    CAllMotifs c;
    c.AnalysisRTExp(string(RESULTS_FOLDER)+"/ALL",argv[1],argv[2]);
    
    return 0;
    
    
    
    time_t t;
    tm* now;
    
    // Start date/time
    t = time(0);   // get time now
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";
    

    //CHumanGenome human;
    //human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
   /* CSignatureAnalysis a;
    
    a.signatures.push_back(CMutationSignature("TCA",2,"T"));
    a.signatures.push_back(CMutationSignature("TCT",2,"T"));
    a.signatures.push_back(CMutationSignature("TCA",2,"G"));
    a.signatures.push_back(CMutationSignature("TCT",2,"G"));
    a.signatures.push_back(CMutationSignature("TCC",2,"T"));
    a.signatures.push_back(CMutationSignature("TCG",2,"T"));
    a.signatures.push_back(CMutationSignature("TCC",2,"G"));
    a.signatures.push_back(CMutationSignature("TCG",2,"G"));
    
    a.ClassifyMutations(&human);

    return 0; */
    
/*
    CMutations m;
    int isHeader = 1;
    m.LoadMutations(CANCER_MUTATIONS, isHeader);
    
    CSignatureAnalysis a;
    a.CalculateExpressionAllMotifs(m, "MOTIFS_EXP/ALLMUTS_in_EXPbins", &human, argv[1], argv[2]);
    a.CalculateTargetsinExpBinAllMotifs("MOTIFS_EXP/ALLMOTIFS_in_EXPbins", &human, argv[1], argv[2]);
    */
 
    
    
    // End date/time
    t = time(0);
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

    return 0;
    
    
    //CAPOBEC a;
    //a.RunAnalysis(argc, argv);
    
    return 0;
}
