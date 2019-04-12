//
//  main.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 07/03/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "options.h"
#include <iostream>
#include "biogenome.hpp"
#include "ghuman.hpp"
#include "apobec.hpp"
#include "mutation.hpp"
#include "dna.hpp"
#include <fstream>
#include <ctime>

int main(int argc, const char * argv[]) {
    
    time_t t;
    tm* now;
    
    // Start date/time
    t = time(0);   // get time now
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";
    
    //
    
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    CAPOBEC a;
    
    a.ClassifyMutations(&human);
    
    /*a.CalculateAPOBECEnrichment(&human);
    
    // Replication timing
    
    a.AnalyzeReplicationTiming(a.apobecMuts,"results_RT_APOBEC.txt");
    a.AnalyzeReplicationTiming(a.otherMuts,"results_RT_OTHER.txt");
    a.CalculateTargetsinRTBins(&human);*/
    
    
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

    if(argc > 1)
        a.CalculateTargetsinRTexpressionBins("TCW_in_RTEXPbins",&human,argv[1],argv[2],1);
    else
       a.CalculateTargetsinRTexpressionBins("TCW_in_RTEXPbins",&human);

    
    // End date/time
    t = time(0);
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

    
    return 0;
}
