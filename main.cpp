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

int main(int argc, const char * argv[]) {
    
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    CAPOBEC a;
    
    a.ClassifyMutations(&human);
    
    /*a.CalculateAPOBECEnrichment(&human);
    
    // Replication timing
    
    a.AnalysisReplicationTiming(a.apobecMuts,"results_RT_APOBEC.txt");
    a.AnalysisReplicationTiming(a.otherMuts,"results_RT_OTHER.txt");

    a.CalculateTargetsinRTBins(&human);*/
    
    
    // Expression
    
    //a.AnalysisExpression(a.apobecMuts, "results_exp_APOBEC.txt");
    //a.AnalysisExpression(a.otherMuts, "results_exp_OTHER.txt");

    if(argc > 1)
        a.CalculateTargetsinExpressionBins(&human,argv[1],argv[2]);
    else
        a.CalculateTargetsinExpressionBins(&human);
    
    return 0;
}
