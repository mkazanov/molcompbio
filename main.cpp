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
#include "dnatest.hpp"
#include "apobec.hpp"
#include "mutation.hpp"
#include "dna.hpp"

int main(int argc, const char * argv[]) {
    
        CMutationSignature s;
        set<string> motifs;
    int cytosineCnt;
    

    
    CAPOBEC a;

    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
 /*   a.ClassifyMutations(&human);

    a.CalculateAPOBECEnrichment(&human);
    
    a.AnalysisReplicationTiming(a.apobecMuts,"results_RT_APOBEC.txt");
    a.AnalysisReplicationTiming(a.otherMuts,"results_RT_OTHER.txt");*/

    a.CalculateTargetsinRTBins(&human);
    
    //a.AnalysisExpression(a.apobecMuts, "results_expression_APOBEC.txt");
    //a.AnalysisExpression(a.otherMuts, "results_expression_OTHER.txt");
    
    return 0;
}
