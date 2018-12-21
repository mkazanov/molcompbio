//
//  main.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 07/03/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <iostream>
#include "biogenome.hpp"
#include "ghuman.hpp"
#include "dnatest.hpp"
#include "apobec.hpp"
#include "mutation.hpp"
#include "dna.hpp"

int main(int argc, const char * argv[]) {
    
    CAPOBEC a;
    a.ClassifyMutations();
    //a.AnalysisReplicationTiming(a.apobecMuts,"results_RT_APOBEC.txt");
    a.AnalysisReplicationTiming(a.otherMuts,"results_RT_OTHER.txt");

    //a.CalculateTargetsinRTBins();
    
    return 0;
}
