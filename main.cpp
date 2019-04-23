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

int main(int argc, const char * argv[]) {
    
    CAllMotifs m;
    m.GenerateMotifsFile(RESULTS_FOLDER+string("/all_motifs.txt"));
    
    m.RunAnalysis("AGC", "MOTIFS");
    exit(0);
    
    CAPOBEC a;
    
    a.RunAnalysis(argc, argv);
    
    return 0;
}
