//
//  allmotifs.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 14/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#ifndef allmotifs_hpp
#define allmotifs_hpp

#include <stdio.h>
#include <string>

using namespace std;

class CAllMotifs{
public:
    void GenerateMotifsFile(string path, int includeComplimentary = 0);
    void AnalysisReplicationTiming(string motif, string dirpath);
    void AnalysisRTExp(string dirpath, string cancer, string sample);
};

#endif /* allmotifs_hpp */
