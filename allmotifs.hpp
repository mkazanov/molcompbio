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
    void GenerateMotifsFile(string path);
    void RunAnalysis(string motif, string dirpath);
};

#endif /* allmotifs_hpp */
