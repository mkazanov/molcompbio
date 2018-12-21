//
//  mutation.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef mutation_hpp
#define mutation_hpp

#define STRLEN_CANCER 4
#define STRLEN_SAMPLE 16
#define STRLEN_CHR 2
#define STRLEN_REFALLELE 1
#define STRLEN_VARALLELE 1

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "mutsignature.hpp"
#include "ghuman.hpp"

using namespace std;

class CMutation{
public:
    char cancer[STRLEN_CANCER+1];
    char sample[STRLEN_SAMPLE+1];
    char chr[STRLEN_CHR+1];
    unsigned long pos;
    string refallele;
    string varallele;
    char isForwardStrand;
    CMutation(string cancer_,
              string sample_,
              string chr_,
              string pos_,
              string refallele,
              string varallele);
    CMutation(){};
};

class CMutations{
public:
    vector<CMutation> mutations;
    unsigned long mutationsCnt;
    void LoadMutations(string path, int isHeader);
    void FilterMutations(CMutations& filteredMutations,
                         vector<CMutationSignature>& signatures,
                         CHumanGenome& human,
                         set<string> cancers,
                         set<string> samples,
                         CMutations* pOtherMutations);
    void SaveToFile(string path);
};

#endif /* mutation_hpp */
