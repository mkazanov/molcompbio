//
//  genehuman.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 09/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef genehuman_hpp
#define genehuman_hpp

#include <string>
#include <set>

using namespace std;

class CHumanGene {
public:
    int chrNum;
    unsigned long startpos;
    unsigned long endpos;
    mutable unsigned long maxendpos;
    char strand;
    unsigned long geneID;
    string geneName;
    CHumanGene(string chr_, string starpos_, string endpos_, string strand_, string info_);
    CHumanGene(int chrNum_, unsigned long pos_);
    bool operator< (const CHumanGene &right) const
    {
        if (chrNum < right.chrNum)
            return true;
        else if (chrNum == right.chrNum)
            return startpos < right.startpos;
        else
            return false;
    }
};

class CHumanGenes {
public:
    set<CHumanGene> genes;
    void LoadGenes(string path);
    void PrepareForSearch();
    int GetGenesByPos(int chrNum, unsigned long pos, vector<CHumanGene>& geneList);
    void SaveToFile(string path);
};

#endif /* genehuman_hpp */
