//
//  expression.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 16/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef expression_hpp
#define expression_hpp

#include <string>
#include <map>
#include <vector>
#include <set>
#include "ghuman.hpp"
#include "genehuman.hpp"

#define EXP_NULLBIN_NOTINGENES  -2
#define EXP_NULLBIN_NOEXPDATA   -1
#define EXP_NULLBIN_CNT 2

using namespace std;

class CExpressionKey
{
public:
    unsigned long geneId;
    string sample;
    CExpressionKey(string geneId_, string sample_);
    CExpressionKey(unsigned long geneId_, string sample_);
    bool operator< (const CExpressionKey &right) const
    {
        if (geneId < right.geneId)
            return true;
        else if (geneId > right.geneId)
            return false;
        else
        {
            if (sample < right.sample)
                return true;
            else if (sample >= right.sample)
                return false;
        }
        return false;
    }
};

class CExpressionBin {
public:
    int binNum;
    double expressionLeft;
    double expressionRight;
    CExpressionBin(int binNum_, double expressionLeft_, double expressionRight_);
};

class CExpression
{
public:
    map<CExpressionKey, double> data;
    void LoadExpression(string path);
    int GetExpression(unsigned long geneId, string sample, double& expressionValue);
    int GetExpressionBinByValue(double expValue, vector<CExpressionBin> bins);
    int GetExpressionBin(string sample, string chr, unsigned long pos, char isForwardMut, CHumanGenes& genes, vector<CExpressionBin>& expBins, int& strand, int& strandInconsistence);
    map<int,unsigned long> CalculateMotifsExpressionBins(vector<CExpressionBin> expBins, CHumanGenes& genes, set<string> motifs, string sample, CHumanGenome* phuman);
};
#endif /* expression_hpp */
