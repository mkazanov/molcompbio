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
#include <unordered_map>
#include <vector>
#include <set>
#include "ghuman.hpp"
#include "genehuman.hpp"
#include <functional>

#define EXP_NULLBIN_NOTINGENES  -2
#define EXP_NULLBIN_NOEXPDATA   -1
#define EXP_NULLBIN_CNT 2

#define EXP_STRAND_PLUS 1
#define EXP_STRAND_MINUS    0
#define EXP_STRAND_NULL -1


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
    bool operator==(const CExpressionKey &right) const
    {
        if (geneId == right.geneId && sample == right.sample)
            return true;
        else
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
    CExpression(){}
    struct hash_pair
    {
        size_t operator()(const CExpressionKey& k) const
        {
            size_t hash1 = hash<unsigned long>()(k.geneId);
            size_t hash2 = hash<string>()(k.sample);
            return hash1 ^ hash2;
        }
    };
    
    unordered_map<CExpressionKey, double, hash_pair> data;
    unordered_map<unsigned long, double> dataSample;
    void LoadExpression(string path, string sample = "");
    int GetExpression(unsigned long geneId, string sample, double& expressionValue);
    int GetExpressionBinByValue(double expValue, vector<CExpressionBin> bins);
    int GetExpressionBin(string sample, string chr, unsigned long pos, char isForwardMut, CHumanGenes& genes, vector<CExpressionBin>& expBins, int& strand, int& strandInconsistence);
    map<int,unsigned long> CalculateMotifsExpressionBins(vector<CExpressionBin> expBins, CHumanGenes& genes, set<string> motifs, string sample, CHumanGenome* phuman);
    static int oppositeStrand(int strand)
    {
        if(strand == EXP_STRAND_PLUS)
            return(EXP_STRAND_MINUS);
        else if(strand == EXP_STRAND_MINUS)
            return(EXP_STRAND_PLUS);
        else
            return(EXP_STRAND_NULL);
    }
};
#endif /* expression_hpp */
