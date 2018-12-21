//
//  apobec.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef apobec_hpp
#define apobec_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include "mutation.hpp"
#include "genehuman.hpp"
#include "expression.hpp"
#include "replicationtime.hpp"

using namespace std;

class CResultsKey {
public:
    string cancer;
    string sample;
    int bin;
    CResultsKey(string cancer_, string sample_, int bin_);
    bool operator< (const CResultsKey &right) const
    {
        if (cancer < right.cancer)
            return true;
        else if (cancer > right.cancer)
            return false;
        else
        {
            if (sample < right.sample)
                return true;
            else if (sample > right.sample)
                return false;
            else
            {
                if(bin < right.bin)
                    return true;
                else
                    return false;
            }
        }
    }

};

class CResultsValue{
public:
    unsigned long mutCnt;
    unsigned long leadingCnt;
    unsigned long laggingCnt;
    unsigned long plusStrandConsistent;
    unsigned long minusStrandConsistent;
    unsigned long plusStrandAll;
    unsigned long minusStrandAll;
    
    CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_);
    CResultsValue(unsigned long mutCnt_, unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll, unsigned long minusStrandAll);
    CResultsValue(){};
};

class CAPOBEC {
public:
    CMutations apobecMuts;
    CMutations otherMuts;
    void ClassifyMutations();
    void AnalysisReplicationTiming(CMutations& muts, string resultsFilename);
    void AnalysisExpression();
    int GetExpressionBin(string sample, string chr, unsigned long pos, char isForwardMut, CHumanGenes& genes, CExpression& exp, vector<CExpressionBin>& expBins, int& strand, int& strandInconsistence);
    void CalculateTargetsinRTBins();
};
    
#endif /* apobec_hpp */
