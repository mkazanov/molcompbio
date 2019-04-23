//
//  apobec.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef signanalysis_hpp
#define signanalysis_hpp

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
    int rtbin;
    int expbin;
    CResultsKey(string cancer_, string sample_, int rtbin_, int expbin_);
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
                if(rtbin < right.rtbin)
                    return true;
                else if (rtbin > right.rtbin)
                    return false;
                else
                {
                    if (expbin < right.expbin)
                        return true;
                    else if (expbin > right.expbin)
                        return false;
                    else
                        return false;
                }
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
    CResultsValue(unsigned long mutCnt_, unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll_, unsigned long minusStrandAll_);
    CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_,unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll_, unsigned long minusStrandAll_);
    CResultsValue(){};
};

class CMapKey {
public:
    string cancer;
    string sample;
    CMapKey(string cancer_, string sample_)
    {
        cancer = cancer_;
        sample = sample_;
    }
    bool operator< (const CMapKey &right) const
    {
        if (cancer < right.cancer)
            return true;
        else if (cancer > right.cancer)
            return false;
        else
        {
            if (sample < right.sample)
                return true;
            else
                return false;
        }
    }
};

class CMapValue {
public:
    unsigned long APOBECmutsCnt;
    unsigned long cytosineMutsCnt;
    double enrichment;
    double enrichmentExcludeTCW;
    CMapValue(unsigned long APOBECmutsCnt_, unsigned long cytosineMutsCnt_, double enrichment_, double enrichmentExcludeTCW_)
    {
        APOBECmutsCnt = APOBECmutsCnt_;
        cytosineMutsCnt = cytosineMutsCnt_;
        enrichment = enrichment_;
        enrichmentExcludeTCW = enrichmentExcludeTCW_;
    }
};

class CSignatureAnalysis {
public:
    vector<CMutationSignature> signatures;
    CMutations signatureMuts;
    CMutations otherMuts;
    void ClassifyMutations(CHumanGenome* phuman_ = NULL);
    void AnalyzeReplicationTiming(CMutations& muts, string resultsFilename);
    void AnalyzeExpression(CMutations& muts, string resultsFilename);
    void AnalyzeRTExpression(CMutations& muts, string resultsFilename);
    int GetExpressionBin(string sample, string chr, unsigned long pos, char isForwardMut, CHumanGenes& genes, CExpression& exp, vector<CExpressionBin>& expBins, int& strand, int& strandInconsistence);
    void CalculateTargetsinRTBins(string motifname, CHumanGenome* phuman_ = NULL, int isSignatureMotif = 1);
    void CalculateTargetsinExpressionBins(string outFilePrefix, CHumanGenome* phuman = NULL, string cancer = "", string sample = "", int isAPOBECmotif = 1);
    void CalculateTargetsinRTexpressionBins(string outFilePrefix, CHumanGenome* phuman = NULL, string cancer = "", string sample = "", int isAPOBECmotif = 1);
    
    void CalculateAPOBECEnrichment(CHumanGenome* phuman_ = NULL);
    set<CCancerSample> LoadCancerSamples(string path);
};
    
#endif /* signanalysis_hpp */
