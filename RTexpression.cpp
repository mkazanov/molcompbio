//
//  RTexpression.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 11/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#include "RTexpression.hpp"
#include "mutsignature.hpp"
#include "replicationtime.hpp"
#include <cstring>

map<pair<int,int>,unsigned long> CRTexpression::CalculateMotifsRTexpressionBins(CReplicationTiming& rt, CExpression exp,  vector<CExpressionBin> expBins, CHumanGenes& genes, set<string> motifs, string sample, CHumanGenome* phuman)
{
    int i,j;
    set<string> motifsall;
    set<string>::iterator si;
    int motiflen;
    map<pair<int,int>,unsigned long> ret;
    
    CMutationSignature msobj;
    char** motifsarr;
    int* strandarr;
    msobj.CheckMotifsNotEmpty(motifs);
    msobj.CheckMotifsSameLength(motifs);
    motifsall = msobj.AddcMotifs(motifs);
    motiflen = (int) (motifsall.begin())->length();
    
    // Prepare char array for copying motifs
    motifsarr = new char*[motifsall.size()];
    strandarr = new int[motifsall.size()];
    for(i=0;i<motifsall.size();i++)
    {
        motifsarr[i] = new char[motiflen];
        if(i<motifs.size())
            strandarr[i] = 1;
        else
            strandarr[i] = 0;
    }
    
    // Copy motifs to char array
    i = 0;
    for(si=motifsall.begin();si!=motifsall.end();si++)
    {
        strncpy(motifsarr[i],(*si).c_str(),(*si).length());
        i++;
    }
    
    int motifsnum;
    motifsnum = (int) motifsall.size();
    unsigned long** cnt;
    int rtBinsSize, expBinsSize;
    rtBinsSize = (int)rt.bins.size() + (int)RT_NULLBIN_CNT;
    expBinsSize = (int)expBins.size() + (int)EXP_NULLBIN_CNT;
    cnt = new unsigned long*[rtBinsSize];
    for(i=0;i<rtBinsSize;i++)
        cnt[i] = new unsigned long[expBinsSize];
    
    for(i=0;i<rtBinsSize;i++)
        for(j=0;j<expBinsSize;j++)
            cnt[i][j] = 0;
    
    int includeCurPos=1;
    CDNAPos pos = CDNAPos(0,0);
    int rtbin, expbin;
    int strand;
    int strandInconsistence;
    for(pos=msobj.NextMotif(CDNAPos(0,0),motifsarr,strandarr,motifsnum,motiflen,phuman,END_GENOME,includeCurPos);
        !pos.isNull();
        pos=msobj.NextMotif(pos,motifsarr,strandarr,motifsnum,motiflen,phuman,END_GENOME))
    {
        rtbin = rt.GetRTBin(pos.chrNum, pos.pos, rt.bins);
        expbin = exp.GetExpressionBin(sample, phuman->chrName[pos.chrNum], pos.pos, 1, genes, expBins, strand, strandInconsistence);
        cnt[rtbin+RT_NULLBIN_CNT][expbin+EXP_NULLBIN_CNT]++;
    }
    
    for(i=0;i<rtBinsSize;i++)
        for(j=0;j<expBinsSize;j++)
            ret.insert(make_pair(make_pair(i-RT_NULLBIN_CNT,j-EXP_NULLBIN_CNT),cnt[i][j]));
    
    return(ret);
}
