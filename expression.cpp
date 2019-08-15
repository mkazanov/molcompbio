//
//  expression.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 16/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "expression.hpp"
#include <fstream>
#include "service.h"
#include <iostream>
#include "mutsignature.hpp"
#include <cstring>

CExpressionKey::CExpressionKey(string geneId_, string sample_)
{
    geneId = str2ul(geneId_);
    sample = sample_;
}

CExpressionKey::CExpressionKey(unsigned long geneId_, string sample_)
{
    geneId = geneId_;
    sample = sample_;
}

CExpressionBin::CExpressionBin(int binNum_, double expressionLeft_, double expressionRight_)
{
    binNum = binNum_;
    expressionLeft = expressionLeft_;
    expressionRight = expressionRight_;
}


void CExpression::LoadExpression(string path, string sample)
{
    cout << "Load expression" << '\n';
    
    string line;
    string geneId,smpl;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        return;
    }
    
    vector<string> flds,flds1;
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            smpl = flds[1].substr(0,16);
            if(sample != "" && smpl != sample)
                continue;
            flds1 = splitd(flds[0],'|');
            geneId = flds1[1];
            if(sample == "")
                data.insert(pair<CExpressionKey,double>(CExpressionKey(geneId,sample),str2d(flds[2])));
            else
                dataSample[str2ul(geneId)] = str2d(flds[2]);
        }
    }
}

int CExpression::GetExpression(unsigned long geneId, string sample, double& expressionValue)
{
    unordered_map<CExpressionKey,double,hash_pair>::iterator it;
    unordered_map<unsigned long,double>::iterator its;
    
    if(!data.empty())
    {
        it = data.find(CExpressionKey(geneId,sample));
        if(it == data.end())
        {
            expressionValue = -1;
            return(0);
        }
        else
        {
            expressionValue = it->second;
            return(1);
        }
    }
    else if (!dataSample.empty())
    {
        its = dataSample.find(geneId);
        if(its == dataSample.end())
        {
            expressionValue = -1;
            return(0);
        }
        else
        {
            expressionValue = its->second;
            return(1);
        }
    }
    else
    {
        cerr << "Error: expression data is empty!" << '\n';
        return(0);
    }
}

int CExpression::GetExpressionBinByValue(double expValue, vector<CExpressionBin> bins)
{
    for(int i=0;i<bins.size();i++)
        if(expValue >= bins[i].expressionLeft && expValue < bins[i].expressionRight)
            return(bins[i].binNum);
    return(-1);
}

int CExpression::GetExpressionBin(string sample, string chr, unsigned long pos, char isForwardMut, CHumanGenes& genes, vector<CExpressionBin>& expBins, int& strand, int& strandInconsistence)
{
    vector<CHumanGene> geneList;
    double expValue, maxexp;
    int strandplus=0,strandminus=0;
    int res;
    int bin;
    int expFound;
    
    genes.GetGenesByPos(CHumanGenome::GetChrNum(string(chr)), pos, geneList);
    
    if(geneList.empty())
    {
        strand = EXP_STRAND_NULL;
        strandInconsistence = 0;
        return(EXP_NULLBIN_NOTINGENES); // mutation not in genes
    }
    maxexp = -100000.0;
    expFound = 0;
    strand = 0;
    for(int i=0;i<geneList.size();i++)
    {
        if((geneList[i].strand == '+' && isForwardMut == 1) || (geneList[i].strand == '-' && isForwardMut == 0))
            strandplus++;
        else if((geneList[i].strand == '-' && isForwardMut == 1) || (geneList[i].strand == '+' && isForwardMut == 0))
            strandminus++;
        res = GetExpression(geneList[i].geneID, sample, expValue);
        if(res)
        {
            if(maxexp < expValue)
            {
                maxexp = expValue;
                if((geneList[i].strand == '+' && isForwardMut == 1) || (geneList[i].strand == '-' && isForwardMut == 0))
                    strand = EXP_STRAND_PLUS;
                else if((geneList[i].strand == '-' && isForwardMut == 1) || (geneList[i].strand == '+' && isForwardMut == 0))
                    strand = EXP_STRAND_MINUS;
                else
                    strand = EXP_STRAND_NULL;
            }
            expFound = 1;
        }
    }
    
    if(strandplus !=0 && strandminus !=0)
        strandInconsistence = 1;
    else
        strandInconsistence = 0;
    
    if(expFound == 0)
        return(EXP_NULLBIN_NOEXPDATA); // mutation in genes, but no expression data
    
    bin = GetExpressionBinByValue(maxexp, expBins);
    return(bin);
}

map<int,unsigned long> CExpression::CalculateMotifsExpressionBins(vector<CExpressionBin> expBins, CHumanGenes& genes, set<string> motifs, string sample, CHumanGenome* phuman)
{
    int i;
    set<string> motifsall;
    set<string>::iterator si;
    int motiflen;
    map<int,unsigned long> ret;
    
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
    unsigned long* cnt;
    cnt = new unsigned long[(int)expBins.size() + (int)EXP_NULLBIN_CNT];
    
    for(i=0;i<((int)expBins.size() + (int)EXP_NULLBIN_CNT);i++)
        cnt[i] = 0;
    
    int includeCurPos=1;
    CDNAPos pos = CDNAPos(0,0);
    int bin;
    int strand;
    int strandInconsistence;
    for(pos=msobj.NextMotif(CDNAPos(0,0),motifsarr,strandarr,motifsnum,motiflen,phuman,END_GENOME,includeCurPos);
        !pos.isNull();
        pos=msobj.NextMotif(pos,motifsarr,strandarr,motifsnum,motiflen,phuman,END_GENOME))
    {
        
        bin = GetExpressionBin(sample, phuman->chrName[pos.chrNum], pos.pos, 1, genes, expBins, strand, strandInconsistence);        
        cnt[bin+EXP_NULLBIN_CNT]++;
    }

    for(i=0;i<((int)expBins.size() + (int)EXP_NULLBIN_CNT);i++)
        ret[i-EXP_NULLBIN_CNT] = cnt[i];
    
    return(ret);
}

