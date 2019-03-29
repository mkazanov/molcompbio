//
//  mutsignature.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "mutsignature.hpp"
#include <set>
#include "options.h"
#include <iostream>

CMutationSignature::CMutationSignature(string motif_, int mutationPos_, string newbase_)
{
    motif = motif_;
    mutationPos = mutationPos_;
    newbase = newbase_;
}

bool CMutationSignature::AnyNewBase()
{
    return((newbase == "X" || newbase == "N"));
}

int CMutationSignature::CheckMotifsNotEmpty(set<string> motifs)
{
    if(motifs.empty())
    {
        cerr << "Error: motifs array is empty." << '\n';
        return(0);
    }
    return(1);
}

int CMutationSignature::CheckMotifsSameLength(set<string> motifs)
{
    set<string>::iterator m;
    unsigned long motiflen;
    
    motiflen = (motifs.begin())->length();
    for(m=motifs.begin();m!=motifs.end();m++)
        if(motiflen != (*m).length())
        {
            cerr << "Error: motifs have different length" << '\n';
            return(0);
        }
    return(1);
}

set<string>  CMutationSignature::AddcMotifs(set<string> motifs)
{
    set<string> motifsall;
    set<string>::iterator m;
    
    motifsall = motifs;
    for(m=motifs.begin();m!=motifs.end();m++)
    {
        motif = CDNA::cDNA(*m);
        motifsall.insert(motif);
    }
    
    return(motifsall);
}

CDNAPos CMutationSignature::NextMotif(CDNAPos pos, char** motifsarr, int motifsnum, int motiflen, CHumanGenome* phuman, int end, int includeCurrentPos)
{
    CDNAPos ret = CDNAPos(-1,0);
    int endChrNum;
    int i,k,n;
    unsigned long j;
    int break2;
    unsigned long startpos;
    
    if(end == END_CHROMOSOME)
        endChrNum = pos.chrNum + 1;
    else
        endChrNum = phuman->chrCnt;
    
    if(includeCurrentPos)
        startpos = pos.pos;
    else
        startpos = pos.pos + 1;
    
    for(i=pos.chrNum;i<endChrNum;i++)
    {
        for(j=startpos;j<(phuman->chrLen[i]-motiflen+1);j++)
        {
            for(k=0;k<motifsnum;k++)
            {
                break2 = 0;
                for(n=0;n<motiflen;n++)
                {
                    if(motifsarr[k][n] == 'X')
                        continue;
                    if(phuman->dna[i][j+n] != motifsarr[k][n])
                    {
                        break2 = 1;
                        break;
                    }
                }
                if(break2)
                    continue;
                ret.chrNum = i;
                ret.pos = j;
                return(ret);
            }
        }
        startpos = 0;
    }
    return(ret);
}

unsigned long CMutationSignature::CountMotifGenome(set<string> motifs, CHumanGenome* phuman)
{
    unsigned long res = 0;
    set<string> motifsall;
    set<string>::iterator si;
    int includeCurPos=1;
    char** motifsarr;
    int motiflen, motifsnum;
    int i;
    
    CheckMotifsNotEmpty(motifs);
    CheckMotifsSameLength(motifs);
    motifsall = AddcMotifs(motifs);
    motiflen = (motifsall.begin())->length();
    motifsnum = motifsall.size();
    
    // Prepare char array for copying motifs
    motifsarr = new char*[motifsall.size()];
    for(i=0;i<motifsall.size();i++)
        motifsarr[i] = new char[motiflen];
    
    // Copy motifs to char array
    i = 0;
    for(si=motifsall.begin();si!=motifsall.end();si++)
    {
        strncpy(motifsarr[i],(*si).c_str(),(*si).length());
        i++;
    }
    
    CDNAPos pos;
    for(pos=NextMotif(CDNAPos(0,0),motifsarr,motifsnum,motiflen,phuman,END_GENOME,includeCurPos);
        !pos.isNull();
        pos=NextMotif(pos,motifsarr,motifsnum,motiflen,phuman,END_GENOME))
        res++;
    
    return(res);
}
