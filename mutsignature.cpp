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

CDNAPos CMutationSignature::NextMotif(CDNAPos pos, set<string>motifsall, CHumanGenome* phuman, int end, int includeCurrentPos)
{
    unsigned long motiflen;
    CDNAPos ret = CDNAPos(-1,0);
    int endChrNum;
    set<string>::iterator m;
    int break2;
    unsigned long startpos;
    
    motiflen = (motifsall.begin())->length();
    
    if(end == END_CHROMOSOME)
        endChrNum = pos.chrNum + 1;
    else
        endChrNum = phuman->chrCnt;
    
    if(includeCurrentPos)
        startpos = pos.pos;
    else
        startpos = pos.pos + 1;
    for(int i=pos.chrNum;i<endChrNum;i++)
    {
        for(unsigned long j=startpos;j<(phuman->chrLen[i]-motiflen+1);j++)
        {
            for(m=motifsall.begin();m!=motifsall.end();m++)
            {
                break2 = 0;
                for(int n=0;n<motiflen;n++)
                {
                    if((*m)[n] == 'X')
                        continue;
                    if(phuman->dna[i][j+n] != (*m)[n])
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
    int includeCurPos=1;
    
    CheckMotifsNotEmpty(motifs);
    CheckMotifsSameLength(motifs);
    motifsall = AddcMotifs(motifs);
    
    CDNAPos pos;
    for(pos=NextMotif(CDNAPos(0,0),motifsall,phuman,END_GENOME,includeCurPos);
        !pos.isNull();
        pos=NextMotif(pos,motifsall,phuman,END_GENOME))
        res++;
    
    return(res);
}
