//
//  replicationtime.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "replicationtime.hpp"
#include "service.h"
#include <fstream>
#include "ghuman.hpp"
#include "dna.hpp"
#include <map>
#include "options.h"
#include <iostream>
#include <algorithm>
#include "mutsignature.hpp"

CReplicationTime::CReplicationTime(string chr_, string startpos_, string endpos_, string RTvalue_)
{
    chrNum = CHumanGenome::GetChrNum(chr_.substr(3));
    startpos = str2ul(startpos_);
    endpos = str2ul(endpos_);
    RTvalue = str2d(RTvalue_);
}

CReplicationTime::CReplicationTime(int chrNum_, unsigned long pos_)
{
    chrNum = chrNum_;
    startpos = pos_;
}

bool CReplicationTime::isRTnull()
{
    if(chrNum == CHR_NULL)
        return true;
    else
        return false;
}

CRTBin::CRTBin(int binNum_, double RTleft_, double RTright_)
{
    binNum = binNum_;
    RTleft = RTleft_;
    RTright = RTright_;
}

void CReplicationTiming::LoadReplicationTiming(string path, int isHeader)
{
    string line;
    clock_t c1,c2;
    
    ifstream f(path);
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }
    
    if(isHeader)
        getline(f, line);
    
    printf("Replication timing loading ...\n");
    int i=0;
    vector<string> flds;
    c1 = clock();
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            RTs.emplace(flds[0],flds[1],flds[2],flds[3]);
            i++;
        }
    }
    c2 = clock();
    printf("Replication timing %i intervals have been loaded\n", RTs.size());
    printf("Executing time: %lu \n", c2 - c1);
    int size = RTs.size();
}

int CReplicationTiming::GetRT(int chrNum, unsigned long pos, double& RTvalue)
{
    set<CReplicationTime>::iterator it;

    CReplicationTime rt(chrNum, pos);
    it = RTs.upper_bound(rt);
    if(it == RTs.begin())
    {
        RTvalue = RT_NULL;
        return(0);
    }
    it--;
    if(it != RTs.end() && it->chrNum == chrNum && pos >= it->startpos && pos <= it->endpos)
    {
        RTvalue = it->RTvalue;
        return (1);
    }
    else
    {
        RTvalue = RT_NULL;
        return(0);
    }
}

int CReplicationTiming::GetRT(int chrNum, unsigned long pos, CReplicationTime& rt)
{
    set<CReplicationTime>::iterator it;
    
    CReplicationTime rtquery(chrNum, pos);
    it = RTs.upper_bound(rtquery);
    if(it == RTs.begin())
        return(0);
    it--;
    if(it != RTs.end() && it->chrNum == chrNum && pos >= it->startpos && pos <= it->endpos)
    {
        rt = (*it);
        return (1);
    }
    else
        return(0);
}

int CReplicationTiming::GetRTBin(double RTvalue, vector<CRTBin> bins)
{
    for(int i=0;i<bins.size();i++)
        if(RTvalue >= bins[i].RTleft && RTvalue < bins[i].RTright)
            return(bins[i].binNum);
    return(RT_NULLBIN_NOBIN);
}

int CReplicationTiming::GetRTBin(int chrNum, unsigned long pos, vector<CRTBin> bins)
{
    double RTvalue;
    int res;
    int bin;
    
    res = GetRT(chrNum,pos,RTvalue);
    if(res)
    {
        bin = GetRTBin(RTvalue,bins);
        return(bin);
    }
    else
        return(RT_NULLBIN_NOVALUE);
}

void CReplicationTiming::ReplicationStrand()
{
    set<CReplicationTime>::iterator it,previt,nextit;
    
    for(it=RTs.begin();it!=RTs.end();++it)
    {
        if(it != RTs.begin())
            previt = prev(it);
        nextit = next(it);
        if(it == RTs.begin() || previt->chrNum != it->chrNum || nextit->chrNum != it->chrNum)
            it->isForward = -1;
        else
        {
            if(nextit->RTvalue - previt->RTvalue > 0)
                it->isForward = 0;
            else if(nextit->RTvalue - previt->RTvalue < 0)
                it->isForward = 1;
            else
                it->isForward = -1;
        }
    }
}

void CReplicationTiming::SaveToFile(string path)
{
    ofstream f;
    f.open(path);
    set<CReplicationTime>::iterator it;

    for(it=RTs.begin();it!=RTs.end();++it)
        f << (*it).chrNum << '\t' << (*it).startpos << '\t' << (*it).endpos << '\t' << (*it).RTvalue << '\t' << (*it).isForward << '\n';
        
    f.close();
}

int CReplicationTiming::CalculateMotifinRTBins(vector<CRTBin> bins, set<string> motifs, string OUT_PATH, CHumanGenome* phuman)
{
    CMutationSignature msobj;
    int bin;
    map<int,unsigned long> results;
    map<int,unsigned long>::iterator it;
    set<string> motifsall;
    long motiflen;
    
    msobj.CheckMotifsNotEmpty(motifs);
    msobj.CheckMotifsSameLength(motifs);
    motifsall = msobj.AddcMotifs(motifs);
    motiflen = (motifsall.begin())->length();

    results[RT_NULLBIN_NOVALUE] = 0;
    results[RT_NULLBIN_NOBIN] = 0;
    for(int i=0;i<bins.size();i++)
        results[bins[i].binNum] = 0;
    
    CDNAPos pos = CDNAPos(0,0);
    int includeCurPos=1;
    int chrNum = 0;
    cout << "Chr:" << pos.chrNum << ", Pos:" << pos.pos << '\n';
    for(pos=msobj.NextMotif(CDNAPos(0,0),motifsall,phuman,END_GENOME,includeCurPos);
        !pos.isNull();
        pos=msobj.NextMotif(pos,motifsall,phuman,END_GENOME))
    {
        bin = GetRTBin(pos.chrNum,pos.pos+(motiflen/2)+1,bins);
        results[bin]++;
        if(chrNum != pos.chrNum)
        {
            cout << "Chr:" << pos.chrNum << ", Pos:" << pos.pos << '\n';
            chrNum = pos.chrNum;
        }
    }
    
    ofstream f;
    f.open(OUT_PATH);
    
    f << "ReplicationBin" << '\t' << "TargetCnt" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first << '\t' << it->second <<'\n';
    
    f.close();
    
    return(1);
}
