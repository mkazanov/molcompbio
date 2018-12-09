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

CReplicationTime CReplicationTiming::GetRT(int chrNum, unsigned long pos)
{
    set<CReplicationTime>::iterator it;
    CReplicationTime RTnull(CHR_NULL,0);

    CReplicationTime rt(chrNum, pos);
    it = RTs.upper_bound(rt);
    it--;
    if(it->chrNum == chrNum && pos >= it->startpos && pos <= it->endpos)
        return (*it);
    else
        return RTnull;
}

int CReplicationTiming::GetRTBin(double RTvalue, vector<CRTBin> bins)
{
    for(int i=0;i<bins.size();i++)
        if(RTvalue >= bins[i].RTleft && RTvalue < bins[i].RTright)
            return(bins[i].binNum);
    return(-1);
}

int CReplicationTiming::GetRTBin(CReplicationTime rt, vector<CRTBin> bins)
{
    int bin;
    
    if(rt.isRTnull())
        bin = -1;
    else
        bin = GetRTBin(rt.RTvalue, bins);

    return(bin);
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
                it->isForward = 1;
            else if(nextit->RTvalue - previt->RTvalue < 0)
                it->isForward = 0;
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

