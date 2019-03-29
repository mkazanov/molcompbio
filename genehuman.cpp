//
//  genehuman.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 09/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "genehuman.hpp"
#include <fstream>
#include "service.h"
#include "ghuman.hpp"
#include <iostream>

CHumanGene::CHumanGene(string chr_, string startpos_, string endpos_, string strand_, string info_)
{
    size_t pos1,pos2;
    
    chrNum = CHumanGenome::GetChrNum(chr_.substr(3));
    startpos = str2ul(startpos_);
    endpos = str2ul(endpos_);
    strand = strand_[0];
    pos1 = info_.find("GeneID:");
    if(pos1 != string::npos)
    {
        pos2 = info_.find(";",pos1);
        geneID = str2ul(info_.substr(pos1+7, pos2-pos1-7));
    }
    pos1 = info_.find("gene=");
    if(pos1 != string::npos)
    {
        pos2 = info_.find(";",pos1);
        geneName = info_.substr(pos1+5, pos2-pos1-5);
    }
}

CHumanGene::CHumanGene(int chrNum_, unsigned long pos_)
{
    chrNum = chrNum_;
    startpos = pos_;
}

CHumanGene::CHumanGene(int chrNum_, unsigned long startpos_, unsigned long endpos_)
{
    chrNum = chrNum_;
    startpos = startpos_;
    endpos = endpos_;
}

void CHumanGenes::LoadGenes(string path)
{
    cout << "Loading genes" << '\n';
    
    string line;
    int chrNum;
    
    ifstream f(path);
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }

    vector<string> flds;
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            if(line[0] == '#')
                continue;
            flds = split(line);
            if(flds[2] != "gene")
                continue;
            chrNum = CHumanGenome::GetChrFromNCBIName(flds[0]);
            if(chrNum < 0 || chrNum > 23)
                continue;
            genes.emplace(flds[0],flds[3],flds[4],flds[6],flds[8]);
        }
    }
    
    f.close();
}

void CHumanGenes::PrepareForSearch()
{
    set<CHumanGene>::iterator it;
    int prevchr;
    unsigned long maxend;

    prevchr = -1;
    maxend = 0;
    for(it=genes.begin();it!=genes.end();++it)
    {
        if(prevchr != (*it).chrNum)
        {
            (*it).maxendpos = 0;
            maxend = (*it).endpos;
        }
        else
        {
            (*it).maxendpos = (maxend < (*it).startpos) ? 0 : maxend;
            maxend = (maxend > (*it).endpos) ? maxend : (*it).endpos;
        }
        prevchr = (*it).chrNum;
    }
}

int CHumanGenes::GetGenesByPos(int chrNum, unsigned long pos, vector<CHumanGene>& geneList)
{
    set<CHumanGene>::iterator it;
    CHumanGene g(chrNum,pos);

    it = genes.upper_bound(g);
    if(it == genes.begin())
        return(0);
    it--;
    while(1)
    {
        if((*it).chrNum == chrNum && (*it).startpos <= pos && (*it).endpos >= pos)
            geneList.push_back((*it));
        if((*it).chrNum != chrNum || pos > (*it).maxendpos)
            break;
        it--;
    }
    
    if(geneList.empty())
        return(0);
    else
        return(1);
}

void CHumanGenes::MergeIntervals(vector<CHumanGene>& ret)
{
    set<CHumanGene>::iterator si;
    CHumanGene prevGene(-1,0);
    unsigned long startpos,endpos;

    for(si=genes.begin();si!=genes.end();si++)
    {
        if((*si).maxendpos == 0 || prevGene.chrNum != (*si).chrNum || prevGene.isNull())
        {
            if(!prevGene.isNull())
                ret.emplace_back(prevGene.chrNum,startpos,endpos);
            startpos = (*si).startpos;
            endpos = (*si).endpos;
        }
        else
        {
            endpos = ((*si).endpos > endpos ? (*si).endpos : endpos);
        }
        prevGene = (*si);
    }
    ret.emplace_back(prevGene.chrNum,startpos,endpos);
}

void CHumanGenes::SaveToFile(string path)
{
    ofstream f;
    f.open(path);
    set<CHumanGene>::iterator it;
    
    for(it=genes.begin();it!=genes.end();++it)
        f << (*it).chrNum << '\t' << (*it).startpos << '\t' << (*it).endpos << '\t' << (*it).maxendpos << '\t' << (*it).strand << '\t' << (*it).geneID << '\t' << (*it).geneName << '\n';
    
    f.close();
}
