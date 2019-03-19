//
//  expression.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 16/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "expression.hpp"
#include <fstream>
#include <vector>
#include "service.h"

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


void CExpression::LoadExpression(string path)
{
    string line;
    string geneId,sample;
    
    ifstream f(path);
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }
    
    vector<string> flds,flds1;
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            flds1 = splitd(flds[0],'|');
            geneId = flds1[1];
            sample = flds[1].substr(0,16);
            data.insert(pair<CExpressionKey,double>(CExpressionKey(geneId,sample),str2d(flds[2])));
        }
    }
}

int CExpression::GetExpression(unsigned long geneId, string sample, double& expressionValue)
{
    map<CExpressionKey,double>::iterator it;
    
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

int CExpression::GetExpressionBin(double expValue, vector<CExpressionBin> bins)
{
    for(int i=0;i<bins.size();i++)
        if(expValue >= bins[i].expressionLeft && expValue < bins[i].expressionRight)
            return(bins[i].binNum);
    return(-1);
}

/*

int CReplicationTiming::CalculateMotifinRTBins(vector<CExpressionBin> bins, vector<string> motifs, string OUT_PATH)
{
    int bin;
    map<int,unsigned long> results;
    map<int,unsigned long>::iterator it;
    vector<string> motifsall;
    unsigned long motiflen,motifhalf;
    int break2;
    string motif;
    
    if(motifs.empty())
    {
        cerr << "Error: motifs array is empty." << '\n';
        return(0);
    }
    
    motiflen = motifs[0].length();
    for(int i=1;i<motifs.size();i++)
        if(motiflen != motifs[i].length())
        {
            cerr << "Error: motifs have different length" << '\n';
            return(0);
        }
    motifhalf = motiflen / 2;
    
    motifsall = motifs;
    for(int i=0;i<motifs.size();i++)
    {
        motif = CDNA::cDNA(motifs[i]);
        if(find(motifs.begin(),motifs.end(),motif) == motifs.end())
            motifsall.push_back(motif);
    }
    
    results[RT_NULLBIN_NOVALUE] = 0;
    results[RT_NULLBIN_NOBIN] = 0;
    for(int i=0;i<bins.size();i++)
        results[bins[i].binNum] = 0;
    
    // Load genome
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    for(int i=0;i<human.chrCnt;i++)
    {
        for(unsigned long j=motifhalf;j<(human.chrLen[i]-(motiflen-motifhalf-1));j++)
        {
            for(int k=0;k<motifsall.size();k++)
            {
                break2 = 0;
                for(int n=0;n<motiflen;n++)
                {
                    if(motifsall[k][n] == 'X')
                        continue;
                    if(human.dna[i][j-motifhalf+n] != motifsall[k][n])
                    {
                        break2 = 1;
                        break;
                    }
                }
                if(break2)
                    continue;
                bin = GetRTBin(i,j+1,bins);
                results[bin]++;
            }
        }
        cout << "chromosome:" << i << '\n';
    }
    
    ofstream f;
    f.open(OUT_PATH);
    
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first << '\t' << it->second <<'\n';
    
    f.close();
    
    return(1);
}

*/
