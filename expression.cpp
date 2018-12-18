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


