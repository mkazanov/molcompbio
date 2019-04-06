//
//  expression.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 16/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef expression_hpp
#define expression_hpp

#include <string>
#include <map>
#include <vector>

using namespace std;

class CExpressionKey
{
public:
    unsigned long geneId;
    string sample;
    CExpressionKey(string geneId_, string sample_);
    CExpressionKey(unsigned long geneId_, string sample_);
    bool operator< (const CExpressionKey &right) const
    {
        if (geneId < right.geneId)
            return true;
        else if (geneId > right.geneId)
            return false;
        else
        {
            if (sample < right.sample)
                return true;
            else if (sample >= right.sample)
                return false;
        }
        return false;
    }
};

class CExpressionBin {
public:
    int binNum;
    double expressionLeft;
    double expressionRight;
    CExpressionBin(int binNum_, double expressionLeft_, double expressionRight_);
};

class CExpression
{
public:
    map<CExpressionKey, double> data;
    void LoadExpression(string path);
    int GetExpression(unsigned long geneId, string sample, double& expressionValue);
    int GetExpressionBin(double expValue, vector<CExpressionBin> bins);
};
#endif /* expression_hpp */
