//
//  RTexpression.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 11/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#ifndef RTexpression_hpp
#define RTexpression_hpp

#include <stdio.h>
#include <map>
#include <vector>
#include <set>
#include "expression.hpp"
#include "replicationtime.hpp"

using namespace std;

class CRTexpression{
public:
    map<pair<int,int>,unsigned long> CalculateMotifsRTexpressionBins(CReplicationTiming& rt, CExpression exp, vector<CRTBin> rtBins,vector<CExpressionBin> expBins, CHumanGenes& genes, set<string> motifs, string sample, CHumanGenome* phuman);
};

#endif /* RTexpression_hpp */
