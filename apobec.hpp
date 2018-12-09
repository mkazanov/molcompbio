//
//  apobec.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef apobec_hpp
#define apobec_hpp

#include <stdio.h>
#include <string>

using namespace std;

class CResultsKey {
public:
    string cancer;
    string sample;
    int RTbin;
    CResultsKey(string cancer_, string sample_, int RTbin_);
    bool operator< (const CResultsKey &right) const
    {
        if (cancer < right.cancer)
            return true;
        else if (cancer > right.cancer)
            return false;
        else
        {
            if (sample < right.sample)
                return true;
            else if (sample > right.sample)
                return false;
            else
            {
                if(RTbin < right.RTbin)
                    return true;
                else
                    return false;
            }
        }
    }

};

class CResultsValue{
public:
    unsigned long mutCnt;
    unsigned long leadingCnt;
    unsigned long laggingCnt;
    CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_);
    CResultsValue(){};
};

void AnalysisReplicationTiming();

#endif /* apobec_hpp */
