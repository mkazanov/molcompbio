//
//  replicationtime.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef replicationtime_hpp
#define replicationtime_hpp

#define RT_NULL 999999

#include <set>
#include <string>

using namespace std;

class CReplicationTime {
public:
    int chrNum;
    unsigned long startpos;
    unsigned long endpos;
    double RTvalue;
    CReplicationTime(string chr_, string startpos_, string endpos_, string RTvalue_);
    CReplicationTime(int chrNum_, unsigned long pos_);
    CReplicationTime(){};
    bool operator< (const CReplicationTime &right) const
    {
        if (chrNum < right.chrNum)
            return true;
        else if (chrNum == right.chrNum)
            return startpos < right.startpos;
        else
            return false;
    }
};

class CRTBin {
public:
    int binNum;
    double RTleft;
    double RTright;
    CRTBin(int binNum_, double RTleft_, double RTright_);
};

class CReplicationTiming {
    set<CReplicationTime> RTs;
public:
    void LoadReplicationTiming(string path, int isHeader);
    double GetRT(int chrNum, unsigned long pos);
    int GetRTBin(int chrNum, unsigned long pos, vector<CRTBin> bins);
};

#endif /* replicationtime_hpp */
