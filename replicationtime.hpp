//
//  replicationtime.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef replicationtime_hpp
#define replicationtime_hpp

#include <set>
#include <string>

#define RT_NULL 999999
#define CHR_NULL -1
#define STRAND_LEADING 1
#define STRAND_LAGGING 0


using namespace std;

class CReplicationTime {
public:
    int chrNum;
    unsigned long startpos;
    unsigned long endpos;
    double RTvalue;
    mutable short isForward;
    CReplicationTime(string chr_, string startpos_, string endpos_, string RTvalue_);
    CReplicationTime(int chrNum_, unsigned long pos_);
    CReplicationTime(){};
    bool isRTnull();
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
    CReplicationTime GetRT(int chrNum, unsigned long pos);
    int GetRTBin(double RTvalue, vector<CRTBin> bins);
    int GetRTBin(CReplicationTime rt, vector<CRTBin> bins);
    void ReplicationStrand();
    void SaveToFile(string path);
};

#endif /* replicationtime_hpp */
