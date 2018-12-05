//
//  dnatest.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 04/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <iostream>
#include "dnatest.hpp"
#include <random>
#include "ghuman.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>

using namespace std;
using namespace std::chrono;

template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

void GenerateRandomDNAPositions()
{
    CHumanGenome g;
    g.InitializeHuman("37", "/Users/mar/BIO/BIODATA/HUMAN/CH37/hs_ref_GRCh37.p5_chr", ".fa", "FASTA", 0);

    random_device rd; // obtain a random number from hardware
    mt19937 eng(rd()); // seed the generator
    uniform_int_distribution<unsigned long> distr(1, 3095677412); // define the range
    unsigned long num;
    unsigned long curlen;
    int chrLen;
    int chr,pos;

    ofstream f;
    f.open("/Users/mar/4/dna_random_positions.txt");
    
    for(int i=0; i<30000000; i++)
    {
        num = distr(eng);
        curlen = 0;
        for(int j=0;j<24;j++)
        {
            chrLen = g.chrLen[j];
            if (num <= (curlen + chrLen))
            {
                pos = (int) (num - curlen);
                chr = j+1;
                f << chr << '\t' << pos << '\n';
                break;
            }
            curlen = curlen + chrLen;
        }
    }
    
    f.close();
}

void ReadRandomDNAPositions()
{
    const int FILE_SIZE = 30000000;
    vector<string> flds;
    int i;
    string line;
    
    struct dnaPosType{
        int chr;
        int pos;
    };
    struct dnaPosType* poss;
    
    CHumanGenome g;
    g.InitializeHuman("37", "/Users/mar/BIO/BIODATA/HUMAN/CH37/hs_ref_GRCh37.p5_chr", ".fa", "FASTA");

    ifstream f;
    f.open("/Users/mar/BIO/PROJECTS/CACHE/DnaTestFiles/dna_random_positions_30M_3.txt");
    
    poss = new dnaPosType[FILE_SIZE];
    
    i = 0;
    while(getline(f, line))
    {
        flds = split(line, '\t');
        poss[i].chr = stoi(flds[0]);
        poss[i].pos = stoi(flds[1]);
        i++;
    }
    f.close();
    
    printf("Start test\n");
    
    char bs;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(i=0;i<FILE_SIZE;i++)
    {
        bs = g.dna[poss[i].chr-1][poss[i].pos-1];
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    
    printf("End test: %i\n", duration);
}

