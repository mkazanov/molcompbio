//
//  mutation.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "mutation.hpp"
#include "dna.hpp"
#include "ghuman.hpp"
#include <fstream>
#include <iostream>
#include "service.h"
#include <chrono>
#include <thread>

CMutation::CMutation(string cancer_,
                     string sample_,
                     string chr_,
                     string pos_,
                     string refallele_,
                     string varallele_)
{
    cancer_.copy(cancer,STRLEN_CANCER);
    cancer[STRLEN_CANCER] = '\0';
    sample_.copy(sample,STRLEN_SAMPLE);
    sample[STRLEN_SAMPLE] = '\0';
    chr_.copy(chr,STRLEN_CHR);
    chr[STRLEN_CHR] = '\0';
    pos = str2ul(pos_.c_str());
    refallele = refallele_;
    varallele = varallele_;
}


void CMutations::LoadMutations(string path, int isHeader)
{
    string line;
    clock_t c1,c2;
    
    ifstream f(path);
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }

    printf("Calculating number of mutations ...\n");
    mutationsCnt = 0;
    while(1)
    {
        getline(f, line);
        if(f.eof())
            break;
        mutationsCnt ++;
    }
    f.clear();
    f.seekg(0, ios::beg);
    
    mutationsCnt -= isHeader;
    
    printf("Allocating space for mutaiton ...\n");
    c1 = clock();
    mutations.reserve(mutationsCnt);
    c2 = clock();
    printf("Executing time: %lu \n", c2 - c1);

    if(isHeader)
        getline(f, line);
    
    printf("Mutations loading ...\n");
    int i=0;
    vector<string> flds;
    c1 = clock();
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            mutations.emplace_back(flds[1],flds[0],flds[2],flds[3],flds[4],flds[5]);
            i++;
        }
    }
    c2 = clock();
    printf("Loaded %i mutations\n", mutationsCnt);
    printf("Executing time: %lu \n", c2 - c1);
    
}

void CMutations::FilterMutations(CMutations& filteredMutations, vector<CMutationSignature>& signatures, CHumanGenome& human,
                                 set<string> cancers, set<string> samples)
{
    int i,j;
    CMutation m;
    CMutationSignature s;
    string mutBase;
    string ss;
    unsigned long startPosShift, endPosShift;
    string gmotif;
    
    for(i=0;i<mutationsCnt;i++)
    {
        for(j=0;j<signatures.size();j++)
        {
            m = mutations[i];
            
            if (string(m.chr) == "M")
                continue;
            if (!cancers.empty() && cancers.find(string(m.cancer)) == cancers.end())
                continue;
            if (!samples.empty() && samples.find(string(m.sample)) == samples.end())
                continue;
            
            s = signatures[j];
            startPosShift = 1 - s.mutationPos;
            endPosShift = s.motif.length() - s.mutationPos;
            mutBase = s.motif[s.mutationPos-1];
            ss = string(m.chr);
            gmotif = human.dnaSubstr(human.GetChrNum(string(m.chr)),m.pos+startPosShift,m.pos+endPosShift);
                        
            if ((m.refallele == mutBase &&
                 (s.AnyNewBase() || (m.varallele == s.newbase) ) &&
                 gmotif == s.motif) ||
                (m.refallele == CDNA::cDNA(mutBase) &&
                 (s.AnyNewBase() || m.varallele == CDNA::cDNA(s.newbase) ) &&
                 gmotif == CDNA::cDNA(s.motif) )
                )
                {
                    filteredMutations.mutations.push_back(m);
                }
                
        }
    }
}

void CMutations::SaveToFile(string path)
{
    ofstream f;
    f.open(path);
    int i;
    
    for(i=0;i<mutations.size();i++)
        f << string(mutations[i].cancer) << '\t' << string(mutations[i].sample) << '\t' <<
        string(mutations[i].chr) << '\t' << ul2str(mutations[i].pos) << '\t' << mutations[i].refallele <<
        '\t' << mutations[i].varallele << '\n';
    
    f.close();
}
