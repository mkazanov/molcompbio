//
//  allmotifs.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 14/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#include "allmotifs.hpp"
#include <set>
#include "dna.hpp"
#include <fstream>
#include "signanalysis.hpp"
#include "options.h"
#include <iostream>

void CAllMotifs::GenerateMotifsFile(string path, int includeComplimentary)
{
    int i,j,k;
    char nuc[4] = {'A','G','C','T'};
    string motif;
    set<string> motifs;
    set<string>::iterator it;
    
    for(i=0;i<4;i++)
        for(j=0;j<4;j++)
            for(k=0;k<4;k++)
            {
                motif = string(1,nuc[i])+string(1,nuc[j])+string(1,nuc[k]);
                if(motifs.find(CDNA::cDNA(motif)) != motifs.end() && includeComplimentary==0)
                    continue;
                motifs.insert(motif);
            }
    
    ofstream f;
    f.open(path.c_str());
    
    for(it=motifs.begin();it!=motifs.end();it++)
        f << (*it) << '\n';
    
    f.close();
}

void CAllMotifs::RunAnalysis(string motif, string dirpath)
{
    time_t t;
    tm* now;
    
    // Start date/time
    t = time(0);   // get time now
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");

    set<char> nuc4,nuc4cur;
    set<char>::iterator it;
    nuc4.insert('A');
    nuc4.insert('C');
    nuc4.insert('G');
    nuc4.insert('T');
    
    nuc4cur = nuc4;
    nuc4cur.erase(motif[1]);
    
    CSignatureAnalysis a;
    
    for(it=nuc4cur.begin();it!=nuc4cur.end();it++)
        a.signatures.push_back(CMutationSignature(motif,2,string(1,*it)));
    
    a.ClassifyMutations(&human);
    a.AnalyzeReplicationTiming(a.signatureMuts,dirpath+"/"+motif+"_3.txt");

    a.signatures.clear();
    a.signatureMuts.mutations.clear();
    a.signatureMuts.cancerSample.clear();
    a.otherMuts.mutations.clear();
    a.otherMuts.mutations.clear();
    
    for(it=nuc4cur.begin();it!=nuc4cur.end();it++)
    {
        a.signatures.push_back(CMutationSignature(motif,2,string(1,*it)));

        a.ClassifyMutations(&human);
        a.AnalyzeReplicationTiming(a.signatureMuts,dirpath+"/"+motif+"_"+string(1,*it)+".txt");
        
        a.signatures.clear();
        a.signatureMuts.mutations.clear();
        a.signatureMuts.cancerSample.clear();
        a.otherMuts.mutations.clear();
        a.otherMuts.mutations.clear();
    }
    
    a.signatures.push_back(CMutationSignature(motif,2,"X"));
    a.CalculateTargetsinRTBins(dirpath+"/"+motif,&human,1);

    // End date/time
    t = time(0);
    now = localtime(&t);
    cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << ' ' << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

}
