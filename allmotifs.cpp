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

void CAllMotifs::AnalysisReplicationTiming(string motif, string dirpath)
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
    
//a.ClassifyMutations(&human);
    a.AnalyzeReplicationTiming(a.signatureMuts,dirpath+"/"+motif+"_3.txt");

    a.signatures.clear();
    a.signatureMuts.mutations.clear();
    a.signatureMuts.cancerSample.clear();
    a.otherMuts.mutations.clear();
    a.otherMuts.mutations.clear();
    
    for(it=nuc4cur.begin();it!=nuc4cur.end();it++)
    {
        a.signatures.push_back(CMutationSignature(motif,2,string(1,*it)));

//a.ClassifyMutations(&human);
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

void CAllMotifs::AnalysisRTExp(string dirpath, string cancer, string sample)
{
    vector<string> cancers;
    vector<string> samples;
    string renamedSample;
    
    //Load human genome
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    // Load mutations
    CMutations m;
    cancers.push_back(cancer);
    samples.push_back(sample);
    m.LoadMutations(CANCER_MUTATIONS_FILEFORMAT, CANCER_MUTATIONS_PATH, cancers, samples, 1, &human);
    m.RenameSamples("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_samples_list.csv",2,3,16);
    renamedSample = m.renamemap[sample];
    renamedSample = renamedSample.substr(0,16);
    
    // Load replication timing
    string path;
    CReplicationTiming rt;
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqHelas3WaveSignalRep1.mybed");
    rt.LoadReplicationTiming(path.c_str(), 0);
    rt.ReplicationStrand();
    
    rt.bins.push_back(CRTBin(0,-100,19.0528));
    rt.bins.push_back(CRTBin(1,19.0528,32.8395));
    rt.bins.push_back(CRTBin(2,32.8395,44.6798));
    rt.bins.push_back(CRTBin(3,44.6798,55.0382));
    rt.bins.push_back(CRTBin(4,55.0382,63.9244));
    rt.bins.push_back(CRTBin(5,63.9244,71.6868));
    rt.bins.push_back(CRTBin(6,71.6868,85.2286));
    
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();
    //genes.SaveToFile("/Users/mar/87/genes.txt");

    CExpression exp;
    exp.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_"+cancer+".txt", renamedSample);
    
    vector<CExpressionBin> expBins;
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));
    
    CSignatureAnalysis a;
    a.RTExpAllMotifs(m, dirpath, rt, exp, expBins, &human, genes, cancer, renamedSample);
    
}

