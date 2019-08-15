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

void CAllMotifs::AnalysisRTExp(string dirpath, string cancer, string sample)
{
    vector<string> cancers;
    vector<string> samples;
    
    //Load human genome
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    // Load mutations
    CMutations m;
    int isHeader = 1;
    cancers.push_back(cancer);
    samples.push_back(sample);
    m.LoadMutations(CANCER_MUTATIONS, isHeader, cancers, samples, 1, &human);
    
    // Load replication timing
    string path;
    CReplicationTiming rtIMR90;
    CReplicationTiming rtMCF7;
    CReplicationTiming rtNHEK;
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed");
    rtIMR90.LoadReplicationTiming(path.c_str(), 0);
    rtIMR90.ReplicationStrand();
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed");
    rtMCF7.LoadReplicationTiming(path.c_str(), 0);
    rtMCF7.ReplicationStrand();
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed");
    rtNHEK.LoadReplicationTiming(path.c_str(), 0);
    rtNHEK.ReplicationStrand();
    
    rtIMR90.bins.push_back(CRTBin(0,-100,13.0766));
    rtIMR90.bins.push_back(CRTBin(1,13.0766,28.3851));
    rtIMR90.bins.push_back(CRTBin(2,28.3851,40.2474));
    rtIMR90.bins.push_back(CRTBin(3,40.2474,52.0254));
    rtIMR90.bins.push_back(CRTBin(4,52.0254,64.7194));
    rtIMR90.bins.push_back(CRTBin(5,64.7194,75.0635));
    rtIMR90.bins.push_back(CRTBin(6,75.0635,90.1735));
    
    rtMCF7.bins.push_back(CRTBin(0,-100,19.8872));
    rtMCF7.bins.push_back(CRTBin(1,19.8872,30.5993));
    rtMCF7.bins.push_back(CRTBin(2,30.5993,40.0364));
    rtMCF7.bins.push_back(CRTBin(3,40.0364,50.556));
    rtMCF7.bins.push_back(CRTBin(4,50.556,60.3904));
    rtMCF7.bins.push_back(CRTBin(5,60.3904,68.9953));
    rtMCF7.bins.push_back(CRTBin(6,68.9953,86.4342));
    
    rtNHEK.bins.push_back(CRTBin(0,-100,22.622));
    rtNHEK.bins.push_back(CRTBin(1,-100,32.1929));
    rtNHEK.bins.push_back(CRTBin(2,-100,41.1202));
    rtNHEK.bins.push_back(CRTBin(3,-100,50.9531));
    rtNHEK.bins.push_back(CRTBin(4,-100,59.6393));
    rtNHEK.bins.push_back(CRTBin(5,-100,67.2079));
    rtNHEK.bins.push_back(CRTBin(6,-100,80.7682));
    
    map<string,CReplicationTiming*> rtmap;
    rtmap["BLCA"] = &rtNHEK;
    rtmap["BRCA"] = &rtMCF7;
    rtmap["HNSC"] = &rtNHEK;
    rtmap["LUAD"] = &rtIMR90;
    rtmap["LUSC"] = &rtIMR90;
    
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();
    //genes.SaveToFile("/Users/mar/87/genes.txt");

    CExpression exp;
    exp.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_"+cancer+".txt", sample);
    
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
    a.RTExpAllMotifs(m, dirpath, rtmap, exp, expBins, &human, genes, cancer, sample);
    
}

