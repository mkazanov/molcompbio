//
//  apobec.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <iostream>

#include "apobec.hpp"
#include "ghuman.hpp"
#include "mutation.hpp"
#include "replicationtime.hpp"
#include <map>
#include <fstream>

CResultsKey::CResultsKey(string cancer_, string sample_, int RTbin_)
{
    cancer = cancer_;
    sample = sample_;
    RTbin = RTbin_;
}

void AnalysisReplicationTiming()
{
    
    // Load genome
    CHumanGenome human;
    human.InitializeHuman("37", "/Users/mar/BIO/BIODATA/HUMAN/CH37/hs_ref_GRCh37.p5_chr", ".fa", "FASTA");

    // Load mutations
    CMutations m;
    int isHeader = 1;
    m.LoadMutations("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations.tsv", isHeader);
    
    // Filter APOBEC mutations
    CMutations apobecMuts;
    vector<CMutationSignature> signatures;
    signatures.emplace_back("TCA",2,"T");
    signatures.emplace_back("TCT",2,"T");
    signatures.emplace_back("TCA",2,"G");
    signatures.emplace_back("TCT",2,"G");

    set<string> cancers;
    set<string> samples;
    cancers.insert("BLCA");
    cancers.insert("BRCA");
    cancers.insert("HNSC");
    cancers.insert("LUAD");
    cancers.insert("LUSC");
    m.FilterMutations(apobecMuts,signatures,human,cancers,samples);
    apobecMuts.SaveToFile("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations_apobec.tsv");
    
    // Load replication timing
    CReplicationTiming rtIMR90;
    CReplicationTiming rtMCF7;
    CReplicationTiming rtNHEK;
    rtIMR90.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed", 0);
    rtMCF7.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed", 0);
    rtNHEK.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed", 0);

    //
    map<CResultsKey, unsigned long> resultsAPOBEC;
    map<CResultsKey, unsigned long>::iterator it;
    CMutation mut;
    int bin;
    
    vector<CRTBin> binsIMR90, binsNHEK, binsMCF7;
    binsIMR90.emplace_back(1,-100,13.0766);
    binsIMR90.emplace_back(2,13.0766,28.3851);
    binsIMR90.emplace_back(3,28.3851,40.2474);
    binsIMR90.emplace_back(4,40.2474,52.0254);
    binsIMR90.emplace_back(5,52.0254,64.7194);
    binsIMR90.emplace_back(6,64.7194,75.0635);
    binsIMR90.emplace_back(7,75.0635,90.1735);
    
    binsMCF7.emplace_back(1,-100,19.8872);
    binsMCF7.emplace_back(2,19.8872,30.5993);
    binsMCF7.emplace_back(3,30.5993,40.0364);
    binsMCF7.emplace_back(4,40.0364,50.556);
    binsMCF7.emplace_back(5,50.556,60.3904);
    binsMCF7.emplace_back(6,60.3904,68.9953);
    binsMCF7.emplace_back(7,68.9953,86.4342);

    binsNHEK.emplace_back(1,-100,22.622);
    binsNHEK.emplace_back(2,-100,32.1929);
    binsNHEK.emplace_back(3,-100,41.1202);
    binsNHEK.emplace_back(4,-100,50.9531);
    binsNHEK.emplace_back(5,-100,59.6393);
    binsNHEK.emplace_back(6,-100,67.2079);
    binsNHEK.emplace_back(7,-100,80.7682);
    
    for(int i=0;i<apobecMuts.mutations.size();i++)
    {
        mut = apobecMuts.mutations[i];
        if (string(mut.cancer) == "LUAD" || string(mut.cancer) == "LUSC")
            bin = rtIMR90.GetRTBin(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, binsIMR90);
        else if (string(mut.cancer) == "BLCA" || string(mut.cancer) == "HNSC")
            bin = rtNHEK.GetRTBin(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, binsNHEK);
        else if (string(mut.cancer) == "BRCA")
            bin = rtMCF7.GetRTBin(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, binsMCF7);

        it = resultsAPOBEC.find(CResultsKey(string(mut.cancer),string(mut.sample),bin));
        if ( it == resultsAPOBEC.end())
            resultsAPOBEC.insert(pair<CResultsKey, unsigned long>(CResultsKey(string(mut.cancer),string(mut.sample),bin), 1));
        else
            it->second++;
    }
    
    ofstream f;
    f.open("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/results_apobec.txt");
    
    for(it=resultsAPOBEC.begin(); it!=resultsAPOBEC.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.RTbin << '\t' << it->second << '\n';

    f.close();
    
}
