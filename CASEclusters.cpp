//
//  CASEclusters.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 07/02/2020.
//  Copyright © 2020 Marat Kazanov. All rights reserved.
//

#include "CASEclusters.hpp"
#include <map>
#include "mutation.hpp"
#include "replicationtime.hpp"
#include "options.h"
#include "expression.hpp"
#include "signanalysis.hpp"

void ProcessClusters()
{
    // Genome
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");

    // Replication timing data
    string path;
    //CReplicationTiming rtIMR90;
    //CReplicationTiming rtMCF7;
    //CReplicationTiming rtNHEK;
    CReplicationTiming rtHELA3;
    //path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed");
    //rtIMR90.LoadReplicationTiming(path.c_str(), 0);
    //rtIMR90.ReplicationStrand();
    //path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed");
    //rtMCF7.LoadReplicationTiming(path.c_str(), 0);
    //rtMCF7.ReplicationStrand();
    //path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed");
    //rtNHEK.LoadReplicationTiming(path.c_str(), 0);
    //rtNHEK.ReplicationStrand();
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqHelas3WaveSignalRep1.mybed");
    rtHELA3.LoadReplicationTiming(path.c_str(), 0);
    rtHELA3.ReplicationStrand();

    map<string,CReplicationTiming*> rtmap;
    //rtmap["BLCA"] = &rtNHEK;
    //rtmap["BRCA"] = &rtMCF7;
    //rtmap["HNSC"] = &rtNHEK;
    //rtmap["LUAD"] = &rtIMR90;
    //rtmap["LUSC"] = &rtIMR90;
    rtmap["CESC"] = &rtHELA3;

    //rtIMR90.bins.push_back(CRTBin(0,-100,13.0766));
    //rtIMR90.bins.push_back(CRTBin(1,13.0766,28.3851));
    //rtIMR90.bins.push_back(CRTBin(2,28.3851,40.2474));
    //rtIMR90.bins.push_back(CRTBin(3,40.2474,52.0254));
    //rtIMR90.bins.push_back(CRTBin(4,52.0254,64.7194));
    //rtIMR90.bins.push_back(CRTBin(5,64.7194,75.0635));
    //rtIMR90.bins.push_back(CRTBin(6,75.0635,90.1735));
    
    //rtMCF7.bins.push_back(CRTBin(0,-100,19.8872));
    //rtMCF7.bins.push_back(CRTBin(1,19.8872,30.5993));
    //rtMCF7.bins.push_back(CRTBin(2,30.5993,40.0364));
    //rtMCF7.bins.push_back(CRTBin(3,40.0364,50.556));
    //rtMCF7.bins.push_back(CRTBin(4,50.556,60.3904));
    //rtMCF7.bins.push_back(CRTBin(5,60.3904,68.9953));
    //rtMCF7.bins.push_back(CRTBin(6,68.9953,86.4342));
    
    //rtNHEK.bins.push_back(CRTBin(0,-100,22.622));
    //rtNHEK.bins.push_back(CRTBin(1,-100,32.1929));
    //rtNHEK.bins.push_back(CRTBin(2,-100,41.1202));
    //rtNHEK.bins.push_back(CRTBin(3,-100,50.9531));
    //rtNHEK.bins.push_back(CRTBin(4,-100,59.6393));
    //rtNHEK.bins.push_back(CRTBin(5,-100,67.2079));
    //rtNHEK.bins.push_back(CRTBin(6,-100,80.7682));
    
    rtHELA3.bins.push_back(CRTBin(0,-100,19.0528));
    rtHELA3.bins.push_back(CRTBin(1,19.0528,32.8395));
    rtHELA3.bins.push_back(CRTBin(2,32.8395,44.6798));
    rtHELA3.bins.push_back(CRTBin(3,44.6798,55.0382));
    rtHELA3.bins.push_back(CRTBin(4,55.0382,63.9244));
    rtHELA3.bins.push_back(CRTBin(5,63.9244,71.6868));
    rtHELA3.bins.push_back(CRTBin(6,71.6868,85.2286));

    // Expression data
    map<string,CExpression*> expmap;
    CExpression blca,brca,hnsc,luad,lusc,cesc;
    
    //blca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BLCA.txt");
    //brca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BRCA.txt");
    //hnsc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_HNSC.txt");
    //luad.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUAD.txt");
    //lusc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUSC.txt");
    cesc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_CESC.txt");
    //expmap["BLCA"] = &blca;
    //expmap["BRCA"] = &brca;
    //expmap["HNSC"] = &hnsc;
    //expmap["LUAD"] = &luad;
    //expmap["LUSC"] = &lusc;
    expmap["CESC"] = &cesc;

    vector<CExpressionBin> expBins;
    
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));
   
    // Load genes
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();

    // Load mutations
    CMutations m;
    vector<string> cancers;
    vector<string> samples;
    m.LoadMutations(FILE_FORMAT_FRIDRIKSSON,
                    "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/cluster/cluster_muts_CESC.csv",
                    cancers,samples,1,&human);
    m.RenameSamples("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_samples_list.csv",2,3,16);
  
    CMutation* mut;
    int bin;
    int strand;
    int strandInconsistence;
    for(int i=0;i<m.mutations.size();i++)
    {
        mut = &m.mutations[i];
     
        rtmap[string(mut->cancer)]->GetRTBinStrand(rtmap,(*mut),bin,strand);
        mut->RTbin = bin;
        mut->RTstrand = strand;
        
        bin = (*expmap[string(mut->cancer)]).GetExpressionBin(string(mut->sample), mut->chr, mut->pos, mut->isForwardStrand, genes, expBins, strand, strandInconsistence);
        mut->EXPbin = bin;
        mut->EXPstrand = strand;
        
    }
    
    m.SaveToFileRTExp("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/cluster/cluster_muts_processed_CESC.csv");
    
}

void MakeApobecMutationsFile()
{
    // Genome
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    
    // Load mutations
    CMutations m;
    vector<string> cancers;
    vector<string> samples;
    m.LoadMutations(FILE_FORMAT_FRIDRIKSSON,
                    "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/cluster/cluster_muts_CESC.csv",
                    cancers,samples,1,&human);
    m.RenameSamples("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_samples_list.csv",2,3,16);

     CSignatureAnalysis a;
     
     a.signatures.push_back(CMutationSignature("TCA",2,"T"));
     a.signatures.push_back(CMutationSignature("TCT",2,"T"));
     a.signatures.push_back(CMutationSignature("TCA",2,"G"));
     a.signatures.push_back(CMutationSignature("TCT",2,"G"));
     a.signatures.push_back(CMutationSignature("TCC",2,"T"));
     a.signatures.push_back(CMutationSignature("TCG",2,"T"));
     a.signatures.push_back(CMutationSignature("TCC",2,"G"));
     a.signatures.push_back(CMutationSignature("TCG",2,"G"));
     
     m.FilterMutations(a.signatureMuts,a.signatures,human,cancers,samples,&a.otherMuts);
     a.signatureMuts.SaveToFile("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/cluster/mutations_apobec_CESC.tsv");

}
