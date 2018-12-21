//
//  apobec.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <iostream>

#include "apobec.hpp"
#include <map>
#include <fstream>
#include "options.h"

CResultsKey::CResultsKey(string cancer_, string sample_, int bin_)
{
    cancer = cancer_;
    sample = sample_;
    bin = bin_;
}


CResultsValue::CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_)
{
    mutCnt = mutCnt_;
    leadingCnt = leadingCnt_;
    laggingCnt = laggingCnt_;
}

CResultsValue::CResultsValue(unsigned long mutCnt_, unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll_, unsigned long minusStrandAll_)
{
    mutCnt = mutCnt_;
    plusStrandConsistent = plusStrandConsistent_;
    minusStrandConsistent = minusStrandConsistent_;
    plusStrandAll = plusStrandAll_;
    minusStrandAll = minusStrandAll_;
}

void CAPOBEC::ClassifyMutations()
{
    // Load genome
    CHumanGenome human;
    human.InitializeHuman("37", "/Users/mar/BIO/BIODATA/HUMAN/CH37/hs_ref_GRCh37.p5_chr", ".fa", "FASTA");
    
    // Load mutations
    CMutations m;
    int isHeader = 1;
    m.LoadMutations("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations.tsv", isHeader);

    // Filter APOBEC mutations
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
    m.FilterMutations(apobecMuts,signatures,human,cancers,samples,&otherMuts);
    //apobecMuts.SaveToFile("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations_apobec.tsv");
    
    // Free human genome
    for(int i=0;i<human.chrCnt;i++)
        delete human.dna[i];
    delete human.dna;
}


void CAPOBEC::AnalysisReplicationTiming(CMutations& muts, string resultsFilename)
{
    
    // Load replication timing
    CReplicationTiming rtIMR90;
    CReplicationTiming rtMCF7;
    CReplicationTiming rtNHEK;
    rtIMR90.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed", 0);
    rtIMR90.ReplicationStrand();
    rtMCF7.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed", 0);
    rtMCF7.ReplicationStrand();
    rtNHEK.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed", 0);
    rtNHEK.ReplicationStrand();

    //
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    CReplicationTime rt;
    CMutation mut;
    CResultsValue rv;
    int bin;
    int strand;
    int res;
    
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
    
    ofstream f1;
    string path1;
    path1 = string(RESULTS_FOLDER)+string("/1.txt");
    f1.open(path1.c_str());

    
    
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];
                
        if (string(mut.cancer) == "LUAD" || string(mut.cancer) == "LUSC")
        {
            res = rtIMR90.GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
            if(res)
                bin = rtIMR90.GetRTBin(rt.RTvalue, binsIMR90);
            else
                bin = RT_NULLBIN_NOVALUE;
        }
        else if (string(mut.cancer) == "BLCA" || string(mut.cancer) == "HNSC")
        {
            res = rtNHEK.GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
            if(res)
                bin = rtNHEK.GetRTBin(rt.RTvalue, binsNHEK);
            else
                bin = RT_NULLBIN_NOVALUE;
        }
        else if (string(mut.cancer) == "BRCA")
        {
            res = rtMCF7.GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
            if(res)
                bin = rtMCF7.GetRTBin(rt.RTvalue, binsMCF7);
            else
                bin = RT_NULLBIN_NOVALUE;
        }

        if((mut.isForwardStrand == 1 && rt.isForward == 1) || (mut.isForwardStrand == 0 && rt.isForward == 0))
            strand = STRAND_LAGGING;
        else if((mut.isForwardStrand == 1 && rt.isForward == 0) || (mut.isForwardStrand == 0 && rt.isForward == 1))
            strand = STRAND_LEADING;
            
        it = results.find(CResultsKey(string(mut.cancer),string(mut.sample),bin));
        if ( it == results.end())
        {
            if (strand == STRAND_LEADING)
                rv = CResultsValue(1,1,0);
            else if (strand == STRAND_LAGGING)
                rv = CResultsValue(1,0,1);
            results.insert(pair<CResultsKey, CResultsValue>(CResultsKey(string(mut.cancer),string(mut.sample),bin), rv));
        }
        else
        {
            it->second.mutCnt++;
            if (strand == STRAND_LEADING)
                it->second.leadingCnt++;
            else if (strand == STRAND_LAGGING)
                it->second.laggingCnt++;
        }
        
        
    }
    
    ofstream f;
    string path;
    path = string(RESULTS_FOLDER)+string("/")+resultsFilename;
    f.open(path.c_str());
    
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.bin << '\t' << it->second.mutCnt << '\t' << it->second.leadingCnt << '\t' << it->second.laggingCnt << '\n';

    f.close();
    
}

int CAPOBEC::GetExpressionBin(string sample, string chr, unsigned long pos, char isForwardMut, CHumanGenes& genes, CExpression& exp, vector<CExpressionBin>& expBins, int& strand, int& strandInconsistence)
{
    vector<CHumanGene> geneList;
    double expValue, maxexp;
    int strandplus=0,strandminus=0;
    int res;
    int ret=0;
    int bin;
    
    genes.GetGenesByPos(CHumanGenome::GetChrNum(string(chr)), pos, geneList);
    if(geneList.empty())
    {
        strand = -1;
        strandInconsistence = 0;
        return(-2); // mutation not in genes
    }
    maxexp = -100000.0;
    ret = 0;
    strand = -1;
    for(int i=0;i<geneList.size();i++)
    {
        if((geneList[i].strand == '+' && isForwardMut == 1) || (geneList[i].strand == '-' && isForwardMut == 0))
            strandplus++;
        else if((geneList[i].strand == '-' && isForwardMut == 1) || (geneList[i].strand == '+' && isForwardMut == 0))
            strandminus--;
        res = exp.GetExpression(geneList[i].geneID, sample, expValue);
        if(res)
        {
            if(maxexp < expValue)
            {
                maxexp = expValue;
                if((geneList[i].strand == '+' && isForwardMut == 1) || (geneList[i].strand == '-' && isForwardMut == 0))
                    strand = 1;
                else if((geneList[i].strand == '-' && isForwardMut == 1) || (geneList[i].strand == '+' && isForwardMut == 0))
                    strand = 0;
            }
            ret = 1;
        }
    }
    
    if(strandplus !=0 && strandminus !=0)
        strandInconsistence = 1;
    else
        strandInconsistence = 0;
    
    if(ret == 0)
        return(-1); // mutation in genes, but no expression data
    
    bin = exp.GetExpressionBin(maxexp, expBins);
    return(bin);
}

void CAPOBEC::AnalysisExpression()
{
    
    CHumanGenes genes;
    genes.LoadGenes("/Users/mar/BIO/BIODATA/HUMAN/CH37/ref_GRCh37.p5_top_level.gff3");
    genes.PrepareForSearch();
    
    CExpression blca,brca,hnsc,luad,lusc;
    
    blca.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_BLCA.txt");
    brca.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_BRCA.txt");
    hnsc.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_HNSC.txt");
    luad.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_LUAD.txt");
    lusc.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_LUSC.txt");

    vector<CExpressionBin> expBins;
    
    expBins.emplace_back(0,-9999999,0);
    expBins.emplace_back(1,0,25);
    expBins.emplace_back(2,25,100);
    expBins.emplace_back(3,100,300);
    expBins.emplace_back(4,300,550);
    expBins.emplace_back(5,550,1000);
    expBins.emplace_back(6,1000,2000);
    expBins.emplace_back(7,2000,99999999999);
    
    
    map<CResultsKey, CResultsValue> resultsAPOBEC;
    map<CResultsKey, CResultsValue>::iterator it;
    CResultsValue rv;
    CMutation mut;
    int strand, strandInconsistence;
    double maxexp;
    int bin;
    for(int i=0;i<apobecMuts.mutations.size();i++)
    {
        mut = apobecMuts.mutations[i];
        if (string(mut.cancer) == "BLCA")
        {
             bin = GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, blca, expBins, strand, strandInconsistence);
        }
        else if (string(mut.cancer) == "BRCA")
        {
            bin = GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, brca, expBins, strand, strandInconsistence);
        }
        else if (string(mut.cancer) == "HNSC")
        {
            bin = GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, hnsc, expBins, strand, strandInconsistence);
        }
        else if (string(mut.cancer) == "LUAD")
        {
            bin = GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, luad, expBins, strand, strandInconsistence);
        }
        else if (string(mut.cancer) == "LUSC")
        {
            bin = GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, lusc, expBins, strand, strandInconsistence);
        }
        
        it = resultsAPOBEC.find(CResultsKey(string(mut.cancer),string(mut.sample),bin));
        if ( it == resultsAPOBEC.end())
        {
            rv = CResultsValue(1,0,0,0,0);
            resultsAPOBEC.insert(pair<CResultsKey, CResultsValue>(CResultsKey(string(mut.cancer),string(mut.sample),bin), rv));
        }
        else
        {
            it->second.mutCnt++;
            if(strand == 1)
            {
                it->second.plusStrandAll++;
                if(strandInconsistence == 0)
                    it->second.plusStrandConsistent++;
            }
            else if(strand == 0)
            {
                it->second.minusStrandAll++;
                if(strandInconsistence == 0)
                    it->second.minusStrandConsistent++;
            }
        }
    }
    
    ofstream f;
    f.open("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/results_expression_apobec.txt");
    
    for(it=resultsAPOBEC.begin(); it!=resultsAPOBEC.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.bin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
    f.close();
}

void CAPOBEC::CalculateTargetsinRTBins()
{
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
    
    CReplicationTiming rtIMR90;
    CReplicationTiming rtMCF7;
    CReplicationTiming rtNHEK;
    rtIMR90.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed", 0);
    rtMCF7.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed", 0);
    rtNHEK.LoadReplicationTiming("/Users/mar/BIO/BIODATA/ReplicationTiming/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed", 0);

    vector<string> motifs;
    motifs.push_back("TCA");
    motifs.push_back("TCT");
 
    rtIMR90.CalculateMotifinRTBins(binsIMR90, motifs, "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/TCW_in_RTbins_IMR90.txt");
    rtMCF7.CalculateMotifinRTBins(binsMCF7, motifs, "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/TCW_in_RTbins_MCF7.txt");
    rtNHEK.CalculateMotifinRTBins(binsNHEK, motifs, "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/TCW_in_RTbins_NHEK.txt");

    motifs.clear();
    motifs.push_back("X");
    
    rtIMR90.CalculateMotifinRTBins(binsIMR90, motifs, "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/ALL_in_RTbins_IMR90.txt");
    rtMCF7.CalculateMotifinRTBins(binsMCF7, motifs, "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/ALL_in_RTbins_MCF7.txt");
    rtNHEK.CalculateMotifinRTBins(binsNHEK, motifs, "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/ALL_in_RTbins_NHEK.txt");

}

