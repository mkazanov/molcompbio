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

void CAPOBEC::ClassifyMutations(CHumanGenome* phuman_)
{
    CHumanGenome* phuman;
    
    // Load genome
    if(phuman_)
        phuman = phuman_;
    else
    {
        phuman = new CHumanGenome();
        phuman->InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    }
        
    // Load mutations
    CMutations m;
    int isHeader = 1;
    m.LoadMutations(CANCER_MUTATIONS, isHeader);

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
    m.FilterMutations(apobecMuts,signatures,(*phuman),cancers,samples,&otherMuts);
    //apobecMuts.SaveToFile("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations_apobec.tsv");
    
    // Free human genome
    if(!phuman_)
    {
        for(int i=0;i<phuman->chrCnt;i++)
            delete phuman->dna[i];
        delete phuman->dna;
    }
}

void CAPOBEC::AnalysisReplicationTiming(CMutations& muts, string resultsFilename)
{
    
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
        else
            strand == STRAND_NULL;
            
        it = results.find(CResultsKey(string(mut.cancer),string(mut.sample),bin));
        if ( it == results.end())
        {
            if (strand == STRAND_LEADING)
                rv = CResultsValue(1,1,0);
            else if (strand == STRAND_LAGGING)
                rv = CResultsValue(1,0,1);
            else
                rv = CResultsValue(1,0,0);
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
    path = string(RESULTS_FOLDER)+string("/")+resultsFilename;
    f.open(path.c_str());
    
    f << "Cancer\tSample\tReplicationBin\tMutationCnt\tLeadingCnt\tLaggingCnt\n";
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
    int bin;
    int expFound;
    
    genes.GetGenesByPos(CHumanGenome::GetChrNum(string(chr)), pos, geneList);
    if(geneList.empty())
    {
        strand = -1;
        strandInconsistence = 0;
        return(-2); // mutation not in genes
    }
    maxexp = -100000.0;
    expFound = 0;
    strand = -1;
    for(int i=0;i<geneList.size();i++)
    {
        if((geneList[i].strand == '+' && isForwardMut == 1) || (geneList[i].strand == '-' && isForwardMut == 0))
            strandplus++;
        else if((geneList[i].strand == '-' && isForwardMut == 1) || (geneList[i].strand == '+' && isForwardMut == 0))
            strandminus++;
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
                else
                    strand = -1;
            }
            expFound = 1;
        }
    }
    
    if(strandplus !=0 && strandminus !=0)
        strandInconsistence = 1;
    else
        strandInconsistence = 0;
    
    if(expFound == 0)
        return(-1); // mutation in genes, but no expression data
    
    bin = exp.GetExpressionBin(maxexp, expBins);
    return(bin);
}

void CAPOBEC::AnalysisExpression(CMutations& muts, string resultsFilename)
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
    
    
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    pair<map<CResultsKey, CResultsValue>::iterator, bool> res;
    CResultsValue rv;
    CMutation mut;
    int strand, strandInconsistence;
    int bin;
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];
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
        
        it = results.find(CResultsKey(string(mut.cancer),string(mut.sample),bin));
        if ( it == results.end())
        {
            rv = CResultsValue(1,0,0,0,0);
            res = results.insert(pair<CResultsKey, CResultsValue>(CResultsKey(string(mut.cancer),string(mut.sample),bin), rv));
            it = res.first;
        }
        else
        {
            it->second.mutCnt++;
        }
        
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
    
    ofstream f;
    string path;
    path = string(RESULTS_FOLDER)+string("/")+resultsFilename;
    f.open(path.c_str());
    
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.bin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
    f.close();
}

void CAPOBEC::CalculateTargetsinRTBins(CHumanGenome* phuman)
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
    
    string path;
    CReplicationTiming rtIMR90;
    CReplicationTiming rtMCF7;
    CReplicationTiming rtNHEK;
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed");
    rtIMR90.LoadReplicationTiming(path.c_str(), 0);
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed");
    rtMCF7.LoadReplicationTiming(path.c_str(), 0);
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed");
    rtNHEK.LoadReplicationTiming(path.c_str(), 0);

    set<string> motifs;
    motifs.insert("TCA");
    motifs.insert("TCT");
 
    path = string(RESULTS_FOLDER)+string("/TCW_in_RTbins_IMR90.txt");
    rtIMR90.CalculateMotifinRTBins(binsIMR90, motifs, path.c_str(), phuman);
    path = string(RESULTS_FOLDER)+string("/TCW_in_RTbins_MCF7.txt");
    rtMCF7.CalculateMotifinRTBins(binsMCF7, motifs, path.c_str(), phuman);
    path = string(RESULTS_FOLDER)+string("/TCW_in_RTbins_NHEK.txt");
    rtNHEK.CalculateMotifinRTBins(binsNHEK, motifs, path.c_str(), phuman);

    motifs.clear();
    motifs.insert("X");
    
    path = string(RESULTS_FOLDER)+string("/ALL_in_RTbins_IMR90.txt");
    rtIMR90.CalculateMotifinRTBins(binsIMR90, motifs, path.c_str(), phuman);
    path = string(RESULTS_FOLDER)+string("/ALL_in_RTbins_MCF7.txt");
    rtMCF7.CalculateMotifinRTBins(binsMCF7, motifs, path.c_str(), phuman);
    path = string(RESULTS_FOLDER)+string("/ALL_in_RTbins_NHEK.txt");
    rtNHEK.CalculateMotifinRTBins(binsNHEK, motifs, path.c_str(), phuman);

}

/*
void CAPOBEC::CalculateTargetsinExpressionBins()
{
    vector<CExpressionBin> expBins;
    
    expBins.emplace_back(0,-9999999,0);
    expBins.emplace_back(1,0,25);
    expBins.emplace_back(2,25,100);
    expBins.emplace_back(3,100,300);
    expBins.emplace_back(4,300,550);
    expBins.emplace_back(5,550,1000);
    expBins.emplace_back(6,1000,2000);
    expBins.emplace_back(7,2000,99999999999);

}

*/

void CAPOBEC::CalculateAPOBECEnrichment(CHumanGenome* phuman_)
{
    CHumanGenome* phuman;
    CMutationSignature s;
    unsigned long cytosineCnt, TCWcnt;
    set<string> motifs;
    ofstream f;
    string path;
    
    class CMapKey {
    public:
        string cancer;
        string sample;
        CMapKey(string cancer_, string sample_)
        {
            cancer = cancer_;
            sample = sample_;
        }
        bool operator< (const CMapKey &right) const
        {
            if (cancer < right.cancer)
                return true;
            else if (cancer > right.cancer)
                return false;
            else
            {
                if (sample < right.sample)
                    return true;
                else 
                    return false;
            }
        }
    };
    
    class CMapValue {
    public:
        unsigned long APOBECmutsCnt;
        unsigned long cytosineMutsCnt;
        double enrichment;
        double enrichmentExcludeTCW;
        CMapValue(unsigned long APOBECmutsCnt_, unsigned long cytosineMutsCnt_, double enrichment_, double enrichmentExcludeTCW_)
        {
            APOBECmutsCnt = APOBECmutsCnt_;
            cytosineMutsCnt = cytosineMutsCnt_;
            enrichment = enrichment_;
            enrichmentExcludeTCW = enrichmentExcludeTCW_;
        }
    };
    
    if(phuman_)
        phuman = phuman_;
    else
    {
        phuman = new CHumanGenome();
        phuman->InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    }
    
    motifs.insert("C");
    cytosineCnt = s.CountMotifGenome(motifs, phuman);
    
    motifs.clear();
    motifs.insert("TCT");
    motifs.insert("TCA");
    TCWcnt = s.CountMotifGenome(motifs, phuman);
    
    cout << "C: " << cytosineCnt << ", TCW: " << TCWcnt << '\n';
    
    if(otherMuts.mutations.size() == 0)
    {
        cerr << "Mutations are not classified." << '\n';
        exit(0);
    }
        
    map<CMapKey, CMapValue> enrichments;
    map<CMapKey, CMapValue>::iterator it;

    for(int i=0;i<apobecMuts.mutations.size();i++)
    {
        it = enrichments.find(CMapKey(apobecMuts.mutations[i].cancer,apobecMuts.mutations[i].sample));
        if(it == enrichments.end())
            enrichments.insert(pair<CMapKey, CMapValue>(CMapKey(apobecMuts.mutations[i].cancer, apobecMuts.mutations[i].sample), CMapValue(1,0,0,0)));
        else
            it->second.APOBECmutsCnt++;
    }
    
    for(int i=0;i<otherMuts.mutations.size();i++)
    {
        if(otherMuts.mutations[i].refallele != "C" && otherMuts.mutations[i].refallele != "G")
            continue;
        it = enrichments.find(CMapKey(otherMuts.mutations[i].cancer,otherMuts.mutations[i].sample));
        if(it == enrichments.end())
            enrichments.insert(pair<CMapKey, CMapValue>(CMapKey(otherMuts.mutations[i].cancer, otherMuts.mutations[i].sample), CMapValue(0,1,0,0)));
        else
            it->second.cytosineMutsCnt++;
    }
    
    path = string(RESULTS_FOLDER)+string("/enrichment.txt");
    f.open(path.c_str());
    double cytTCWratio;
    
    cytTCWratio = (double) cytosineCnt / TCWcnt;
    
    f << "cancer" << '\t' << "sample" << '\t' << "APOBEC_mutaions" << '\t' << "Cytosine_mutations" << '\t' << "Enrichment" << '\t' << "Enrichment_exclude_TCW" << '\n';
    for(it=enrichments.begin(); it!=enrichments.end(); ++it)
    {
        cout << it->first.cancer << '\t' << it->first.sample << '\t' << it->second.APOBECmutsCnt << '\t' << it->second.cytosineMutsCnt << '\t' << it->second.enrichment << '\n';
        it->second.enrichment = ((double)it->second.APOBECmutsCnt / ((double)it->second.cytosineMutsCnt + (double)it->second.APOBECmutsCnt)) * cytTCWratio;
        it->second.enrichmentExcludeTCW = ((double)it->second.APOBECmutsCnt / ((double)it->second.cytosineMutsCnt)) * cytTCWratio;
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->second.APOBECmutsCnt << '\t' << it->second.cytosineMutsCnt << '\t' << it->second.enrichment << '\t' << it->second.enrichmentExcludeTCW << '\n';
    }
        
    f.close();
    
    // Free human genome
    if(!phuman_)
    {
        for(int i=0;i<phuman->chrCnt;i++)
            delete phuman->dna[i];
        delete phuman->dna;
    }
}

