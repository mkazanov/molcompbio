  //
//  apobec.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <iostream>

#include "signanalysis.hpp"
#include <map>
#include <fstream>
#include "options.h"
#include <ctime>
#include <cstring>
#include "RTexpression.hpp"
#include "dna.hpp"

CResultsKey::CResultsKey(string motif_, string cancer_, string sample_, int rtbin_, int expbin_)
{
    motif = motif_;
    cancer = cancer_;
    sample = sample_;
    rtbin = rtbin_;
    expbin = expbin_;
}

CResultsKey::CResultsKey(string cancer_, string sample_, int rtbin_, int expbin_)
{
    motif = "X";
    cancer = cancer_;
    sample = sample_;
    rtbin = rtbin_;
    expbin = expbin_;
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

CResultsValue::CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_,unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll_, unsigned long minusStrandAll_)
{
    mutCnt = mutCnt_;
    leadingCnt = leadingCnt_;
    laggingCnt = laggingCnt_;
    plusStrandConsistent = plusStrandConsistent_;
    minusStrandConsistent = minusStrandConsistent_;
    plusStrandAll = plusStrandAll_;
    minusStrandAll = minusStrandAll_;
}

void CSignatureAnalysis:: ClassifyMutations(CHumanGenome* phuman_)
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
    
    set<string> cancers;
    set<string> samples;
    cancers.insert("BLCA");
    cancers.insert("BRCA");
    cancers.insert("HNSC");
    cancers.insert("LUAD");
    cancers.insert("LUSC");
    m.FilterMutations(signatureMuts,signatures,(*phuman),cancers,samples,&otherMuts);
    signatureMuts.SaveToFile("/Users/mar/63/mutations_apobec.tsv");
    
    signatureMuts.GetUniqueCancersSamples();
    
    // Free human genome
    if(!phuman_)
    {
        for(int i=0;i<phuman->chrCnt;i++)
            delete phuman->dna[i];
        delete phuman->dna;
    }
}

void CSignatureAnalysis::AnalyzeReplicationTiming(CMutations& muts, string resultsFilename)
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
    
    map<string,CReplicationTiming*> rtmap;
    rtmap["BLCA"] = &rtNHEK;
    rtmap["BRCA"] = &rtMCF7;
    rtmap["HNSC"] = &rtNHEK;
    rtmap["LUAD"] = &rtIMR90;
    rtmap["LUSC"] = &rtIMR90;

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
    
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];
                
        res = rtmap[string(mut.cancer)]->GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
        if(res)
            bin = rtmap[string(mut.cancer)]->GetRTBin(rt.RTvalue, (rtmap[string(mut.cancer)]->bins));
        else
            bin = RT_NULLBIN_NOVALUE;

        if((mut.isForwardStrand == 1 && rt.isForward == 1) || (mut.isForwardStrand == 0 && rt.isForward == 0))
            strand = STRAND_LAGGING;
        else if((mut.isForwardStrand == 1 && rt.isForward == 0) || (mut.isForwardStrand == 0 && rt.isForward == 1))
            strand = STRAND_LEADING;
        else
            strand = STRAND_NULL;
            
        it = results.find(CResultsKey(string(mut.cancer),string(mut.sample),bin,EXP_NULLBIN_NOEXPDATA));
        if ( it == results.end())
        {
            if (strand == STRAND_LEADING)
                rv = CResultsValue(1,1,0);
            else if (strand == STRAND_LAGGING)
                rv = CResultsValue(1,0,1);
            else
                rv = CResultsValue(1,0,0);
            results.insert(pair<CResultsKey, CResultsValue>(CResultsKey(string(mut.cancer),string(mut.sample),bin,EXP_NULLBIN_NOEXPDATA), rv));
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
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.rtbin << '\t' << it->second.mutCnt << '\t' << it->second.leadingCnt << '\t' << it->second.laggingCnt << '\n';

    f.close();
    
}

void CSignatureAnalysis::AnalyzeExpression(CMutations& muts, string resultsFilename)
{
    
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();
    
    map<string,CExpression*> expmap;
    CExpression blca,brca,hnsc,luad,lusc;
    
    blca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BLCA.txt");
    brca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BRCA.txt");
    hnsc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_HNSC.txt");
    luad.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUAD.txt");
    lusc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUSC.txt");
    expmap["BLCA"] = &blca;
    expmap["BRCA"] = &brca;
    expmap["HNSC"] = &hnsc;
    expmap["LUAD"] = &luad;
    expmap["LUSC"] = &lusc;
    
    vector<CExpressionBin> expBins;
    
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));
    
    
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
        bin = (*expmap[string(mut.cancer)]).GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, expBins, strand, strandInconsistence);
        
        it = results.find(CResultsKey(string(mut.cancer),string(mut.sample),RT_NULLBIN_NOBIN,bin));
        if ( it == results.end())
        {
            rv = CResultsValue(1,0,0,0,0);
            res = results.insert(pair<CResultsKey, CResultsValue>(CResultsKey(string(mut.cancer),string(mut.sample),RT_NULLBIN_NOBIN,bin), rv));
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
    
    f << "Cancer" << '\t' << "Sample" << '\t' << "ExpressionBin" << '\t' << "MutationCnt" << '\t' << "PlusStrandConsistent" << '\t' << "MinusStrandConsistent" << '\t' << "PlusStrandAll" << '\t' << "MinusStrandAll" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.expbin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
    f.close();
}

void CSignatureAnalysis::CalculateTargetsinRTBins(string dir_motifname, CHumanGenome* phuman, int isSignatureMotif)
{
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
    
    set<string> motifs;
    if(isSignatureMotif)
    {
        vector<CMutationSignature>::iterator s;
        for(s=signatures.begin();s!=signatures.end();s++)
            motifs.insert(s->motif);
 
        path = string(RESULTS_FOLDER)+string("/"+dir_motifname+"_in_RTbins_IMR90.txt");
        rtIMR90.CalculateMotifinRTBins(motifs, path.c_str(), phuman);
        path = string(RESULTS_FOLDER)+string("/"+dir_motifname+"_in_RTbins_MCF7.txt");
        rtMCF7.CalculateMotifinRTBins(motifs, path.c_str(), phuman);
        path = string(RESULTS_FOLDER)+string("/"+dir_motifname+"_in_RTbins_NHEK.txt");
        rtNHEK.CalculateMotifinRTBins(motifs, path.c_str(), phuman);
        
    }
    else
    {

        motifs.clear();
        motifs.insert("X");
    
        path = string(RESULTS_FOLDER)+string("/ALL_in_RTbins_IMR90.txt");
        rtIMR90.CalculateMotifinRTBins(motifs, path.c_str(), phuman);
        path = string(RESULTS_FOLDER)+string("/ALL_in_RTbins_MCF7.txt");
        rtMCF7.CalculateMotifinRTBins(motifs, path.c_str(), phuman);
        path = string(RESULTS_FOLDER)+string("/ALL_in_RTbins_NHEK.txt");
        rtNHEK.CalculateMotifinRTBins(motifs, path.c_str(), phuman);
    }

}


set<CCancerSample> CSignatureAnalysis::LoadCancerSamples(string path)
{
    set<CCancerSample> ret;
    string line;
    
    ifstream f(path.c_str());
    if (!f.is_open())
        return(ret);

    vector<string> flds;
    getline(f, line);
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            ret.insert(CCancerSample(flds[0],flds[1]));
        }
    }
    
    f.close();
    
    return(ret);
}

// IMPORTANT: this function contrary to CalculateTargetsinExpBinAllMotifs calculates statistics for complementary motifs
void CSignatureAnalysis::CalculateExpressionAllMotifs(CMutations& muts, string outFilePrefix, CHumanGenome* phuman, string cancer, string sample)
{
    set<string> cancers;
    set<string> samples;
    cancers.insert(cancer);
    samples.insert(sample);
    
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch()    ;
    
    CExpression exp;
    exp.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_"+cancer+".txt");
    
    vector<CExpressionBin> expBins;
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));

    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    CResultsValue rv;
    pair<map<CResultsKey, CResultsValue>::iterator, bool> res;
    
    CMutation mut;
    int strand, strandInconsistence;
    int bin;
    int chrNum;
    string motif;
    int tmp=0;
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];
        if (string(mut.chr) == "M")
            continue;
        if (!cancers.empty() && cancers.find(string(mut.cancer)) == cancers.end())
            continue;
        if (!samples.empty() && samples.find(string(mut.sample)) == samples.end())
            continue;
        
        tmp++;
        bin = exp.GetExpressionBin(string(mut.sample), string(mut.chr), mut.pos, 1, genes, expBins, strand, strandInconsistence);
        chrNum = CHumanGenome::GetChrNum(string(mut.chr));
        motif = string(1,phuman->dna[chrNum][mut.pos-2]) + string(1,phuman->dna[chrNum][mut.pos-1]) + string(1,phuman->dna[chrNum][mut.pos]);
 
        it = results.find(CResultsKey(motif,string(cancer),string(sample),RT_NULLBIN_NOBIN,bin));
        if ( it == results.end())
        {
            rv = CResultsValue(1,0,0,0,0);
            res = results.insert(pair<CResultsKey,CResultsValue>(CResultsKey(motif,string(cancer),string(sample),RT_NULLBIN_NOBIN,bin), rv));
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
        
        motif = CDNA::cDNA(motif);
        it = results.find(CResultsKey(motif,string(cancer),string(sample),RT_NULLBIN_NOBIN,bin));
        if ( it == results.end())
        {
            rv = CResultsValue(1,0,0,0,0);
            res = results.insert(pair<CResultsKey,CResultsValue>(CResultsKey(motif,string(cancer),string(sample),RT_NULLBIN_NOBIN,bin), rv));
            it = res.first;
        }
        else
        {
            it->second.mutCnt++;
        }
        
        if(strand == 0)
        {
            it->second.plusStrandAll++;
            if(strandInconsistence == 0)
                it->second.plusStrandConsistent++;
        }
        else if(strand == 1)
        {
            it->second.minusStrandAll++;
            if(strandInconsistence == 0)
                it->second.minusStrandConsistent++;
        }
    }
    
    cout << "Mutnum: " << tmp << '\n';
    
    ofstream f;
    string path;
    if(cancer == "" && sample == "")
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+".txt";
    else
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+"_"+cancer+"_"+sample+".txt";
    f.open(path.c_str());
    
    f << "Motif" << '\t' << "Cancer" << '\t' << "Sample" << '\t' << "ExpressionBin" << '\t' << "MutationCnt" << '\t' << "PlusStrandConsistent" << '\t' << "MinusStrandConsistent" << '\t' << "PlusStrandAll" << '\t' << "MinusStrandAll" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.motif << '\t' << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.expbin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
    f.close();
    
}

// IMPORTANT: this function contrary to CalculateExpressionAllMotifs DO NOT calculates statistics for complementary motifs
void CSignatureAnalysis::CalculateTargetsinExpBinAllMotifs(string outFilePrefix, CHumanGenome* phuman, string cancer, string sample)
{
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();
    
    CExpression exp;
    exp.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_"+cancer+".txt");

    vector<CExpressionBin> expBins;
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));

    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    CResultsValue rv;
    pair<map<CResultsKey, CResultsValue>::iterator, bool> res;
    
    int i,j;
    int bin;
    int strand;
    int strandInconsistence;
    string motif,cmotif;
    for(i=0;i<phuman->chrCnt;i++)
        for(j=1;j<(phuman->chrLen[i]-1);j++)
        {
            if(phuman->dna[i][j] == 'N' || phuman->dna[i][j-1] == 'N' || phuman->dna[i][j+1] == 'N')
                continue;
            
            bin = exp.GetExpressionBin(sample, phuman->chrName[i], j+1, 1, genes, expBins, strand, strandInconsistence);
            
            motif = string(1,phuman->dna[i][j-1]) + string(1,phuman->dna[i][j]) + string(1,phuman->dna[i][j+1]);

            it = results.find(CResultsKey(motif,string(cancer),string(sample),RT_NULLBIN_NOBIN,bin));
            if ( it == results.end())
            {
                rv = CResultsValue(1,0,0,0,0);
                res = results.insert(pair<CResultsKey,CResultsValue>(CResultsKey(motif,string(cancer),string(sample),RT_NULLBIN_NOBIN,bin), rv));
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
    if(cancer == "" && sample == "")
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+".txt";
    else
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+"_"+cancer+"_"+sample+"_onestrand.txt";
    f.open(path.c_str());
    
    f << "Motif" << '\t' << "Cancer" << '\t' << "Sample" << '\t' << "ExpressionBin" << '\t' << "MutationCnt" << '\t' << "PlusStrandConsistent" << '\t' << "MinusStrandConsistent" << '\t' << "PlusStrandAll" << '\t' << "MinusStrandAll" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.motif << '\t' << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.expbin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
    f.close();
}

void CSignatureAnalysis::CalculateTargetsinExpressionBins(string outFilePrefix, CHumanGenome* phuman, string cancer, string sample, int isSignatureMotif)
{
    int bin;
    int i;
    
    //////// Load data
    
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();

    CExpression blca,brca,hnsc,luad,lusc;
    map<string,CExpression*> expmap;
    
    blca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BLCA.txt");
    brca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BRCA.txt");
    hnsc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_HNSC.txt");
    luad.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUAD.txt");
    lusc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUSC.txt");
    expmap["BLCA"] = &blca;
    expmap["BRCA"] = &brca;
    expmap["HNSC"] = &hnsc;
    expmap["LUAD"] = &luad;
    expmap["LUSC"] = &lusc;
    
    vector<CExpressionBin> expBins;
    
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));

    set<string> motifs;
    if(isSignatureMotif)
    {
        vector<CMutationSignature>::iterator s;
        for(s=signatures.begin();s!=signatures.end();s++)
            motifs.insert(s->motif);
    }
    else
        motifs.insert("X");
    
    ///// Processing
    
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    
    set<CCancerSample>::iterator s;
    set<CCancerSample> cancerSample;
    CResultsValue rv;
    pair<map<CResultsKey, CResultsValue>::iterator, bool> res;
    
    if(cancer == "" && sample == "")
        cancerSample = signatureMuts.cancerSample;
    else
        cancerSample.insert(CCancerSample(cancer,sample));

    for(s=cancerSample.begin();s!=cancerSample.end();s++)
        for(bin=(-EXP_NULLBIN_CNT);bin<(int)expBins.size();bin++)
        {
            rv = CResultsValue(0,0,0,0,0);
            results.insert(pair<CResultsKey, CResultsValue>(CResultsKey((*s).cancer,(*s).sample,RT_NULLBIN_NOBIN,bin), rv));
            cout << "bin:" << bin << "\n";
        }
    
    //set<CCancerSample> doneCancerSamples;
    //doneCancerSamples = LoadCancerSamples(string(RESULTS_FOLDER)+"/TCW_in_EXPbins.txt");
    map<int,unsigned long> binsResults;
    for(s=cancerSample.begin();s!=cancerSample.end();s++)
    {
        cout << "Cancer:" << (*s).cancer << ", Sample:" << (*s).sample << '\n';
        //if(doneCancerSamples.find((*s)) != doneCancerSamples.end())
        //    continue;
        binsResults = (*expmap[(*s).cancer]).CalculateMotifsExpressionBins(expBins, genes, motifs, (*s).sample, phuman);
        for(i=(-EXP_NULLBIN_CNT);i<(int)expBins.size();i++)
            results[CResultsKey((*s).cancer,(*s).sample,RT_NULLBIN_NOBIN,i)].mutCnt = binsResults[i];
    }
    
    
    // Save current results to file
    ofstream f;
    string path;
    if(cancer == "" && sample == "")
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+".txt";
    else
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+"_"+cancer+"_"+sample+".txt";
    f.open(path.c_str());
    
    cout << results.size() << "\n";
    f << "Cancer" << '\t' << "Sample" << '\t' << "ExpressionBin" << '\t' <<  "TargetCnt" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.expbin << '\t' << it->second.mutCnt <<'\n';
        
    f.close();
}



void CSignatureAnalysis::CalculateAPOBECEnrichment(CHumanGenome* phuman_)
{
    CHumanGenome* phuman;
    CMutationSignature s;
    unsigned long cytosineCnt, TCWcnt;
    set<string> motifs;
    ofstream f;
    string path;
    
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
    vector<CMutationSignature>::iterator si;
    for(si=signatures.begin();si!=signatures.end();si++)
        motifs.insert(si->motif);
    TCWcnt = s.CountMotifGenome(motifs, phuman);
    
    cout << "C: " << cytosineCnt << ", TCW: " << TCWcnt << '\n';
    
    if(otherMuts.mutations.size() == 0)
    {
        cerr << "Mutations are not classified." << '\n';
        exit(0);
    }
        
    map<CMapKey, CMapValue> enrichments;
    map<CMapKey, CMapValue>::iterator it;

    for(int i=0;i<signatureMuts.mutations.size();i++)
    {
        it = enrichments.find(CMapKey(signatureMuts.mutations[i].cancer,signatureMuts.mutations[i].sample));
        if(it == enrichments.end())
            enrichments.insert(pair<CMapKey, CMapValue>(CMapKey(signatureMuts.mutations[i].cancer, signatureMuts.mutations[i].sample), CMapValue(1,0,0,0)));
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

void CSignatureAnalysis::AnalyzeRTExpression(CMutations& muts, string resultsFilename)
{
    
    // Load replication timing
    string path;
    CReplicationTiming rtIMR90;
    CReplicationTiming rtMCF7;
    CReplicationTiming rtNHEK;
    int rtBinsSize;
    
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed");
    rtIMR90.LoadReplicationTiming(path.c_str(), 0);
    rtIMR90.ReplicationStrand();
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed");
    rtMCF7.LoadReplicationTiming(path.c_str(), 0);
    rtMCF7.ReplicationStrand();
    path = string(REPLICATION_TIMING_FOLDER)+string("/wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed");
    rtNHEK.LoadReplicationTiming(path.c_str(), 0);
    rtNHEK.ReplicationStrand();
    
    map<string,CReplicationTiming*> rtmap;
    rtmap["BLCA"] = &rtNHEK;
    rtmap["BRCA"] = &rtMCF7;
    rtmap["HNSC"] = &rtNHEK;
    rtmap["LUAD"] = &rtIMR90;
    rtmap["LUSC"] = &rtIMR90;

    vector<CRTBin> binsIMR90, binsNHEK, binsMCF7;
    binsIMR90.push_back(CRTBin(0,-100,13.0766));
    binsIMR90.push_back(CRTBin(1,13.0766,28.3851));
    binsIMR90.push_back(CRTBin(2,28.3851,40.2474));
    binsIMR90.push_back(CRTBin(3,40.2474,52.0254));
    binsIMR90.push_back(CRTBin(4,52.0254,64.7194));
    binsIMR90.push_back(CRTBin(5,64.7194,75.0635));
    binsIMR90.push_back(CRTBin(6,75.0635,90.1735));
    
    binsMCF7.push_back(CRTBin(0,-100,19.8872));
    binsMCF7.push_back(CRTBin(1,19.8872,30.5993));
    binsMCF7.push_back(CRTBin(2,30.5993,40.0364));
    binsMCF7.push_back(CRTBin(3,40.0364,50.556));
    binsMCF7.push_back(CRTBin(4,50.556,60.3904));
    binsMCF7.push_back(CRTBin(5,60.3904,68.9953));
    binsMCF7.push_back(CRTBin(6,68.9953,86.4342));
    
    binsNHEK.push_back(CRTBin(0,-100,22.622));
    binsNHEK.push_back(CRTBin(1,-100,32.1929));
    binsNHEK.push_back(CRTBin(2,-100,41.1202));
    binsNHEK.push_back(CRTBin(3,-100,50.9531));
    binsNHEK.push_back(CRTBin(4,-100,59.6393));
    binsNHEK.push_back(CRTBin(5,-100,67.2079));
    binsNHEK.push_back(CRTBin(6,-100,80.7682));
    
    rtBinsSize = (int)binsIMR90.size();
    
    map<string,vector<CRTBin>*> rtbinsmap;
    rtbinsmap["BLCA"] = &binsNHEK;
    rtbinsmap["BRCA"] = &binsMCF7;
    rtbinsmap["HNSC"] = &binsNHEK;
    rtbinsmap["LUAD"] = &binsIMR90;
    rtbinsmap["LUSC"] = &binsIMR90;

    // Load expression data
    
    CHumanGenes genes;
    genes.LoadGenes(HUMAN_GENES);
    genes.PrepareForSearch();
    
    map<string,CExpression*> expmap;
    CExpression blca,brca,hnsc,luad,lusc;
    
    blca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BLCA.txt");
    brca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BRCA.txt");
    hnsc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_HNSC.txt");
    luad.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUAD.txt");
    lusc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUSC.txt");
    expmap["BLCA"] = &blca;
    expmap["BRCA"] = &brca;
    expmap["HNSC"] = &hnsc;
    expmap["LUAD"] = &luad;
    expmap["LUSC"] = &lusc;
    
    vector<CExpressionBin> expBins;
    
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));

    
    // Processing mutations
    
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    set<CCancerSample>::iterator s;
    int rtbin,expbin;
    CResultsValue rv;
    
    for(s=signatureMuts.cancerSample.begin();s!=signatureMuts.cancerSample.end();s++)
        for(rtbin=(-RT_NULLBIN_CNT);rtbin<rtBinsSize;rtbin++)
            for(expbin=(-EXP_NULLBIN_CNT);expbin<(int)expBins.size();expbin++)
            {
                rv = CResultsValue(0,0,0,0,0,0,0);
                results.insert(pair<CResultsKey, CResultsValue>(CResultsKey((*s).cancer,(*s).sample,rtbin,expbin), rv));
            }
    
    CReplicationTime rt;
    CMutation mut;
    int res;
    int repStrand,transStrand,strandInconsistence;
    
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];

        res = rtmap[string(mut.cancer)]->GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
        if(res)
            rtbin = rtIMR90.GetRTBin(rt.RTvalue, (*rtbinsmap[string(mut.cancer)]));
        else
            rtbin = RT_NULLBIN_NOVALUE;
        
        if((mut.isForwardStrand == 1 && rt.isForward == 1) || (mut.isForwardStrand == 0 && rt.isForward == 0))
            repStrand = STRAND_LAGGING;
        else if((mut.isForwardStrand == 1 && rt.isForward == 0) || (mut.isForwardStrand == 0 && rt.isForward == 1))
            repStrand = STRAND_LEADING;
        else
            repStrand = STRAND_NULL;

    
        expbin = (*expmap[string(mut.cancer)]).GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, expBins, transStrand, strandInconsistence);
        
        it = results.find(CResultsKey(string(mut.cancer),string(mut.sample),rtbin,expbin));
        if (it != results.end())
        {
            it->second.mutCnt++;
            if (repStrand == STRAND_LEADING)
                it->second.leadingCnt++;
            else if (repStrand == STRAND_LAGGING)
                it->second.laggingCnt++;

            if(transStrand == 1)
            {
                it->second.plusStrandAll++;
                if(strandInconsistence == 0)
                    it->second.plusStrandConsistent++;
            }
            else if(transStrand == 0)
            {
                it->second.minusStrandAll++;
                if(strandInconsistence == 0)
                    it->second.minusStrandConsistent++;
            }
            
        }
        else
            cerr << "Element of the result map not found\n";
        
    }
    
    ofstream f;
    path = string(RESULTS_FOLDER)+string("/")+resultsFilename;
    f.open(path.c_str());
    
    f << "Cancer" << '\t' << "Sample" << '\t' << "ReplicationBin" << '\t' << "ExpressionBin" << '\t' << "MutationCnt" << '\t' << "LeadingCnt" << '\t' << "LaggingCnt" << '\t' << "PlusStrandConsistent" << '\t' << "MinusStrandConsistent" << '\t' << "PlusStrandAll" << '\t' << "MinusStrandAll" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.rtbin << '\t' << it->first.expbin << '\t' << it->second.mutCnt << '\t' << it->second.leadingCnt << '\t' << it->second.laggingCnt << '\t' << it->second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll << '\n';
    
    f.close();
    
}

void CSignatureAnalysis::CalculateTargetsinRTexpressionBins(string outFilePrefix, CHumanGenome* phuman, string cancer, string sample, int isSignatureMotif)
{
    int i,j;
    
    // Load replication data
    
    int rtBinsSize;
    
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
    
    vector<CRTBin> binsIMR90, binsNHEK, binsMCF7;
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
    
    rtBinsSize = (int)rtIMR90.bins.size();
    
    map<string,CReplicationTiming*> rtmap;
    rtmap["BLCA"] = &rtNHEK;
    rtmap["BRCA"] = &rtMCF7;
    rtmap["HNSC"] = &rtNHEK;
    rtmap["LUAD"] = &rtIMR90;
    rtmap["LUSC"] = &rtIMR90;
    
    //////// Load expression data
    
    CHumanGenes genes;
    genes.LoadGenes(string(HUMAN_GENES));
    genes.PrepareForSearch();
    
    CExpression blca,brca,hnsc,luad,lusc;
    map<string,CExpression*> expmap;
    
    blca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BLCA.txt");
    brca.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_BRCA.txt");
    hnsc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_HNSC.txt");
    luad.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUAD.txt");
    lusc.LoadExpression(string(EXPRESSION_FOLDER)+"/unpivot_expression_LUSC.txt");
    expmap["BLCA"] = &blca;
    expmap["BRCA"] = &brca;
    expmap["HNSC"] = &hnsc;
    expmap["LUAD"] = &luad;
    expmap["LUSC"] = &lusc;

    vector<CExpressionBin> expBins;
    
    expBins.push_back(CExpressionBin(0,-9999999,0));
    expBins.push_back(CExpressionBin(1,0,25));
    expBins.push_back(CExpressionBin(2,25,100));
    expBins.push_back(CExpressionBin(3,100,300));
    expBins.push_back(CExpressionBin(4,300,550));
    expBins.push_back(CExpressionBin(5,550,1000));
    expBins.push_back(CExpressionBin(6,1000,2000));
    expBins.push_back(CExpressionBin(7,2000,99999999999));
    
    set<string> motifs;
    if(isSignatureMotif)
    {
        vector<CMutationSignature>::iterator s;
        for(s=signatures.begin();s!=signatures.end();s++)
            motifs.insert(s->motif);
    }
    else
        motifs.insert("X");
    
    set<CCancerSample>::iterator s;
    set<CCancerSample> cancerSample;
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    
    if(cancer == "" && sample == "")
        cancerSample = signatureMuts.cancerSample;
    else
        cancerSample.insert(CCancerSample(cancer,sample));

    CRTexpression rtexp;
    map<pair<int,int>,unsigned long> binsResults;
    for(s=cancerSample.begin();s!=cancerSample.end();s++)
    {
        cout << "Cancer:" << (*s).cancer << ", Sample:" << (*s).sample << '\n';

        binsResults = rtexp.CalculateMotifsRTexpressionBins((*rtmap[(*s).cancer]),(*expmap[(*s).cancer]),expBins,genes,motifs,(*s).sample, phuman);
    
        for(i=(-RT_NULLBIN_CNT);i<rtBinsSize;i++)
            for(j=(-EXP_NULLBIN_CNT);j<(int)expBins.size();j++)
                results[CResultsKey((*s).cancer,(*s).sample,i,j)].mutCnt = binsResults[make_pair(i,j)];
    }

    // Save results to file
    ofstream f;
    if(cancer == "" && sample == "")
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+".txt";
    else
        path = string(RESULTS_FOLDER)+"/"+outFilePrefix+"_"+cancer+"_"+sample+".txt";
    f.open(path.c_str());

    f << "Cancer" << '\t' << "Sample" << '\t' << "ReplicationBin" << '\t' << "ExpressionBin" << '\t' << "MutationCnt" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.rtbin << '\t' << it->first.expbin << '\t' << it->second.mutCnt << '\n';
    
    f.close();

}

CRTexpMapKey2::CRTexpMapKey2(string motif_, int RTbin_, int expbin_, int RTstrand_, int senseStrand_)
{
    motif = motif_;
    RTbin = RTbin_;
    expbin = expbin_;
    RTstrand = RTstrand_;
    senseStrand = senseStrand_;
}

CRTexpMapKey::CRTexpMapKey(string motif_, int RTbin_, int expbin_, int RTstrand_, int senseStrand_, char mutbase_)
{
    motif = motif_;
    RTbin = RTbin_;
    expbin = expbin_;
    RTstrand = RTstrand_;
    senseStrand = senseStrand_;
    mutbase = mutbase_;
}

void CSignatureAnalysis::RTExpAllMotifs(CMutations& muts, string dirpath, map<string,CReplicationTiming*> rtmap,     CExpression exp, vector<CExpressionBin> expBins, CHumanGenome* phuman, CHumanGenes genes, string cancer, string sample)
{
    int res;
    int RTbin;
    int expbin;
    int RTstrand;
    int expStrand, strandInconsistence;
    int chrNum;
    int rtBinsSize;
    string motif;
    clock_t t1,t2,t3,t4,t5,st4,st5;
    char nuc[4] = {'A','G','C','T'};
    set<string> motifs;
    set<string>::iterator mit;
    int i,j,k;

    CReplicationTime rt;
    
    map<CRTexpMapKey, unsigned long> results,results2strands;
    map<CRTexpMapKey, unsigned long>::iterator it,it2;
    
    for(i=0;i<4;i++)
        for(j=0;j<4;j++)
            for(k=0;k<4;k++)
            {
                motif = string(1,nuc[i])+string(1,nuc[j])+string(1,nuc[k]);
                motifs.insert(motif);
            }
    
    rtBinsSize = rtmap[cancer]->bins.size();
    
    CMutation mut;
    t1 = clock();
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];
        
        // RT
        res = rtmap[string(mut.cancer)]->GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
        if(res)
            RTbin = rtmap[string(mut.cancer)]->GetRTBin(rt.RTvalue, (rtmap[string(mut.cancer)]->bins));
        else
            RTbin = RT_NULLBIN_NOVALUE;
        
        if((mut.isForwardStrand == 1 && rt.isForward == 0) || (mut.isForwardStrand == 0 && rt.isForward == 1))
            RTstrand = STRAND_LEADING;
        else if((mut.isForwardStrand == 1 && rt.isForward == 1) || (mut.isForwardStrand == 0 && rt.isForward == 0))
            RTstrand = STRAND_LAGGING;
        else
            RTstrand = STRAND_NULL;
        
        // Expression
        expbin = exp.GetExpressionBin(string(mut.sample), string(mut.chr), mut.pos, 1, genes, expBins, expStrand, strandInconsistence);

        // Motif
        chrNum = CHumanGenome::GetChrNum(string(mut.chr));
        motif = string(1,phuman->dna[chrNum][mut.pos-2]) + string(1,phuman->dna[chrNum][mut.pos-1]) + string(1,phuman->dna[chrNum][mut.pos]);
        
        cout << mut.chr << '\t' << mut.pos << '\t' << motif << '\t' << RTbin << '\t' << RTstrand << '\t' << expbin << '\t' << expStrand << '\n';

        // write results to map
        it = results.find(CRTexpMapKey(motif,RTbin,expbin,RTstrand,expStrand,mut.varallele[0]));
        if(it != results.end())
            it->second++;
        else
            results.insert(pair<CRTexpMapKey, unsigned long>(CRTexpMapKey(motif,RTbin,expbin,RTstrand,expStrand,mut.varallele[0]), 1));
    }
  
    /*ofstream f;
    string path;
    path = dirpath + "/RTEXP_MUT_" + sample + "_onestrand.txt";
    f.open(path.c_str());
    for(it=results.begin();it!=results.end();it++)
    f << it->first.motif << '\t' << it->first.RTbin << '\t' << it->first.RTstrand << '\t' << it->first.expbin << '\t' << it->first.senseStrand << '\t' << it->first.mutbase << '\t' << it->second << '\n';
    f.close();*/
    
    results2strands = results;
    
    string cmotif;
    char cvarallele;
    int opRTstrand;
    int opEXPstrand;
    for(it=results.begin();it!=results.end();it++)
    {
        cmotif = CDNA::cDNA(it->first.motif);
        cvarallele = CDNA::cDNA(string(1,it->first.mutbase))[0];
        opRTstrand = CReplicationTiming::oppositeStrand(it->first.RTstrand);
        opEXPstrand = CExpression::oppositeStrand(it->first.senseStrand);
        it2 = results2strands.find(CRTexpMapKey(cmotif,it->first.RTbin,it->first.expbin,opRTstrand,opEXPstrand,cvarallele));
        if(it2 != results2strands.end())
            it2->second += it->second;
        else
            results2strands.insert(pair<CRTexpMapKey, unsigned long>(CRTexpMapKey(cmotif,it->first.RTbin,it->first.expbin,opRTstrand,opEXPstrand,cvarallele), it->second));
    }
    
    ofstream f;
    string path;
    path = dirpath + "/RTEXP_MUT_" + sample + ".txt";
    f.open(path.c_str());
    for(it=results2strands.begin();it!=results2strands.end();it++)
        f << it->first.motif << '\t' << it->first.RTbin << '\t' << it->first.RTstrand << '\t' << it->first.expbin << '\t' << it->first.senseStrand << '\t' << it->first.mutbase << '\t' << it->second << '\n';
    f.close();
    
    // Targets
    
    t2 = clock();
 
    map<CRTexpMapKey2, unsigned long> resultsT,results2strandsT;
    map<CRTexpMapKey2, unsigned long>::iterator itT,it2T;
    
    for(mit=motifs.begin();mit!=motifs.end();mit++)
        for(RTbin=(-RT_NULLBIN_CNT);RTbin<rtBinsSize;RTbin++)
            for(expbin=(-EXP_NULLBIN_CNT);expbin<(int)expBins.size();expbin++)
                for(i=-1;i<2;i++)
                    for(j=-1;j<2;j++)
                        resultsT.insert(pair<CRTexpMapKey2, unsigned long>(CRTexpMapKey2((*mit),RTbin,expbin,i,j), 0));
    
    for(i=0;i<phuman->chrCnt;i++)
    {
        cout << "Chromosome: " << i << '\n';
        st4 = 0; st5 = 0;
       for(j=1;j<(phuman->chrLen[i]-1);j++)
        {
            if(!(CDNA::inACGT(phuman->dna[i][j]) &&
               CDNA::inACGT(phuman->dna[i][j-1]) &&
               CDNA::inACGT(phuman->dna[i][j+1])))
                continue;
            
            // RT
            t4 = clock();
            res = rtmap[cancer]->GetRT(i, j+1, rt);
            if(res)
                RTbin = rtmap[cancer]->GetRTBin(rt.RTvalue, (rtmap[cancer]->bins));
            else
                RTbin = RT_NULLBIN_NOVALUE;
            
            if(rt.isForward == 0)
                RTstrand = STRAND_LEADING;
            else if(rt.isForward == 1)
                RTstrand = STRAND_LAGGING;
            else
                RTstrand = STRAND_NULL;
            
            st4 += (clock()-t4);

            // Expression
            t5 = clock();
            expbin = exp.GetExpressionBin(sample, phuman->chrName[i], j+1, 1, genes, expBins, expStrand, strandInconsistence);
            st5 += (clock()-t5);
            
            
            // Motif
            motif = string(1,phuman->dna[i][j-1]) + string(1,phuman->dna[i][j]) + string(1,phuman->dna[i][j+1]);

            // write results to map
            itT = resultsT.find(CRTexpMapKey2(motif,RTbin,expbin,RTstrand,expStrand));
            if(itT != resultsT.end())
                itT->second++;
            else
                cerr << "Error: map key not found";

        }
        cout << "RT: " << st4 << "\nExpression: " << st5 << '\n';
    }
    
    results2strandsT = resultsT;

    for(itT=resultsT.begin();itT!=resultsT.end();itT++)
    {
        cmotif = CDNA::cDNA(itT->first.motif);
        opRTstrand = CReplicationTiming::oppositeStrand(itT->first.RTstrand);
        opEXPstrand = CExpression::oppositeStrand(itT->first.senseStrand);
        it2T = results2strandsT.find(CRTexpMapKey2(cmotif,itT->first.RTbin,itT->first.expbin,opRTstrand,opEXPstrand));
        if(it2T != results2strandsT.end())
            it2T->second += itT->second;
        else
            cerr << "Error: map key not found";
    }

    path = dirpath + "/RTEXP_TRG_" + sample + ".txt";
    f.open(path.c_str());
    for(itT=results2strandsT.begin();itT!=results2strandsT.end();itT++)
        f << itT->first.motif << '\t' << itT->first.RTbin << '\t' << itT->first.RTstrand << '\t' << itT->first.expbin << '\t' << itT->first.senseStrand << '\t' << itT->second << '\n';
    f.close();
    
    t3 = clock();
    
    cout << "Mutations: " << t2 - t1 << "\n";
    cout << "Targets: " << t3 - t2 << "\n";
    
}
