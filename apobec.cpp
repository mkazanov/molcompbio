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
#include <ctime>
#include <cstring>
#include "RTexpression.hpp"

CResultsKey::CResultsKey(string cancer_, string sample_, int rtbin_, int expbin_)
{
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
    signatures.push_back(CMutationSignature("TCA",2,"T"));
    signatures.push_back(CMutationSignature("TCT",2,"T"));
    signatures.push_back(CMutationSignature("TCA",2,"G"));
    signatures.push_back(CMutationSignature("TCT",2,"G"));
    
    set<string> cancers;
    set<string> samples;
    cancers.insert("BLCA");
    cancers.insert("BRCA");
    cancers.insert("HNSC");
    cancers.insert("LUAD");
    cancers.insert("LUSC");
    m.FilterMutations(apobecMuts,signatures,(*phuman),cancers,samples,&otherMuts);
    //apobecMuts.SaveToFile("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations_apobec.tsv");
    
    apobecMuts.GetUniqueCancersSamples();
    
    // Free human genome
    if(!phuman_)
    {
        for(int i=0;i<phuman->chrCnt;i++)
            delete phuman->dna[i];
        delete phuman->dna;
    }
}

void CAPOBEC::AnalyzeReplicationTiming(CMutations& muts, string resultsFilename)
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
    
    map<string,vector<CRTBin>*> rtbinsmap;
    rtbinsmap["BLCA"] = &binsNHEK;
    rtbinsmap["BRCA"] = &binsMCF7;
    rtbinsmap["HNSC"] = &binsNHEK;
    rtbinsmap["LUAD"] = &binsIMR90;
    rtbinsmap["LUSC"] = &binsIMR90;
    
    
    for(int i=0;i<muts.mutations.size();i++)
    {
        mut = muts.mutations[i];
                
        res = rtmap[string(mut.cancer)]->GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
        if(res)
            bin = rtmap[string(mut.cancer)]->GetRTBin(rt.RTvalue, (*rtbinsmap[string(mut.cancer)]));
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

void CAPOBEC::AnalyzeExpression(CMutations& muts, string resultsFilename)
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
    
    f << "Cancer" << '\t' << "Sample" << '\t' << "MutationCnt" << '\t' << "PlusStrandConsistent" << '\t' << "MinusStrandConsistent" << '\t' << "PlusStrandAll" << '\t' << "MinusStrandAll" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.expbin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
    f.close();
}

void CAPOBEC::CalculateTargetsinRTBins(CHumanGenome* phuman)
{
    vector<CRTBin> binsIMR90, binsNHEK, binsMCF7;
    binsIMR90.push_back(CRTBin(1,-100,13.0766));
    binsIMR90.push_back(CRTBin(2,13.0766,28.3851));
    binsIMR90.push_back(CRTBin(3,28.3851,40.2474));
    binsIMR90.push_back(CRTBin(4,40.2474,52.0254));
    binsIMR90.push_back(CRTBin(5,52.0254,64.7194));
    binsIMR90.push_back(CRTBin(6,64.7194,75.0635));
    binsIMR90.push_back(CRTBin(7,75.0635,90.1735));
    
    binsMCF7.push_back(CRTBin(1,-100,19.8872));
    binsMCF7.push_back(CRTBin(2,19.8872,30.5993));
    binsMCF7.push_back(CRTBin(3,30.5993,40.0364));
    binsMCF7.push_back(CRTBin(4,40.0364,50.556));
    binsMCF7.push_back(CRTBin(5,50.556,60.3904));
    binsMCF7.push_back(CRTBin(6,60.3904,68.9953));
    binsMCF7.push_back(CRTBin(7,68.9953,86.4342));
    
    binsNHEK.push_back(CRTBin(1,-100,22.622));
    binsNHEK.push_back(CRTBin(2,-100,32.1929));
    binsNHEK.push_back(CRTBin(3,-100,41.1202));
    binsNHEK.push_back(CRTBin(4,-100,50.9531));
    binsNHEK.push_back(CRTBin(5,-100,59.6393));
    binsNHEK.push_back(CRTBin(6,-100,67.2079));
    binsNHEK.push_back(CRTBin(7,-100,80.7682));
    
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


set<CCancerSample> CAPOBEC::LoadCancerSamples(string path)
{
    set<CCancerSample> ret;
    string line;
    int chrNum;
    
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


void CAPOBEC::CalculateTargetsinExpressionBins(string outFilePrefix, CHumanGenome* phuman, string cancer, string sample, int isAPOBECmotif)
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
    if(isAPOBECmotif)
    {
        motifs.insert("TCA");
        motifs.insert("TCT");
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
        cancerSample = apobecMuts.cancerSample;
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



void CAPOBEC::CalculateAPOBECEnrichment(CHumanGenome* phuman_)
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

void CAPOBEC::AnalyzeRTExpression(CMutations& muts, string resultsFilename)
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
    
    for(s=apobecMuts.cancerSample.begin();s!=apobecMuts.cancerSample.end();s++)
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

void CAPOBEC::CalculateTargetsinRTexpressionBins(string outFilePrefix, CHumanGenome* phuman, string cancer, string sample, int isAPOBECmotif)
{
    int i,j;
    
    // Load replication data
    
    int rtBinsSize;
    
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
    if(isAPOBECmotif)
    {
        motifs.insert("TCA");
        motifs.insert("TCT");
    }
    else
        motifs.insert("X");
    
    set<CCancerSample>::iterator s;
    set<CCancerSample> cancerSample;
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    
    if(cancer == "" && sample == "")
        cancerSample = apobecMuts.cancerSample;
    else
        cancerSample.insert(CCancerSample(cancer,sample));

    CRTexpression rtexp;
    map<pair<int,int>,unsigned long> binsResults;
    for(s=cancerSample.begin();s!=cancerSample.end();s++)
    {
        cout << "Cancer:" << (*s).cancer << ", Sample:" << (*s).sample << '\n';

        binsResults = rtexp.CalculateMotifsRTexpressionBins((*rtmap[(*s).cancer]),(*expmap[(*s).cancer]),(*rtbinsmap[string((*s).cancer)]),expBins,genes,motifs,(*s).sample, phuman);
    
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
