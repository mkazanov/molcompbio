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
    
    map<string,CExpression*> expmap;
    CExpression blca,brca,hnsc,luad,lusc;
    
    blca.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_BLCA.txt");
    brca.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_BRCA.txt");
    hnsc.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_HNSC.txt");
    luad.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_LUAD.txt");
    lusc.LoadExpression("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/expression/unpivot_expression_LUSC.txt");
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
        bin = GetExpressionBin(string(mut.sample), mut.chr, mut.pos, mut.isForwardStrand, genes, (*expmap[string(mut.cancer)]), expBins, strand, strandInconsistence);
        
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
    
    f << "Cancer" << '\t' << "Sample" << '\t' << "PlusStrandConsistent" << '\t' << "MinusStrandConsistent" << '\t' << "PlusStrandAll" << '\t' << "MinusStrandAll" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.bin << '\t' << it->second.mutCnt << '\t' << it-> second.plusStrandConsistent << '\t' << it->second.minusStrandConsistent << '\t' << it->second.plusStrandAll << '\t' << it->second.minusStrandAll <<'\n';
    
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


void CAPOBEC::CalculateTargetsinExpressionBins(CHumanGenome* phuman, string cancer, string sample)
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
    motifs.insert("TCA");
    motifs.insert("TCT");
    set<string> motifsall;
    set<string>::iterator si;
    int motiflen;
    
    CMutationSignature msobj;
    char** motifsarr;
    msobj.CheckMotifsNotEmpty(motifs);
    msobj.CheckMotifsSameLength(motifs);
    motifsall = msobj.AddcMotifs(motifs);
    motiflen = (motifsall.begin())->length();
    
    // Prepare char array for copying motifs
    motifsarr = new char*[motifsall.size()];
    for(i=0;i<motifsall.size();i++)
        motifsarr[i] = new char[motiflen];
    
    // Copy motifs to char array
    i = 0;
    for(si=motifsall.begin();si!=motifsall.end();si++)
    {
        strncpy(motifsarr[i],(*si).c_str(),(*si).length());
        i++;
    }
    
    ///// Processing
    
    map<CResultsKey, CResultsValue> results;
    map<CResultsKey, CResultsValue>::iterator it;
    
    CDNAPos pos = CDNAPos(0,0);
    int includeCurPos=1;
    set<CCancerSample>::iterator s;
    set<CCancerSample> cancerSample;
    int strand;
    int strandInconsistence;
    CResultsValue rv;
    pair<map<CResultsKey, CResultsValue>::iterator, bool> res;
    
    clock_t time,time1=0,time2=0, time3=0;
    
    if(cancer == "" && sample == "")
        cancerSample = apobecMuts.cancerSample;
    else
        cancerSample.insert(CCancerSample(cancer,sample));
    
    for(s=cancerSample.begin();s!=cancerSample.end();s++)
        for(bin=-2;bin<8;bin++)
        {
            rv = CResultsValue(0,0,0,0,0);
            results.insert(pair<CResultsKey, CResultsValue>(CResultsKey((*s).cancer,(*s).sample,bin), rv));

        }
    
    int motifsnum;
    motifsnum = motifsall.size();
    unsigned long cnt[10];
    
    //set<CCancerSample> doneCancerSamples;
    //doneCancerSamples = LoadCancerSamples(string(RESULTS_FOLDER)+"/TCW_in_EXPbins.txt");
    for(s=cancerSample.begin();s!=cancerSample.end();s++)
    {
        cout << "Cancer:" << (*s).cancer << ", Sample:" << (*s).sample << '\n';
        //if(doneCancerSamples.find((*s)) != doneCancerSamples.end())
        //    continue;
        for(i=0;i<10;i++)
            cnt[i] = 0;
        time = clock();
        for(pos=msobj.NextMotif(CDNAPos(0,0),motifsarr,motifsnum,motiflen,phuman,END_GENOME,includeCurPos);
            !pos.isNull();
            pos=msobj.NextMotif(pos,motifsarr,motifsnum,motiflen,phuman,END_GENOME))
        {
            time3 += (clock() - time);
            
            time = clock();
            bin = GetExpressionBin((*s).sample, phuman->chrName[pos.chrNum], pos.pos, 1, genes, (*expmap[(*s).cancer]), expBins, strand, strandInconsistence);
            time1 += (clock()-time);
            
            time = clock();
            cnt[bin+2]++;
            time2 += (clock()-time);
            
            time = clock();
        }
        for(i=0;i<10;i++)
            results[CResultsKey((*s).cancer,(*s).sample,i-2)].mutCnt = cnt[i];
        
        // Save current results to file
        ofstream f;
        string path;
        if(cancer == "" && sample == "")
            path = string(RESULTS_FOLDER)+"/TCW_in_EXPbins.txt";
        else
            path = string(RESULTS_FOLDER)+"/TCW_in_EXPbins_"+cancer+"_"+sample+".txt";
        f.open(path.c_str());
        
        f << "Cancer" << '\t' << "Sample" << '\t' << "ExpressionBin" << '\t' <<  "TargetCnt" << '\n';
        for(it=results.begin(); it!=results.end(); ++it)
            f << it->first.cancer << '\t' << it->first.sample << '\t' << it->first.bin << '\t' << it->second.mutCnt <<'\n';
        
        f.close();
        
        cout << "Time1:" << time1 << ", time2:" << time2 << ", time3: " << time3 << '\n';
    }
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

