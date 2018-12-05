//
//  g_human.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <string>
#include "ghuman.hpp"

using namespace std;

// Human genome

void CHumanGenome::InitializeHuman(string version_, string pathPrefix, string pathSuffix, string format, int isLoadFiles)
{
    version = version_;
    chrCnt = 24;
    chrName =  new string[chrCnt] {"1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
    dna = new char*[chrCnt];
    
    if (version == "37")
    {
        chrLen = new int[chrCnt] {
            249250621,
            243199373,
            198022430,
            191154276,
            180915260,
            171115067,
            159138663,
            146364022,
            141213431,
            135534747,
            135006516,
            133851895,
            115169878,
            107349540,
            102531392,
            90354753,
            81195210,
            78077248,
            59128983,
            63025520,
            48129895,
            51304566,
            155270560,
            59373566};
    }
    else if (version == "38")
    {
        chrLen = new int[chrCnt] {
            248956422,
            242193529,
            198295559,
            190214555,
            181538259,
            170805979,
            159345973,
            145138636,
            138394717,
            133797422,
            135086622,
            133275309,
            114364328,
            107043718,
            101991189,
            90338345,
            83257441,
            80373285,
            58617616,
            64444167,
            46709983,
            50818468,
            156040895,
            57227415};
    }
    
    if (isLoadFiles)
    {
        Read2Memory(pathPrefix, pathSuffix, format);
    }
}


