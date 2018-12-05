//
//  dna.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef dna_hpp
#define dna_hpp

#include <string>
#include <algorithm>

using namespace std;

class CDNA{
public:
    static string cDNA(string dna)
    {
        string ret;
        
        ret = dna;
        transform(ret.begin(),ret.end(),ret.begin(),::toupper);
        reverse(ret.begin(),ret.end());
        replace(ret.begin(),ret.end(),'A','X');
        replace(ret.begin(),ret.end(),'C','Y');
        replace(ret.begin(),ret.end(),'G','Z');
        replace(ret.begin(),ret.end(),'T','A');
        replace(ret.begin(),ret.end(),'X','T');
        replace(ret.begin(),ret.end(),'Y','G');
        replace(ret.begin(),ret.end(),'Z','C');
        
        return(ret);
    }
};

#endif /* dna_hpp */
