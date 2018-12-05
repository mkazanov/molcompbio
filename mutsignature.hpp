//
//  mutsignature.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef mutsignature_hpp
#define mutsignature_hpp

#include <string>
#include <stdio.h>

using namespace std;

class CMutationSignature {
public:
    string motif;
    int mutationPos;
    string newbase; //'X' or 'N' - any base
    CMutationSignature(string motif_, int mutationPos_, string newbase_);
    CMutationSignature(){};
    bool AnyNewBase();
};

#endif /* mutsignature_hpp */
