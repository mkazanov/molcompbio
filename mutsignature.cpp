//
//  mutsignature.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "mutsignature.hpp"

CMutationSignature::CMutationSignature(string motif_, int mutationPos_, string newbase_)
{
    motif = motif_;
    mutationPos = mutationPos_;
    newbase = newbase_;
}

bool CMutationSignature::AnyNewBase()
{
    return((newbase == "X" || newbase == "N"));
}
