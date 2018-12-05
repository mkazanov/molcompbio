//
//  g_human.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef ghuman_hpp
#define ghuman_hpp

#include <stdio.h>
#include "biogenome.hpp"
#include "bioread.hpp"
#include "service.h"

class CHumanGenome : public CGenome
{
public:
    void InitializeHuman(string version, string pathPrefix, string pathSuffix, string format, int isLoadFiles = 1);
    static int GetChrNum(string chr)
    {
        int ret;
        
        if (chr == "X")
            ret = 22;
        else if (chr == "Y")
            ret = 23;
        else
            ret = str2i(chr) - 1;
        
        return(ret);
    }
};



#endif /* ghuman_hpp */
