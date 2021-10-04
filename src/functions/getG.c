#include "math.h"
DAT getG(DAT tw, DAT tg){
    DAT g = 0.0;
    if(tw<=tg){
        g = 9.81*tw*((DAT)1.0/(DAT)tg);
    }
    else{
        g = 9.81;
    }
    return g;
}