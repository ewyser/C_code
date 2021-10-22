#include "math.h"
#include "../include/whichType.h"
int whichType(DAT xn, DAT xmin, DAT xmax, DAT h){
    int type=0;
    if(xn==xmin || xn==xmax){type = 1;}
    else if(xn==xmin+h){type = 2;}
    else if(xn>=xmin+2.0*h && xn<=xmax-2.0*h){type = 3;}
    else if(xn==xmax-h){type = 4;}
    else{type = 0;} 
    return(type);
}