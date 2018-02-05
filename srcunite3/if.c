// calc if, Hyperpolarization-activated current
// checked ++

#include "sa.h"

REAL
if_(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ena, Ek; 
  REAL ay, by, ty, yss, ifna, ifk;
  // g** should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      // if ions!=const, Ena, Ek should be calculated in the prog
      Ena = RTF*log(nao/C->nai); // here nao and nai =const; Ena defined once
      Ek = RTF*log(ko/C->ki); // here ko and ki =const; Ek defined once
    }

    ay = exp(-(S->E + 78.91)/26.62);
    by = exp((S->E + 75.13)/21.25);

    ty = 1./ (ay + by);
    yss = ty * ay;
    
    Sn->y = yss - (yss-S->y)*exp(-ht/ty);	// solve eq

    ifna = C->gfna*S->y*(S->E-Ena);
    ifk = C->gfk*S->y*(S->E-Ek);

    return ifna+ifk;

} /** if_ **/
