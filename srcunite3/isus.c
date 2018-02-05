// calc Isus(=Ikp) -- very much like Ito; for speed may be joined to ito
// checked ++

#include "sa.h"

REAL
isus(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ek; 
  REAL tr, rss;
  // gsus should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      Ek = RTF* log(ko/C->ki);	// here ko and ki =const; Ek defined once
      printf("### Ek=%g\n", Ek);
    }

  rss = 1./(exp(-(S->E-10.93)/19.7) +1);

  tr = 15.59e-3/( 1.037*exp((S->E+30.61)*0.09) + 0.369*exp(-0.12*(S->E+23.84)) ) +2.98e-3;


  Sn->r = rss - (rss-S->r)*exp(-ht/tr);		//solve eq

  return C->gsus*S->r*(S->E-Ek);

} /** isus **/
