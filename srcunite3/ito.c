// calc Ito - transient
// checked see notes for tq!!!

#include "sa.h"

REAL
ito(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ek; 
  REAL tq, qss, tr, rss;
  // gto should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      Ek = RTF* log(ko/C->ki);	// here ko and ki =const; Ek defined once
      printf("### Ek=%g\n", Ek);
    }

  qss = 1./ (exp((S->E + 59.37)/13.1) + 1);

  // as in the paper --- likely an ERROR here 
  //  tq = 10.1e-3 + 65.17e-3/( 0.57*exp(-0.08*(S->E+49))) + 0.24e-4*exp(0.1*(S->E+50.93));

  // in the .f file
  tq = 10.1e-3 + 65.17e-3/( 0.5686*exp(-0.08161*(S->E+49.0))+0.7174*exp(0.2719*(S->E+50.93)) );

  
 
  Sn->q = qss - (qss-S->q)*exp(-ht/tq);		// solve eq


  rss = 1./(exp(-(S->E-10.93)/19.7) +1);

  tr = 15.59e-3/( 1.037*exp((S->E+30.61)*0.09) + 0.369*exp(-0.12*(S->E+23.84)) ) +2.98e-3;


  // avoid solving it twice (in isus)
  //  Sn->r = rss - (rss-S->r)*exp(-ht/tr);		//solve eq

  return C->gto*S->q*S->r*(S->E-Ek);

} /** ito **/
