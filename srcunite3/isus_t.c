// calc Isus(=Ikp) -- very much like Ito; for speed may be joined to ito
// checked ++
// use tables to speed-up 
#include "sa.h"

// use nonconst nai, ki

////////// tau and ss 
void
tss_r( REAL ee, REAL *tau, REAL *ss )
{
  *ss  = 1./(exp(-(ee-10.93)/19.7) +1);
  *tau = 15.59e-3/( 1.037*exp((ee+30.61)*0.09) + 0.369*exp(-0.12*(ee+23.84)) ) +2.98e-3;
}///tss_r

///////////////////////////////////////////////////////////////////////
REAL
isus_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C )
{
  static int first=1; 
  static REAL Ek; 
  static REAL tr[NTABLE],rss[NTABLE];
  REAL tr_t, rss_t;
  // gsus should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      //      Ek = RTF* log(ko/ki);	// here ko and ki =const; Ek defined once
      printf("### Ek=%g\n", Ek);
      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_r( ee, tr+i, rss+i );
	}//fori
      }//fill
    }//forfirst

   // estimate table values
  {
    register int iii=T.n1;
    tr_t  = tr[iii] + T.frac*(tr[iii+1]-tr[iii]);
    rss_t = rss[iii] + T.frac*(rss[iii+1]-rss[iii]);
  }
 
    Sn->r = rss_t - (rss_t - S->r)*exp(-ht/tr_t);	// solve eq


  return C->gsus*S->r*(S->E-C->EK);

} /** isus_t **/
