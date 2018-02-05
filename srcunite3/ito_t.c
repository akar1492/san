// calc Ito - transient
// checked see notes for tq!!!
// use tables to speed-up 
#include "sa.h"

// use nonconst nai, ki


////////// tau and ss 
void
tss_q( REAL ee, REAL *tau, REAL *ss )
{
  *ss = 1./ (exp((ee + 59.37)/13.1) + 1);
  // as in the paper --- likely an ERROR here 
  //  tq = 10.1e-3 + 65.17e-3/( 0.57*exp(-0.08*(ee+49))) + 0.24e-4*exp(0.1*(ee+50.93));

  // in the .f file
  *tau = 10.1e-3 + 65.17e-3/( 0.5686*exp(-0.08161*(ee+49.0))+0.7174*exp(0.2719*(ee+50.93)) );

}///tss_q

#if 0
// this func was defined in isus
////////// tau and ss 
void
tss_r( REAL ee, REAL *tau, REAL *ss )
{
  *ss = 1./(exp(-(ee-10.93)/19.7) +1);
  *tau = 15.59e-3/( 1.037*exp((ee+30.61)*0.09) + 0.369*exp(-0.12*(ee+23.84)) ) +2.98e-3;

}///tss_r
#endif

///////////////////////////////////////////////////////////////////////
REAL
ito_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C )
{
  static int first=1; 
  static REAL Ek; 
  static REAL tq[NTABLE],qss[NTABLE];
  REAL tq_t, qss_t;
  static REAL tr[NTABLE],rss[NTABLE];
  REAL tr_t, rss_t;

  // gto should be here, but it proved to be varied from center out

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
	  tss_q( ee, tq+i, qss+i );
	  tss_r( ee, tr+i, rss+i );
	}//fori
      }//fill
    }//forfirst

  // estimate table values
  {
    register int iii=T.n1;
    tq_t  = tq[iii] + T.frac*(tq[iii+1]-tq[iii]);
    qss_t = qss[iii] + T.frac*(qss[iii+1]-qss[iii]);
    tr_t  = tr[iii] + T.frac*(tr[iii+1]-tr[iii]);
    rss_t = rss[iii] + T.frac*(rss[iii+1]-rss[iii]);
  }

  Sn->q = qss_t - (qss_t - S->q)*exp(-ht/tq_t);	// solve eq



  // avoid solving it twice (in isus)
  //   Sn->r = rss_t - (rss_t - S->r)*exp(-ht/tr_t);	// solve eq

  return C->gto*S->q*S->r*(S->E-C->EK);

} /** ito_t **/
