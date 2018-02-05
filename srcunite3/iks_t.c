// calc iKs, slow delayed rectifying
// checked ++
// use tables to speed-up 

// use nonconst nai, ki

#include "sa.h"

////////// tau and ss 
void
tss_xs( REAL ee, REAL *tau, REAL *ss )
{
  REAL axs, bxs;

  axs = 14./(exp(-(ee-40)/9.) + 1);
  bxs = exp(-(ee)/45.);

  *tau = 1./(axs+bxs);
  *ss  = axs* *tau;
}///tss_xs

///////////////////////////////////////////////////////////////////////
REAL
iks_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C )
{
  static int first=1; 
  static REAL Eks; 
  static REAL txs[NTABLE],xsss[NTABLE];
  REAL txs_t, xsss_t;
  float pnak = 0.12;


  // gkr should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    { 
      //      Eks = RTF*log((ko+pnak*nao)/(ki + pnak*nai)); // here nao and ki =const; Eks defined once
      // printf("### Eks=%g\n", Eks);

      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_xs( ee, txs+i, xsss+i );
	}//fori
      }//fill
    }//forfirst

  // estimate table values
  {
    register int iii=T.n1;
    txs_t  = txs[iii] + T.frac*(txs[iii+1]-txs[iii]);
    xsss_t = xsss[iii] + T.frac*(xsss[iii+1]-xsss[iii]);
  }
 
    Sn->xs = xsss_t - (xsss_t - S->xs)*exp(-ht/txs_t);	// solve eq

    // asume nai, ki -- nonconst 
      Eks = RTF*log((ko+pnak*nao)/(C->ki + pnak*C->nai)); 

    return C->gks*S->xs*S->xs*(S->E-Eks);

} /** iks_t **/
