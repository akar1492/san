// calc iKs, slow delayed rectifying
// checked ++


#include "sa.h"

REAL
iks(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Eks; 
  REAL axs, bxs, txs, xss;
  // gkr should be here, but it proved to be varied from center out
  float pnak = 0.12;

  // first time
  for( ; first; first=0 )
    {
      //      Eks = RTF*log((ko+pnak*nao)/(ki + pnak*nai)); // here nao and ki =const; Eks defined once
      // printf("### Eks=%g\n", Eks);
    }

    axs = 14./(exp(-(S->E-40)/9.) + 1);
    bxs = exp(-(S->E)/45.);

    txs = 1./(axs+bxs);
    xss = axs*txs;
     //RA
    //    printf("+++%g %g %g %g\n", xss, txs, axs, bxs );
   
    Sn->xs = xss - (xss-S->xs)*exp(-ht/txs);	// solve eq

    // asume nai, ki -- nonconst 
    Eks = RTF*log((ko+pnak*nao)/(C->ki + pnak*C->nai)); 

    return C->gks*S->xs*S->xs*(S->E-Eks);

} /** iks **/
