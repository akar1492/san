// calc ina
//--new1

#include "sa.h"

REAL
ina(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ena; 
  REAL mss, tm, fna, h, h1ss, th1, h2ss, th2;

  // first time
  for( ; first; first=0 )
    {
      Ena = RTF* log(nao/C->nai);	// here nao and nai =const; Ena defined once
      printf("### Ena=%g\n", Ena);
    }

  
  mss = pow(1./(1+ exp(-S->E/5.46)) , 0.333333333 );

  tm = 6.247e-4/ (0.832*exp(-0.335*(S->E+56.7))+0.627*exp(0.082*(S->E+65.01)))+4e-5;

  Sn->m = mss - (mss-S->m)*exp(-ht/tm);	// solve eq

  //??? we follow article; are there misprints?
  h1ss = h2ss = 1./(1+exp((S->E+66.1)/6.4));

  th1 = 3.717e-6*exp(-0.2815*(S->E+17.1)) / (1+3.732e-3*exp(-0.3426*(S->E+37.76))) + 5.977e-4;

  th2 = 3.186e-8*exp(-0.6219*(S->E+18.8)) / (1+7.189e-5*exp(-0.6683*(S->E+34.07))) + 3.556e-3;
  
  Sn->h1 = h1ss - (h1ss-S->h1)*exp(-ht/th1);	// solve eq
  Sn->h2 = h2ss - (h2ss-S->h2)*exp(-ht/th2);	// solve eq

  fna = 9.52e-2 *exp(-6.3e-2*(S->E+34.4)) / (1.+1.66*exp(-0.225*(S->E+63.7))) + 8.69e-2;

  h = (1-fna)*S->h1 + fna*S->h2;

  // mintau
  //    mintau = MIN( mintau, th1 );
  //  mintau = MIN( mintau, th2 );
  //  mintau = MIN( mintau, tm );


  return C->gna*S->m*S->m*S->m*h*nao*FRD/RTF* (exp((S->E-Ena)/RTF)-1)/(exp(S->E/RTF)-1) *S->E;

} /** ina **/
