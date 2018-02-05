 // Ip = Inak 
// checked ++

// use nonconst nai, ki

#include "sa.h"

REAL
ip(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1;
  static REAL b1, b2, kmna, kmk, coeff;

 // first time
  for( ; first; first=0 )
    {
      kmna = 5.64;
      kmk = 0.621;

      // NOTE: b1, b2 may be calulated just once  
      //      b1 = nai/(kmna+nai);
      //      b2 = ko/(kmk+ko);
      //      coeff = ipss*b1*b1*b1*b2*b2; 
    }//forfirst

      b1 = C->nai/(kmna+C->nai);
      b2 = ko/(kmk+ko);

      coeff = C->ipss*b1*b1*b1*b2*b2; 

  return coeff* 1.6/(1.5+exp(-(S->E+60)/40));

} /** ip **/
