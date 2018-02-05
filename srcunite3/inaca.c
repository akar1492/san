// calc Inaca
// checked ++

// use nonconst nai, ki

#include "sa.h"

REAL
inaca(struct State *S, struct State *Sn, REAL ht, struct Cpar *C, struct Caintra_state *Ca )
{
  static int first=1; 
  static REAL  gnaca, dnaca; 
  static REAL nai3, nao3;
  // g** should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      // consts
      gnaca = 0.5;
      dnaca = 1e-4;

      // NOTE: to speed-up, nai3, nao3 may be calculated just once  
      //      nai3 = nai*nai*nai;
      //      nao3 = nao*nao*nao;
    }//forfirst

  //      printf("###### nai, nai3 = %g %g\n", nai, nai3 );

  nai3 = C->nai*C->nai*C->nai;
  nao3 = nao*nao*nao;

  return C->knaca*( nai3*cao  * exp(S->E*0.03743*gnaca) - nao3*Ca->cai * exp(S->E*0.03743f*(gnaca-1)))  / (1+ dnaca*(nao3*Ca->cai + nai3*cao));


} /** inaca **/
