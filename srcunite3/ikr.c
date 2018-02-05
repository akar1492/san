// calc i K r, rapid delayed rectifying
// checked ++


#include "sa.h"

REAL
ikr(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ek; 
  REAL pafss, tpaf, passs, tpas, piiss, tpii, fkr, pa;
  // gkr should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      Ek = RTF* log(ko/C->ki);	// here nao and ki =const; Ek defined once
      printf("### Ek = %g\n", Ek);
    }

    pafss = passs = 1./(exp(-(S->E+14.2)/10.6) + 1);

    tpaf = 1./( exp((S->E-9)/15.9) *37.2 + exp(-(S->E-9)/22.5) *0.96 );

    tpas = 1./( exp((S->E-9)/17) *4.2 + exp(-(S->E-9)/21.6) *0.15);

    Sn->paf = pafss - (pafss-S->paf)*exp(-ht/tpaf);	// solve eq
    Sn->pas = passs - (passs-S->pas)*exp(-ht/tpas);	// solve eq


    piiss = 1./(exp((S->E+18.6)/10.1) + 1);
    tpii = 0.002;

    Sn->pii = piiss - (piiss-S->pii)*exp(-ht/tpii);	// solve eq

    fkr = 0.4;
    pa = (1-fkr)*S->paf + fkr*S->pas;

 
    return C->gkr*pa*S->pii*(S->E-Ek);

} /** ikr **/
