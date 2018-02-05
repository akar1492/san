// calc T-type iCa
// checked ++

#include "sa.h"

REAL
icat(struct State *S, struct State *Sn, REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ecat; 
  REAL adt, bdt, tdt, dtss, aft, bft, tft, ftss; 
  // gca should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      Ecat = 45;
      printf("### Ecat = %g\n", Ecat);
    }

    adt = exp((S->E + 26.3f) /30) *1068;
    bdt = exp(-(S->E + 26.3f) /30) *1068;
    tdt = 1./(adt+bdt);
    dtss =  1./(exp(-(S->E+37)/6.8) +1);

    Sn->dt = dtss - (dtss-S->dt)*exp(-ht/tdt);	// solve eq

    aft = exp(-(S->E + 71.7)/83.3) *15.3;
    bft = exp((S->E + 71.7) /15.38) *15;
    tft = 1./(aft + bft);
    ftss = 1./(exp((S->E + 71.)/9.) +1.);

    Sn->ft = ftss - (ftss-S->ft)*exp(-ht/tft);	// solve eq

    return C->gcat*S->ft*S->dt*(S->E-Ecat);

} /** icat **/
