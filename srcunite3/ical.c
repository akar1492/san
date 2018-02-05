// calc L-type iCa
// checked; see cooments; ++

#include "sa.h"

REAL
ical(struct State *S, struct State *Sn,REAL ht, struct Cpar *C )
{
  static int first=1; 
  static REAL Ecal; 
  REAL adl, bdl, tdl, dlss, afl, bfl, tfl, flss; 
  // gca should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      Ecal = 46.4; //mV
      printf("### Ecal = %g\n", Ecal);
      //      printf("===%g  %g %g\n", cao, cai, log(cao/cai));
    }

  //!!! adl, bdl, afl checked against exp exceptions
#define ETOL 1e-4 	// tolerance for E


  if( fabs(S->E+35) <ETOL ) 
    adl = 14.19*2.5-42.45*S->E /(exp(-0.208*S->E )-1);
  else if( fabs(S->E) <ETOL )
    adl=-14.19*(S->E +35)/(exp(-(S->E +35)/2.5)-1)+42.45*0.208;
  else	// normal
    adl = -14.19*(S->E+35)/(exp(-(S->E+35)/2.5)-1.) - 42.45*S->E/(exp(-0.208*S->E)-1.); 


  // check for zeros
  if( fabs(S->E-5) <ETOL )
    bdl = 5.71/0.4;
  else		//normal
    bdl = 5.71*(S->E-5.)/(exp(0.4*(S->E-5))-1); 
  
  tdl = 1./(adl+bdl);

  dlss = 1./(1+exp(-(S->E+23.1)/6));

  //  printf("**%g %g %g\n",  S->E, dlss, tdl );

  Sn->dl = dlss - (dlss-S->dl)*exp(-ht/tdl);	// solve eq

  //check for zeros
  if( fabs(S->E +28) <ETOL )
    afl = 3.12*4;
  else		// normal
    afl = 3.12*(S->E+28)/(exp((S->E+28)*0.25) -1);

  bfl = 25./(1+exp(-(S->E+28)*0.25)); //

  tfl = 1./(afl+bfl);

  // check for reality
  if( tfl<ht || tdl<ht ) printf( "!!!ERROR: ht, tdl, tfl = %g %g %g\n", ht, tdl, tfl );

  flss = 1./(1+exp((S->E+45)*0.2));

  Sn->fl = flss - (flss-S->fl)*exp(-ht/tfl);	// solve eq

    // mine
  //    printf("***%g %g %g %g %g %g\n", adl, bdl, afl, bfl, S->fl, S->dl  );

  // mintau
  mintau = MIN(mintau, tdl );
  mintau = MIN(mintau, tfl );


  return C->gcal*( S->fl*S->dl + 0.006/(1+exp(-(S->E+14.1)/6)) )*(S->E-Ecal);
} /** ical **/
