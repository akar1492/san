// calc T-type iCa
// checked ++
// use tables to speed-up 

#include "sa.h"

////////// tau and ss 
void
tssACH_dt( REAL ee, REAL *tau, REAL *ss )
{
  REAL adt, bdt;

    adt = exp((ee + 26.3f) /30) *1068;
    bdt = exp(-(ee + 26.3f) /30) *1068;
    
    *tau = 1./(adt+bdt);
    *ss  = 1./(exp(-(ee+37)/6.8) +1);
}///tss_dt

////////// tau and ss 
void
tssACH_ft( REAL ee, REAL *tau, REAL *ss )
{
  REAL aft, bft;

    aft = exp(-(ee + 71.7)/83.3) *15.3;
    bft = exp((ee + 71.7) /15.38) *15;

    *tau = 1./(aft + bft);
    *ss  = 1./(exp((ee + 71.)/9.) +1.);

}////tss_ft

///////////////////////////////////////////////////////////////////////
REAL
icatach_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C )
{
  static int first=1; 
  static REAL Ecat; 
  static REAL tdt[NTABLE],dtss[NTABLE], tft[NTABLE],ftss[NTABLE]; 
  REAL tdt_t, dtss_t, tft_t, ftss_t; 

  static REAL gach; 	// factor [0..1] to show deprression Icat(!) by ACH

  // gca should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      int i;

      Ecat = 45;
      printf("### Ecat = %g\n", Ecat);

      //fill tables
      for( i=-1; ++i<NTABLE; )
	{
	  REAL ee;
	  
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);

	  // fill tabs
	  tssACH_dt( ee, tdt+i, dtss+i );
	  tssACH_ft( ee, tft+i, ftss+i );
	}//fori
    }//forfirst

  // outside "first"!!! -RA-
      // ACH depression
      // accord to Zhang's .f code Icat(!) as well as Ical depends on cent/peri
/* #ifdef CENTRAL */
/* 	 gach = 1. - 0.14* ach/(1.2e-7 +ach); //center */
/* #elif defined PERIPHERAL */
/* 	 gach = 1. - 0.56* ach/(1.2e-7 +ach); //peri */
/* #endif */
	 gach = 1. - C->kachicat* ach/(1.2e-7 +ach); 

  // estimate table values
  {
    register int iii=T.n1;
    tdt_t  = tdt[iii] + T.frac*(tdt[iii+1]-tdt[iii]);
    dtss_t = dtss[iii] + T.frac*(dtss[iii+1]-dtss[iii]);
    tft_t  = tft[iii] + T.frac*(tft[iii+1]-tft[iii]);
    ftss_t = ftss[iii] + T.frac*(ftss[iii+1]-ftss[iii]);
  }

  Sn->dt = dtss_t - (dtss_t - S->dt)*exp(-ht/tdt_t);	// solve eq
  Sn->ft = ftss_t - (ftss_t - S->ft)*exp(-ht/tft_t);	// solve eq

  return  gach*C->gcat*S->ft*S->dt*(S->E-Ecat);

} /** icatACH **/
