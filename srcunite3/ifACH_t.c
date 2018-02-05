// calc if, Hyperpolarization-activated current
// checked ++
// use tables to speed-up
// Shift activation curve due to ACH

// use nonconst nai, ki

#include "sa.h"
#include <math.h>

////////// tau and ss 
void
tssACH_y( REAL ee, REAL *tau, REAL *ss )
{
  REAL ay, by;
  REAL achshft;

  // ACH consts for If
  //#define smax 	(-7.2) 	// /mv/
#define smax 	(-7.5) 	// /mv/  (Zhang's .f)
#define nf	0.69	// exponent
#define k05	1.26e-8	// /M/ ACh
  
  // Note: to speed-up m.b. 
  // 1) make achshft static or define(if ach=const)
  // 2) assign buffer for pow(ach,nf)
  achshft = smax* pow(ach,nf) /( pow(k05,nf) + pow(ach,nf) );


  ay = exp(-(ee + 78.91 -achshft)/26.62);
  by = exp( (ee + 75.13 -achshft)/21.25);

  *tau = 1./ (ay + by);
  *ss = *tau * ay;
}///tss_y

///////////////////////////////////////////////////////////////////////
REAL
ifach_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C , struct Is *I)
{
  static int first=1; 
  static REAL Ena, Ek;
  static REAL ty[NTABLE],yss[NTABLE]; 
  REAL ty_t, yss_t, ifna, ifk;
  // g** should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    {
      //      Ena = RTF*log(nao/nai); // here nao and nai =const; Ena defined once
      //      Ek = RTF*log(ko/ki); // here ko and ki =const; Ek defined once

      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tssACH_y( ee, ty+i, yss+i );
	}//fori
      }//filltables
    }//forfirst
  ///////////////////////////////////////////////////


  // estimate table values
/*   { */
/*     register int iii=T.n1; */
/*     ty_t  = ty[iii] + T.frac*(ty[iii+1]-ty[iii]); */
/*     yss_t = yss[iii] + T.frac*(yss[iii+1]-yss[iii]); */
/*   } */
 
  // outside "first"!!! -RA-
  //!!!! do not use table; calc directly!!!-- to account nonconst ACh
  tssACH_y( S->E, &ty_t, &yss_t );
    
    Sn->y = yss_t - (yss_t - S->y)*exp(-ht/ty_t);	// solve eq

    /// fifna, fifk -- extern vars to keep track of these currents separately
    I->fifna = ifna = C->gfna*S->y*(S->E - C->ENA);
    I->fifk  = ifk  = C->gfk*S->y*(S->E - C->EK);

    return ifna+ifk;

} /** ifach_t **/
