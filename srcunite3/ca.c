/* Intracellular Ca dynamics as in: 
   [Dynamical description of sinoatrial node pacemaking: 
    improved mathematical model for primary pacemaker cell.
    Yasutaka Kurata, Ichiro Hisatome, Sunao Imanishi, and Toshishige Shibamoto
    Am J Physiol Heart Circ Physiol 283: H2074­H2101, 2002;]
*/
// Rubin, September 2003
#include <stdio.h>

#include "sa.h"
#include "ca.h"

/* Return: Cai
   args: Cai_old, ht, ical, icat, inaca
*/


REAL ca_intra(REAL ht, REAL icatotal, struct Cpar *C, struct Caintra_state *Ca )
{
  REAL caret_new; 				// return value
  REAL dftc, dftmc, dftmm, dfcmi, dfcms, dfcq;	// time derivatives
  REAL jcadiff, jrel, jup, jtr; 		// fluxes
  static int first=1;


  ht *= 1000; 	// time in ms (was s)

  //NOTE! we use cai from the argument
  //  Ca.casub = cai;	
  Ca->cai = C->cai;	

  // Ca++ buffering -- derivatives on t
  dftc 	= kftc*Ca->cai*(1.-Ca->ftc)-kbtc*Ca->ftc;
  dfcmi = kfcm*Ca->cai*(1.-Ca->fcmi)-kbcm*Ca->fcmi;
  dfcms = kfcm*Ca->casub*(1.-Ca->fcms)-kbcm*Ca->fcms;
  dfcq 	= kfcq*Ca->carel*(1.-Ca->fcq)-kbcq*Ca->fcq;
  dftmc = kftmc*Ca->cai*(1.-Ca->ftmc-Ca->ftmm)-kbtmc*Ca->ftmc;
  dftmm = kftmm*mgi*(1.-Ca->ftmc-Ca->ftmm)-kbtmm*Ca->ftmm; 

///// TEMP!!!!test
/*   dftc 	=0;  */
/*   dfcmi =0; */
/*   dfcms =0; */
/*   dfcq 	=0; */
/*   dftmc =0; */
/*   dftmm =0; */

  // Ca++ fluxes: diffusion and because of SR
  jcadiff = (Ca->casub-Ca->cai)/tdiffca;
  jrel 	= prel*(Ca->carel-Ca->casub)/(1.+(krel*krel/(Ca->casub*Ca->casub)));
  jup	= pup/(1.+kup/Ca->cai);
  jtr	= (Ca->caup-Ca->carel)/ttr;

  // Now, Ca concentrations
  // NOTE! we substitute new values directly to the struct!

  Ca->cai = Ca->cai + ht*( (jcadiff*C->vsub-jup*C->vup)/C->vi - (cmtot*dfcmi + tctot*dftc + tmctot*dftmc) );

  Ca->casub	= Ca->casub + ht*( (   -1000*(icatotal)/(2.*FRD)  + jrel*C->vrel)/C->vsub -jcadiff - cmtot*dfcms );

  Ca->carel 	= Ca->carel + ht*( jtr - jrel  -cqtot*dfcq  );

  Ca->caup 	= Ca->caup + ht*( jup - jtr*C->vrel/C->vup );


  // fractional occupancies
  Ca->ftc  	= Ca->ftc  + ht* dftc;
  Ca->fcmi  	= Ca->fcmi + ht* dfcmi;
  Ca->fcms  	= Ca->fcms + ht* dfcms;
  Ca->fcq  	= Ca->fcq  + ht* dfcq;
  Ca->ftmc  	= Ca->ftmc + ht* dftmc;
  Ca->ftmm  	= Ca->ftmm + ht* dftmm;

#if 1
  // save ca data
  {
    extern REAL t;
    static REAL tgt=0.;
    static FILE *foca;
    static int first2=1;

    // open file once
    for( ;first2;first2=0) foca = fopen( "ca.dat", "w" );

    if( t>=tgt )
      {
	tgt += 1e-3;	// [s]
	fprintf( foca, "%g\t%g %g %g %g %g %g %g %g %g %g   %g\t%g\t%g\t%g\n", t,  
		 Ca->ftc, Ca->ftmc, Ca->ftmm, Ca->fcmi, Ca->fcms, Ca->fcq, Ca->cai, Ca->casub, Ca->caup, Ca->carel,  icatotal,  jrel*C->vrel/C->vsub,-jcadiff, -cmtot*dfcms  );

      }
  }//save_ca_date
#endif

  //  return casub_new;
  return Ca->cai;
} /** ca_intra **/


