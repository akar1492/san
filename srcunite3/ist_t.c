// ist -- Sustained inward current
//
// Am J Physiol Heart Circ Physiol 283: H2074­H2101, 2002;
// Dynamical description of sinoatrial node pacemaking: improved mathematical model for primary pacemaker cell
// YASUTAKA KURATA, ICHIRO HISATOME, SUNAO IMANISHI, and TOSHISHIGE SHIBAMOTO
//
// use tables to speed-up 

#include "sa.h"


////////// tau and ss 
void
tss_qa( REAL ee, REAL *tau, REAL *ss )
{
  REAL a, b;

  a = 1/( 0.15*exp(-ee/11) + 0.2*exp(-ee/700) );

  b = 1/( 16*exp(ee/8) + 15*exp(ee/50) );
    
  *tau  = 1./(a+b);
  *ss   = 1./(1+exp(-(ee+57)/5));
}///tss_qa

////////// tau and ss 
void
tss_qi( REAL ee, REAL *tau, REAL *ss )
{
  REAL a, b;

  a = 0.1504/(3100*exp(ee/13)+700*exp(ee/70));
  b = 0.1504/(95*exp(-ee/10)+50*exp(-ee/700)) + 0.000229/(1+exp(-ee/5));
    
  *tau  = 1./(a+b);
  *ss   = a* *tau;
}///tss_qi


///////////////////////////////////////////////////////////////////////
REAL
ist_t(struct State *S, struct State *Sn, REAL ht, struct Table T, struct Cpar *C )
{
  static int first=1; 
  static REAL Est; 
  static REAL qat[NTABLE],qass[NTABLE], qit[NTABLE],qiss[NTABLE];
  REAL qat_t, qass_t, qit_t, qiss_t;
  // gkr should be here, but it proved to be varied from center out

  // first time
  for( ; first; first=0 )
    { 
      float pnak = 0.12;
      Est = 37.4;	// CONST
      printf("### Est=%g\n", Est);

      //fill tables
      { int i; REAL ee;
      for( i=-1; ++i<NTABLE; )
	{
	  ee = ETMIN+ (ETMAX-ETMIN)*i/(NTABLE-1);
	  // tabs
	  tss_qa( ee, qat+i, qass+i );
	  tss_qi( ee, qit+i, qiss+i );
	}//fori
      }//fill
    }//forfirst

  // estimate table values
  {
    register int iii=T.n1;
    qat_t  = qat[iii] + T.frac*(qat[iii+1]-qat[iii]);
    qass_t = qass[iii] + T.frac*(qass[iii+1]-qass[iii]);
    qit_t  = qit[iii] + T.frac*(qit[iii+1]-qit[iii]);
    qiss_t = qiss[iii] + T.frac*(qiss[iii+1]-qiss[iii]);
  }
 
    Sn->qa = qass_t - (qass_t - S->qa)*exp(-ht/qat_t);	// solve eq
    Sn->qi = qiss_t - (qiss_t - S->qi)*exp(-ht/qit_t);	// solve eq

#define gst (0.015 *C->cm*1e6) 		// [ns/pF]*[muF*1e6] 
    return gst*S->qa*S->qi*(S->E -Est);

} /** ist_t **/
