// main module to cal SA (central)
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "sa.h"

#define TABLE 1
#define ACH 1 
#define WRITEIC 0               // (at the end) write new IC to a file);
#define IONS 1

// extern subs
extern REAL ibgna(struct State *, struct State *, REAL, struct Cpar *);
extern REAL ibgca(struct State *, struct State *, REAL, struct Cpar *, struct Caintra_state *);
extern REAL ibgk(struct State *, struct State *, REAL, struct Cpar *);
extern REAL ical_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL icalach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL icat_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL icatach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL ikr_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL iks_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL ikach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL if_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Is *);
extern REAL ifach_t(struct State *, struct State *, REAL, struct Table, struct Cpar *, struct Is *);
extern REAL ina_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL isus_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL ito_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);
extern REAL ist_t(struct State *, struct State *, REAL, struct Table, struct Cpar *);

extern REAL ical(struct State *, struct State *, REAL, struct Cpar *);
extern REAL icat(struct State *, struct State *, REAL, struct Cpar *);
extern REAL if_(struct State *, struct State *, REAL, struct Cpar *);
extern REAL ikr(struct State *, struct State *, REAL, struct Cpar *);
extern REAL iks(struct State *, struct State *, REAL, struct Cpar *);
extern REAL ina(struct State *, struct State *, REAL, struct Cpar *);
extern REAL inaca(struct State *, struct State *, REAL, struct Cpar *, struct Caintra_state *);
extern REAL ip(struct State *, struct State *, REAL, struct Cpar *);
extern REAL isus(struct State *, struct State *, REAL, struct Cpar *);
extern REAL ito(struct State *, struct State *, REAL, struct Cpar *);

extern REAL icap(struct Cpar *, struct Caintra_state *);
extern REAL ca_intra( REAL, REAL, struct Cpar *, struct Caintra_state *);
extern unsigned long clocks( char );

REAL basic_pars( REAL vthr, REAL vo1, REAL vn1, REAL t, REAL ht );
void settype(struct Cpar *);

REAL t, ht=0.01e-3;		// [s]

REAL itotal=0.;			// total current = sum all currents

struct State Stbuf;	// buffer state  



REAL istim;			// [na] stimul current -- variable
 
REAL trun=50;//30;		// [s] run time 


void settype(struct Cpar *C){
  C->cm=cmC+C->ctype*(cmP-cmC);
  C->gna=gnaC+C->ctype*(gnaP-gnaC);
  C->gto=gtoC+C->ctype*(gtoP-gtoC);
  C->gsus=gsusC+C->ctype*(gsusP-gsusC);
  C->gkr=gkrC+C->ctype*(gkrP-gkrC);
  C->gks=gksC+C->ctype*(gksP-gksC);
  C->gfna=gfnaC+C->ctype*(gfnaP-gfnaC);
  C->gfk=gfkC+C->ctype*(gfkP-gfkC);
  C->gbna=gbnaC+C->ctype*(gbnaP-gbnaC);
  C->gbca=gbcaC+C->ctype*(gbcaP-gbcaC);
  C->gbk=gbkC+C->ctype*(gbkP-gbkC);
  C->ipss=(ipssC+C->ctype*(ipssP-ipssC));
  C->knaca=knacaC+C->ctype*(knacaP-knacaC);
  C->icapmax=icapmaxC+C->ctype*(icapmaxP-icapmaxC);

  C->kachical=kachicalC+C->ctype*(kachicalP-kachicalC);
  C->kachicat=kachicatC+C->ctype*(kachicatP-kachicatC);
  C->kachikach=kachikachC+C->ctype*(kachikachP-kachikachC);

  C->vc=(0.11*C->cm*1e6);
  C->vrel=(0.0012*C->vc);
  C->vup=(0.0116*C->vc);
  C->vsub=(0.01*C->vc);
  C->vi=(0.46*C->vc-C->vsub);
	
  C->gcal=gcalC+C->ctype*(gcalP-gcalC);
  C->gcat=gcatC+C->ctype*(gcatP-gcatC);
  
}
//////////////////////////////////////////////////////////////////////////
// Total Current
//////////////////////////////////////////////////////////////////////////
REAL
fitotal(struct State *St, struct State *Stn,REAL ht, struct Cpar *Cp, struct Is *I, struct Table T, struct Caintra_state *Ca)
{
  I->fibgna = ibgna(St, Stn, ht, Cp);
  I->fibgca = ibgca(St, Stn, ht, Cp,Ca);
  I->fibgk = ibgk(St, Stn, ht, Cp);
  if(1/* t>0 */) clocks('s');

#if TABLE  

#if ACH
  I->fical = icalach_t(St, Stn, ht, T, Cp);
  I->ficat = icat_t(St, Stn, ht, T, Cp);
  I->fif = ifach_t(St, Stn, ht, T, Cp, I);
  I->fikach = ikach_t(St, Stn, ht, T, Cp);
#else //noACH
  I->fical = ical_t(St, Stn, ht, T, Cp);
  I->ficat = icat_t(St, Stn, ht, T, Cp);
  I->fif = if_t(St, Stn, ht, T, Cp, I);
#endif //ACH
  I->fikr = ikr_t(St, Stn, ht, T, Cp);
  I->fiks = iks_t(St, Stn, ht, T, Cp);
  I->fisus = isus_t(St, Stn, ht, T, Cp);
  I->fito = ito_t(St, Stn, ht, T, Cp);
  I->fina = ina_t(St, Stn, ht, T, Cp);
#endif //TABLE

   if(1/* t>0 */) clocks('f');

  I->finaca = inaca(St, Stn, ht, Cp,Ca);
  I->fip = ip(St, Stn, ht, Cp);
  I->ficap = icap(Cp,Ca);
//printf("!!!!!!!!!fibgna=%g\n fibgca=%g\n fibgk=%g\n fical=%g\n ficat=%g\n fif=%g\n fikr=%g\n fiks=%g\n fina=%g\n finaca=%g\n fip=%g\n fisus=%g\n fito=%g\n fikach=%g\n fist=%g\n ficap=%g\n", I->fibgna,I->fibgca,I->fibgk,I->fical,I->ficat,I->fif,I->fikr,I->fiks,I->fina,I->finaca,I->fip,I->fisus,I->fito,I->fikach,I->fist,I->ficap);
  return I->fibgna+I->fibgca+I->fibgk+I->fical+I->ficat+I->fif+I->fikr+I->fiks+I->fina+I->finaca+I->fip+I->fisus+I->fito +I->fikach +I->fist +I->ficap;
}
///////////////////////////////////////////////////////////////////////////
// Func( t, E ) is an rh for eq dE/dt = Func
// t = ht
// E = E
///////////////////////////////////////////////////////////////////////////
REAL Func( REAL t, struct State *Stn, struct State *St, struct Cpar *Cp, struct Is *I, struct Table T,  struct Caintra_state *Ca){
  Stbuf = *St; 		// state struct: start from same IC
  return -(istim+ fitotal(&Stbuf, Stn, t, Cp, I, T,Ca))/Cp->cm;  	// RH part +istim
}/** Func **/

///////////////////////////////////////////////////////////////////////////
// calculate Nernst potentials using intra- and extra- concentrations 
///////////////////////////////////////////////////////////////////////////
REAL
potentials()
{
  Cp1.ENA 	= RTF* log(nao/Cp1.nai);
  Cp1.EK 		= RTF* log(ko/Cp1.ki);
  Cp1.ECA =	0.5*RTF*log(cao/Cp1.cai); // zca=2;


} /** potentials **/

///////////////////////////////////////////////////////////////////////////
main(int argc, char **argv){

   Cp1.ctype=atof(argv[1]);
   settype(&Cp1);
//   printf("#ctype=%g\n", ctype );
//     trun=atof(argv[2]);
//   printf("#trun=%g\n", trun );
//  ctype=1.0;


//    printf("###cm\t%g\n gna\t%g\n gto\t%g\n gsus\t%g\n gkr\t%g\n gks\t%g\n gfna\t%g\n gfk\t%g\n gbna\t%g\n gbca\t%g\n gbk\t%g\n ipss\t%g\n knaca\t%g\n icapmax\t%g\n kachical\t%g\n kachicat\t%g\n kachikach\t%g\n", cm, gna, gto, gsus, gkr, gks, gfna, gfk, gbna, gbca, gbk, ipss, knaca, icapmax, kachical, kachicat, kachikach);
  


  ach=0e-8; //2.5e-8;			// consts here


//  printf("###gcal,gcat= %g %g\n", gcal, gcat );

  {
    FILE *fin = fopen( "state.dat", "r" );	// recent in this dir
    if(!fin) exit(puts("!!! Cannot open IC file"));
    fread( &Stn1, sizeof(struct State), 1, fin );


      fread( &Ca1, sizeof(struct Caintra_state), 1, fin );

      Cp1.cai = Ca1.cai;				// saved value
      fread( &Cp1.nai, sizeof(REAL), 1, fin );
      fread( &Cp1.ki, sizeof(REAL), 1, fin );
//	Cp1.ki=140;
   fclose(fin);}


 

  for( t=0; t<trun; t+=ht ){
      potentials();				// Nernst potentials
      basic_pars( -30, St1.E, Stn1.E, t, ht );	// check fronts, periods etc.



      {static float tgt=0.;
	if( t>=tgt ){
	       tgt += 0.001;	// [s]
		printf("%g\t%g\t%g %g %g %g %g %g %g %g %g %g %g %g %g\t%g %g %g %g\n", t, Stn1.E, I1.fibgna+I1.fibgca+I1.fibgk,I1.fical,I1.ficat,I1.fif,I1.fikr,I1.fiks,I1.fina,I1.finaca,I1.fip,I1.fisus,I1.fito,I1.fikach,I1.fist,I1.ficap,  Cp1.cai, Cp1.ki, Cp1.nai );
}}//output
	    
	{St1 = Stn1; mintau = 1e33;}
#if WRITEIC
	{static float twrite=0.;
		if( t>=twrite ){
			twrite+=100;
			FILE *fout = fopen( "state.dat", "w" );
			fwrite( &Stn1, sizeof(struct State), 1, fout);
			fwrite( &Ca1, sizeof(struct Caintra_state), 1, fout);
			fwrite( &Cp1.nai, sizeof(REAL), 1, fout );
			fwrite( &Cp1.ki, sizeof(REAL), 1, fout );
			fclose(fout);
		}
	}
#endif //WRITEIC

#if TABLE
      {
	int ii; double dd;
	dd = (St1.E-ETMIN)/(ETMAX-ETMIN)*(NTABLE-1); 
	ii = (int)dd;
	
	T1.n1 = ii;
	if( T1.n1<0 || T1.n1>=NTABLE ) fprintf(stderr, "!!!T1.n1=%d\n", T1.n1);
	T1.frac = dd-ii;
      }
    

#endif //TABLE
      // Euler 


       Stn1.E = St1.E + ht*(Func(ht,&Stn1, &St1, &Cp1, &I1, T1, &Ca1));
       Cp1.cai = ca_intra(ht, I1.fical+I1.ficat-2*I1.finaca+I1.ficap+I1.fibgca, &Cp1, &Ca1);
#if IONS
       Cp1.nai += ht* -1e6*(I1.fina+3.*I1.finaca+3.*I1.fip+I1.fibgna+I1.fifna)/(FRD*Cp1.vi); 
       Cp1.ki  += ht* -1e6*(-2.*I1.fip+I1.fikr+I1.fiks+I1.fikach+I1.fito+I1.fisus+I1.fifk+I1.fibgk)/(FRD*Cp1.vi);
#endif


    }//for_t
	
#if WRITEIC
	FILE *fout = fopen( "state.dat", "w" );
	fwrite( &Stn1, sizeof(struct State), 1, fout);
	fwrite( &Ca1, sizeof(struct Caintra_state), 1, fout);
	fwrite( &Cp1.nai, sizeof(REAL), 1, fout );
	fwrite( &Cp1.ki, sizeof(REAL), 1, fout );
	fclose(fout);
#endif //WRITEIC
  clocks('p');
}//main
