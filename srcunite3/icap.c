/* icap as in Demir etc.
 */
//Rubin, December 2004
// icapmax -> variable
//Rubin, October 2007


#include "sa.h"

//#define icapmax 0.02869 /*nA*/ //ORIG
//#define icapmax (cm* 0.6*0.0011*1e6)   /*nA*/ //test
//#define icapmax 0.08 /*nA*/ //test for ko=8.1

REAL
icap(struct Cpar *C, struct Caintra_state *Ca )
{
  //  return icapmax/(1+(0.002/cai)); //test
  return C->icapmax/(1+(0.0004/Ca->casub));
}/** icap **/
