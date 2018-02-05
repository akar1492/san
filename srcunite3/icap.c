/* icap as in Demir etc.
 */
//Rubin, December 2004
// icapmax -> variable
//Rubin, October 2007


#include "sa.h"

#define icapmax 0.08 /*nA*/

REAL
icap(struct Cpar *C, struct Caintra_state *Ca )
{
  return icapmax/(1+(0.002/Ca->cai));
}/** icap **/
