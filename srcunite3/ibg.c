// calc Ibg, background Na, K, Ca 

#include "sa.h"

//////// Na ////////////
REAL
ibgna(struct State *S, struct State *Sn, REAL ht,  struct Cpar *C )
{
  return C->gbna*(S->E -C->ENA);

} /** ibgna **/


//////// K //////////////
REAL
ibgk(struct State *S, struct State *Sn, REAL ht,  struct Cpar *C )
{
  return C->gbk*(S->E -C->EK);

} /** ibgk **/


///////// Ca /////////////
REAL
ibgca(struct State *S, struct State *Sn, REAL ht,  struct Cpar *C, struct Caintra_state *Ca )
{
  return C->gbca*(S->E -C->ECA);

} /** ibgca **/
