#include "allass.h"

void newentry2()
{
/*Make new entry in map structure */ 
   if(g_p12==NULL)
     {
      map_startPtr2=(mapPtr2)malloc(sizeof(mapkb2));
      map_startPtr2->nextPtr=NULL;
      g_p12=map_startPtr2;
      } 
      else
      {
      g_p12->nextPtr=(mapPtr2)malloc(sizeof(mapkb2));
      g_p12=g_p12->nextPtr;
      } 
}
