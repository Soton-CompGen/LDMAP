#include "allass.h"
void fill_gai(int n11, int n12, int n13, int n21, int n22, int n23, int n31, int n32, int n33, double kb1, double kb2, char locus1[20], char locus2[20], double freq1, double freq2)
{
double n;
n=n11+n12+n13+n21+n22+n23+n31+n32+n33;

if(g_aiPtr==NULL)
  {
      aistartPtr=(aisPtr)malloc(sizeof(ais));
      aistartPtr->nextPtr=NULL;
      g_aiPtr=aistartPtr;
      g_aiPtr->kb1=kb1;
      g_aiPtr->kb2=kb2;
      g_aiPtr->n=n;
      g_aiPtr->aitab[0]=n11; 
      g_aiPtr->aitab[1]=n12; 
      g_aiPtr->aitab[2]=n13; 
      g_aiPtr->aitab[3]=n21; 
      g_aiPtr->aitab[4]=n22; 
      g_aiPtr->aitab[5]=n23; 
      g_aiPtr->aitab[6]=n31; 
      g_aiPtr->aitab[7]=n32; 
      g_aiPtr->aitab[8]=n33; 
      strcpy(g_aiPtr->locus1, locus1);
      strcpy(g_aiPtr->locus2, locus2);
      g_aiPtr->freq1=freq1;
      g_aiPtr->freq2=freq2;
  }


else  /* If not the first record in the ais structure */
  {
      g_aiPtr->nextPtr=(aisPtr)malloc(sizeof(ais));
      g_aiPtr=g_aiPtr->nextPtr;
      g_aiPtr->kb1=kb1;
      g_aiPtr->kb2=kb2;
      g_aiPtr->n=n;
      g_aiPtr->aitab[0]=n11; 
      g_aiPtr->aitab[1]=n12; 
      g_aiPtr->aitab[2]=n13; 
      g_aiPtr->aitab[3]=n21; 
      g_aiPtr->aitab[4]=n22; 
      g_aiPtr->aitab[5]=n23; 
      g_aiPtr->aitab[6]=n31; 
      g_aiPtr->aitab[7]=n32; 
      g_aiPtr->aitab[8]=n33; 
      strcpy(g_aiPtr->locus1, locus1);
      strcpy(g_aiPtr->locus2, locus2);
      g_aiPtr->freq1=freq1;
      g_aiPtr->freq2=freq2;

  }

  }
