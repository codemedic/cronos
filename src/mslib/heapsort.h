#ifndef HEAPSORT_H
#define HEAPSORT_H

#include "abrand.h"

template<class T> void sort_objects(int n, T *object,int *indices)
  {
    int l,j,ir,i;
    T tempobj;
    int tempindex;

    if (n<2) return;		// doesn't work if only one observation
    // but we don't need to sort then anyway

    if (indices!=NULL)
      for (i=0 ; i<n ; ++i)
	indices[i] = i;     // start off with identity permutation

    l = (n >> 1)+1;
    ir = n;

    for (;;) {
      if (l>1)		// i.e. we are still in hiring phase
	{ 
	  --l;
	  if (indices!=NULL)
	    tempindex = indices[l-1];
	  tempobj = object[l-1];
	}
      else {
	tempobj=object[ir-1];
	if (indices!=NULL)
	  tempindex=indices[ir-1];
	object[ir-1]=object[0];
	if (indices!=NULL)
	  indices[ir-1]=indices[0];
	if (--ir == 1) {	// we are finished
	  object[0]=tempobj;
	  if (indices!=NULL)
	    indices[0]=tempindex;
	  return;
	}
      }
      i=l;		
      j=l<<1;
      while (j <= ir)
	{
	  if ((j < ir) && (object[j-1]<object[j])) 
	    ++j;	// j is better underling
	  if (tempobj<object[j-1]) {		// demote tempobj
	    object[i-1]=object[j-1];
	    if (indices!=NULL)
	      indices[i-1]=indices[j-1];
	    j += (i=j);
	  }
	  else j=ir+1;			// finished sifting
	}
      object[i-1]=tempobj;
      if (indices!=NULL)
	indices[i-1]=tempindex;
    }
  }



#endif
