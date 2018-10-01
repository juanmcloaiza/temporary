#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/*
 * This is nothing more than a trivial binary-interpolating search routine on arrays of double.
 * In the case you need a search routine for whatever data type you must modify it.
 */

int INLINE_FUNC getindex(double *table, int lower, int greatest, double *tval, int *index)
{
  double vh, maxval, minval;

  int found, half, upi, downi, res;

  downi = 0;

  maxval = table[greatest];
  minval = table[lower];

  upi = greatest;
  found = 0;

  if(*tval >= maxval)
    {
      *index = greatest;
      /**tval = maxval;*/
      res = -1;
    }
  else if(*tval <= minval)
    {
      *index = lower;
      /**tval = minval;*/
      res = -2;
    }
  else
    {
      while(found != 1)
	{
	  /* non-interpolated binary guess          
	   * half = (upi + downi) / 2;
	   */

	  /* linearly interpolate median index.
	   * this work if values are equally linearly distributed as is our case.
	   */

	  /*if((half = downi + (int)((*tval - table[downi]) * (upi - downi) / (table[upi] - table[downi]))) != downi) */
	  if((half = (upi + downi) / 2) != downi)
	    {
	      vh = table[half];

	      if(*tval < vh)
		upi = half;
	      else if(*tval > vh)
		downi = half;
	      else
		{
		  downi = half;
		  found = 1;
		}
	      if((upi - downi) == 1)
		found = 1;
	    }
	  else
	    found = 1;
	}

      res = 0;
      *index = downi;
    }
  return (res);
}
