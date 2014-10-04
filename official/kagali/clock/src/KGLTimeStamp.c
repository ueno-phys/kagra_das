/*
 *  Copyright (C) 2014 Koh Ueno
 *
 *
 *  This program is used to roughly estimate the time for a computational process.
 *  You need to call this function twice: just before and after the process,
 *  and calculate the difference of the two outputs (in units of seconds).
 *
 *  If you edit this code, add the date and your name (also your e-mail address)
 *  to the following list and describe the difference from the current version.
 *
 *  Creation Date and Author: 
 *  2014-06-26 ; First version by Koh Ueno (ueno@vega.ess.sci.osaka-u.ac.jp)
 *
 */

#include <kagali/KGLStdlib.h>
#include <time.h>

void KGLTimeStamp( //begin{proto}
    clock_t clock0
    ) //end{proto}
{
  printf("                                                                  %f\n",
	 (double)(clock()-clock0)/CLOCKS_PER_SEC);
  return;
}
