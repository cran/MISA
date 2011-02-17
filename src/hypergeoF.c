/*********************************************************
* hypergeoF.c	Last modified 02/04/11                   *
* Imported from BAS 0.93                                 *
*                                                        *
* Description: Support functions for bayesglm_fit().     *
*                                                        *
* Author: Merlise Clyde                                  *
*********************************************************/

extern double hyp2f1(double, double, double, double);

void hypergeometric2F1(a, b, c, x, y)
     double *y, *x, *a, *b, *c;
{

 *y = hyp2f1(*a, *b, *c, *x);
}

 
