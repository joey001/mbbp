/*
 * utils.hpp
 *
 *  Created on: Mar 7, 2016
 *      Author: zhou
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_
#include <stdio.h>
/*Define the structure*/
typedef struct ST_RandAccessList{
	int *vlist;
	int *vpos;
	int vnum;
	int capacity;
}RandAccessList;

/***************utils*********************/
RandAccessList* ral_init(int capacity);
void ral_add(RandAccessList *ral, int vid);
void ral_delete(RandAccessList *ral, int vid);
void ral_clear(RandAccessList *ral);
void ral_release(RandAccessList *ral);
int cmpfunc (const void * a, const void * b);
void ral_showList(RandAccessList *ral, FILE *f);

void shuffle(int *randlist, int len, int randlen);


#endif /* UTILS_HPP_ */
