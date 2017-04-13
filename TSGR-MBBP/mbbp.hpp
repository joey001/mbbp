/*
 * mbbp.hpp
 *
 *  Created on: Dec 3, 2015
 *      Author: zhou
 */

#ifndef MBBP_HPP_
#define MBBP_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <queue>

using namespace std;

#define MAX_VAL 999999999
#define PNUM_PARTIE 2

#define REDUCE_FUNCTION_ON 0

/*Define the structure*/
typedef struct ST_RandAccessList{
	int *vlist;
	int *vpos;
	int vnum;
	int capacity;
}RandAccessList;

typedef struct ST_RT_LS_Search{
	//clock_t start_time;
	clock_t best_end_time;
	int best_iter;
	int best_up_size;	/*Objective function*/
	int best_down_size;
	int best_balanced_size;
	int *best_clique;
}RT_LS;

extern void copy_RT_LS(RT_LS *dest, RT_LS *src);
void feasible_local_search(int max_iter, int max_seconds, RT_LS *result_ls, int *total_iter, int *total_sec);
void enum_search(int lbound, int max_second, RT_LS *rt, int *total_iter, int* total_seconds, int* is_opt);
/***************utils*********************/
RandAccessList* ral_init(int capacity);
void ral_add(RandAccessList *ral, int vid);
void ral_delete(RandAccessList *ral, int vid);
void ral_clear(RandAccessList *ral);
void ral_release(RandAccessList *ral);
int cmpfunc (const void * a, const void * b);
void ral_showList(RandAccessList *ral, FILE *f);

void shuffle(int *randlist, int len, int randlen);


#endif /* MBBP_HPP_ */
