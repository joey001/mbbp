/*
 * ral.cpp
 *
 *  Created on: Mar 7, 2016
 *      Author: zhou
 */


#include "utils.hpp"
#include <stdlib.h>
#include <string.h>
/*****************RandAcessList*****************************/
RandAccessList* ral_init(int capacity){
	RandAccessList *ral = new RandAccessList;
	ral->vlist = new int[capacity];
	ral->vpos = new int[capacity];
	ral->vnum = 0;
	ral->capacity = capacity;
	for (int i = 0; i< capacity;i++){
		ral->vpos[i] = capacity;
	}
	return ral;
}

void ral_add(RandAccessList *ral, int vid){
//	assert(ral->vpos[vid] >= ral->vnum);
	ral->vlist[ral->vnum] = vid;
	ral->vpos[vid] = ral->vnum;
	ral->vnum++;
}
void ral_delete(RandAccessList *ral, int vid){
//	assert(ral->vpos[vid] < ral->vnum);
	int last_id = ral->vlist[ral->vnum - 1];
	int id_pos = ral->vpos[vid];
	ral->vlist[id_pos] = last_id;
	ral->vpos[last_id] = id_pos;
	ral->vnum--;
//	ral->vpos[vid] = ral->vnum; /*It is not obligatory*/
}

void ral_clear(RandAccessList *ral){
	ral->vnum = 0;
}
void ral_release(RandAccessList *ral){
	delete[] ral->vlist;
	delete[] ral->vpos;
	delete ral;
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

void ral_showList(RandAccessList *ral, FILE *f){
	fprintf(f, "Total %d: ",ral->vnum);
	int *tmp_lst= new int[ral->capacity];
	memcpy(tmp_lst, ral->vlist, ral->vnum * sizeof(int));
	qsort(tmp_lst, ral->vnum, sizeof(int), cmpfunc);
	for(int i = 0;i < ral->vnum; i++){
		fprintf(f, "%d ", tmp_lst[i]);
	}
	fprintf(f, "\n");
}
