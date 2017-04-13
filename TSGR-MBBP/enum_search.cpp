/*
 * mbc_bb.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: zhou
 */


#include <set>
#include <vector>
#include <algorithm>
#include "mbbp.hpp"
using namespace std;

extern int sg_comp_id;
extern int sg_num_vertex[PNUM_PARTIE];
extern int sg_num_sum_ve;
extern int *sg_v_edge_count;
extern int **sg_v_adj_vertex;

static int best_balance_size;
static vector<int> best_solution;
static clock_t start_time;
static clock_t find_time;
static int allowed_seconds;
static int node_cnt;
static int find_node;
static int optima;

static void print_vector(vector<int>& vec){
	int max_print_size = min(20, (int)vec.size());
	printf(" [%d] : ", (int)vec.size());
	for (int i = 0; i < max_print_size; i++){
			printf("%d ", vec[i]);
	}
	if (max_print_size == (int)(vec.size()))
		printf("\n");
	else
		printf("...\n");
}

static void bbclique(vector<int>& X, vector<int>& Y, int xsize, int ysize, vector<int>& cursol){
	node_cnt++;
	if (X.size() == 0){
		if (xsize > best_balance_size){
			best_balance_size = xsize;
			find_time = clock();
			find_node = node_cnt;
			best_solution.clear();
			best_solution.resize(cursol.size());
			copy(cursol.begin(), cursol.end(), best_solution.begin());
		}
		return;
	}
	if ((float)(clock() - start_time) / CLOCKS_PER_SEC > allowed_seconds){
		optima = 0;
		return;
	}
	while (!X.empty()){
		if (xsize + (int)X.size() <= best_balance_size){
			break;
		}
		int v = X.back();
		X.pop_back();
		vector<int> Y_new;
		for (int i = sg_v_edge_count[v]-1; i >= 0; i--){
			int vadj = sg_v_adj_vertex[v][i];
			vector<int>::iterator rt_find = find(Y.begin(), Y.end(), vadj);
			if (rt_find != Y.end()){
				Y_new.push_back(vadj);
			}
		}
		cursol.push_back(v);
		vector<int> X_new(X);
		if (xsize < ysize){
			bbclique(X_new, Y_new, xsize+1, ysize, cursol);
		}else{ //xsize == ysize
			bbclique(Y_new, X_new, ysize, xsize+1, cursol);
		}
		cursol.pop_back();
	}
}

void enum_search(int lbound, int max_second, RT_LS *rt, int *total_nodes, int* total_seconds, int *is_opt){
	vector<int> up_set;
	vector<int> down_set;
	vector<int> sol;


	start_time = clock();
	find_time = clock();
	allowed_seconds = max_second;
	node_cnt = 0;
	find_node = -1;
	best_balance_size = lbound;
	best_solution.clear();
	optima = 1;

	for (int i = sg_num_vertex[0]-1; i >= 0; i--){
		up_set.push_back(i);
	}
	for (int i = sg_num_sum_ve - 1; i >= sg_num_vertex[0]; i--){
		down_set.push_back(i);
	}
	bbclique(up_set, down_set, 0, 0, sol);

	rt->best_balanced_size = best_balance_size;
	rt->best_up_size = best_balance_size;
	rt->best_down_size = best_balance_size;
	rt->best_iter = find_node;
	rt->best_end_time = find_time;
	copy(best_solution.begin(), best_solution.end(), rt->best_clique);
	*total_seconds = (float)(clock() - start_time) / CLOCKS_PER_SEC;
	*total_nodes = node_cnt;
	*is_opt = optima;
	//debug
//	for (int i = 0; i < best_solution.size(); i++){
//		printf("%d ",best_solution[i]);
//	}
//	printf("\n");
}
