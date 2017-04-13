/*
 * feasible_ls.cpp
 *
 *  Created on: Jan 24, 2016
 *      Author: zhou
 */
#include <set>
#include <vector>
#include <algorithm>
#include "mbbp.hpp"
using namespace std;

//extern int sg_comp_id;
//extern int sg_is_complement;
extern int crt_num_vertex[PNUM_PARTIE];
extern int crt_num_sum_ve;
extern int *crt_v_edges_count;
extern int **crt_v_adj_vertex;
extern float para_tt;

static RandAccessList *cur_X_sol;
static RandAccessList *cur_Y_sol;
static RandAccessList *sol_X_nbs; // the neighbors of X
static RandAccessList *sol_Y_nbs;// the neighbors of Y

//static int *tmp_lst;
static int *is_in_c;
static int *deg;
//static int *mark_sol;
//static int n_diff;
//static int *frequency;

static int *tabu_add;
static int iter;

/*Record the best solution, hit time*/
static clock_t fs_start_time;
static clock_t best_end_time;
static int best_iter;
static int best_up_size;	/*Objective function*/
static int best_down_size;
static int best_balanced_size;
static int *best_clique;

#define SG_Balanced_Cur_Size() (min(cur_X_sol->vnum, cur_Y_sol->vnum))

static void init_data(){
	cur_X_sol = ral_init(crt_num_sum_ve);
	cur_Y_sol = ral_init(crt_num_sum_ve);

	sol_X_nbs = ral_init(crt_num_sum_ve);
	sol_Y_nbs = ral_init(crt_num_sum_ve);
	is_in_c = new int[crt_num_sum_ve];
	deg = new int[crt_num_sum_ve];

	tabu_add = new int[crt_num_sum_ve];
	memset(tabu_add, 0, sizeof(int) * crt_num_sum_ve);
	iter = 0;

	//tmp_lst = new int[num_sum_ve];

	fs_start_time = clock();
	best_clique = new int[crt_num_sum_ve];

	memset(is_in_c, 0, sizeof(int) * crt_num_sum_ve);
	memset(deg, 0, sizeof(int) * crt_num_sum_ve);

//	mark_sol = new int[crt_num_sum_ve];
//	memset(mark_sol, 0, sizeof(int) *crt_num_sum_ve);
//	n_diff = 0;
//	frequency = new int[crt_num_sum_ve];
//	memset(frequency, 0, sizeof(int) * crt_num_sum_ve);
}


static void release_data(){

	ral_release(cur_X_sol);
	ral_release(cur_Y_sol);
	ral_release(sol_X_nbs);
	ral_release(sol_Y_nbs);

	delete[] is_in_c;
	delete[] deg;
	//delete[] tmp_lst;

	delete[] best_clique;
	delete[] tabu_add;
//	delete[] mark_sol;
//	delete[] frequency;
}

static void print_vector(vector<int>* vec){
	int max_print_size = min(10, (int)vec->size());
	printf(" [%d] : ", vec->size());
	for (int i = 0; i < max_print_size; i++){
			printf("%d ", vec->at(i));
	}
	if (max_print_size == (int)vec->size())
		printf("\n");
	else
		printf("...\n");
}

static void record_current(){
	//assert(cur_down_sol->vnum == cur_X_sol->vnum);
	best_end_time = clock();
	best_iter = iter;
	best_up_size = cur_X_sol->vnum;
	best_down_size = cur_Y_sol->vnum;
	best_balanced_size = SG_Balanced_Cur_Size();
	memcpy(best_clique, cur_X_sol->vlist, cur_X_sol->vnum * sizeof(int));
	memcpy(best_clique + cur_X_sol->vnum, cur_Y_sol->vlist, cur_Y_sol->vnum * sizeof(int));
}


static void add_vertex(int v){
	assert(!is_in_c[v]);
	is_in_c[v] = 1;
	if (v < crt_num_vertex[0]){
		ral_add(cur_X_sol, v);
		if (deg[v] > 0)
			ral_delete(sol_Y_nbs, v);
	}else{
		ral_add(cur_Y_sol, v);
		if (deg[v] > 0)
			ral_delete(sol_X_nbs, v);
	}

	for (int i = 0; i < crt_v_edges_count[v]; i++){
		int vadj = crt_v_adj_vertex[v][i];
		deg[vadj]++;
		if (!is_in_c[vadj] && deg[vadj] == 1){
			if (v < crt_num_vertex[0])
				ral_add(sol_X_nbs, vadj);
			else
				ral_add(sol_Y_nbs, vadj);
		}
	}
}

static void remove_vertex(int v){
	assert(is_in_c[v]);
	is_in_c[v] = 0;
	if (v < crt_num_vertex[0]){
		ral_delete(cur_X_sol, v);
		if (deg[v] > 0)
			ral_add(sol_Y_nbs, v);
	}else{
		ral_delete(cur_Y_sol, v);
		if (deg[v] > 0)
			ral_add(sol_X_nbs, v);
	}
	for (int i = 0; i < crt_v_edges_count[v]; i++){
		int vadj = crt_v_adj_vertex[v][i];
		deg[vadj]--;
		if (!is_in_c[vadj] && deg[vadj] == 0){
			if (v < crt_num_vertex[0])
				ral_delete(sol_X_nbs, vadj);
			else
				ral_delete(sol_Y_nbs, vadj);
		}
	}
}

static void print_cur_sol(){
	printf("up side:");
	for (int i = 0; i < cur_X_sol->vnum; i++){
		printf("%d ", cur_X_sol->vlist[i]);
	}
	printf("\ndown side:");
	for (int i = 0; i < cur_Y_sol->vnum; i++){
		printf("%d ", cur_Y_sol->vlist[i]);
	}
	printf("\n");

//	printf("Degree:")
//	for (int v = 0; v < num_sum_ve; v++){
//		printf("%d[%d] ",v, deg[v]);
//	}
//	printf("\n");
}

static void restart_ls(){
	ral_clear(cur_X_sol);
	ral_clear(cur_Y_sol);
	ral_clear(sol_X_nbs);
	ral_clear(sol_Y_nbs);
	memset(is_in_c, 0, sizeof(int) * crt_num_sum_ve);
	memset(deg, 0, sizeof(int)* crt_num_sum_ve);
//	memset(mark_sol, 0, sizeof(int) *crt_num_sum_ve);
//	n_diff = 0;
}

/*Generate init balanced solution*/
static void balanced_init_solution(){
	RandAccessList *cand_lst = NULL;
	RandAccessList *con_lst = NULL;
	int part = rand() % 2;
	vector<int> vec_add;

	assert(crt_num_vertex[0] > 0 && crt_num_vertex[1] > 0);
	int startv = -1;
	if (part == 0){
		startv = rand() % crt_num_vertex[0];
	}else{
		startv = crt_num_vertex[0] + rand() % crt_num_vertex[1];
	}
	add_vertex(startv);
	while(1){
		//switch to another side;
		part = 1 - part;
		con_lst = part == 0? cur_Y_sol: cur_X_sol;
		cand_lst = part == 0 ? sol_Y_nbs: sol_X_nbs;
		vec_add.clear();
		for (int i = 0; i < cand_lst->vnum; i++){
			int v = cand_lst->vlist[i];
			if (deg[v] == con_lst->vnum){
				vec_add.push_back(v);
			}
		}
		if (vec_add.empty())
			break;
		int vrand = vec_add[rand()%vec_add.size()];
		add_vertex(vrand);
	}
//	for (int i = 0; i < cur_X_sol->vnum; i++){
//		mark_sol[cur_X_sol->vlist[i]] = 1;
//	}
//	for (int i = 0; i < cur_Y_sol->vnum; i++){
//		mark_sol[cur_Y_sol->vlist[i]] = 1;
//	}
//	n_diff = 0;

}


static void push_with_tabu(int v, int tt_base){
	vector<int> vec_drop;
	RandAccessList *con_lst = v < crt_num_vertex[0] ? cur_Y_sol: cur_X_sol;
	RandAccessList *ve_lst =  v < crt_num_vertex[0] ? cur_X_sol: cur_Y_sol;

	assert(!is_in_c[v]);
	add_vertex(v);
//	if (mark_sol[v] != 1){
//		n_diff++;
//	}
	tabu_add[v] = iter;
//	printf("Push %d ",v);
	if (deg[v] < con_lst->vnum){
		int idx = 0;
		while (idx < con_lst->vnum){
			int ve = con_lst->vlist[idx];
			if (deg[ve] == ve_lst->vnum - 1){
				remove_vertex(ve);
//				if (mark_sol[v] == 0){
//					n_diff--;
//				}
				tabu_add[ve] = iter + max(7, (int)(para_tt * (rand()%tt_base)));
			}else if (deg[ve] == ve_lst->vnum){
				idx++;
			}else{
				printf("ERROR in push %d, iter %d\n", ve, iter);
			}
		}
	}
//	printf("\n");
}

void decompose(RandAccessList *nbs_lst, RandAccessList *con_lst, RandAccessList *ve_lst, vector<int> &M1,
		vector<int> &M2, vector<int> &M3){
	for (int i = 0; i < nbs_lst->vnum; i++){
		int v = nbs_lst->vlist[i];
		int n_conflict = con_lst->vnum - deg[v];
		int delta = 0;
		if (ve_lst->vnum > con_lst->vnum){
			delta = -n_conflict;
		}else{
			delta = min(1, con_lst->vnum - ve_lst->vnum - n_conflict);
		}
		/*We only interested in the add and switch case*/
		if (n_conflict <= 1){
			if (delta == 1 && (tabu_add[v] <= iter || SG_Balanced_Cur_Size() + 1 > best_balanced_size)){
				M1.push_back(v);
			}else if (tabu_add[v] <= iter){
				M2.push_back(v);
			}
		}
	}
}

void recover_balance(){
	RandAccessList *longer_lst = cur_X_sol->vnum > cur_Y_sol->vnum ? cur_X_sol: cur_Y_sol;
	RandAccessList *shorter_lst = cur_X_sol->vnum < cur_Y_sol->vnum ? cur_X_sol: cur_Y_sol;
	if (cur_X_sol->vnum == cur_Y_sol->vnum )
		return;
	while (longer_lst->vnum > shorter_lst->vnum){
		int vdrop = longer_lst->vlist[rand() % longer_lst->vnum];
		remove_vertex(vdrop);
		tabu_add[vdrop] = iter + max(7, (int)(para_tt * (rand() % longer_lst->vnum)));
	}
}
/*Entrance of the local search algorithm*/
void feasible_local_search(int max_iter, int max_seconds, RT_LS *result_ls, int *total_iter, int *total_sec){
	int percent_max_iter = max(max_iter / 100, 1);
	int stop = 0;
	vector<int> vec_M1;
	vector<int> vec_M2;
	vector<int> vec_M3;
	int local_best = 0;
	//debug
//	printf("sg size %d, max iter %d\n", num_sum_ve, max_iter);
	assert(result_ls != NULL);
	init_data();
	balanced_init_solution();
	if (SG_Balanced_Cur_Size() > best_balanced_size){
		record_current();
	}
	while (!stop){
		vec_M1.clear();
		vec_M2.clear();
		vec_M3.clear();
		int randidx = -1;
		int vpush = -1;
		decompose(sol_X_nbs, cur_X_sol, cur_Y_sol, vec_M1, vec_M2, vec_M3);
		decompose(sol_Y_nbs, cur_Y_sol, cur_X_sol, vec_M1, vec_M2, vec_M3);
//		printf("ITER %d: Sol X:%d Y:%d\n",iter, cur_X_sol->vnum, cur_Y_sol->vnum);
//		printf("M1: %d, M2: %d, M3: %d, Total %d\n", vec_M1.size(),
//				vec_M2.size(), vec_M3.size(), sol_X_nbs->vnum + sol_Y_nbs->vnum);

		if (!vec_M1.empty()){
			randidx = rand() % vec_M1.size();
			vpush = vec_M1[randidx];
//			printf("In iter %d: Choose %d, from M1\n\n",iter, vpush);
			push_with_tabu(vpush,vec_M1.size());
		}else if (!vec_M2.empty()){
			randidx = rand() % vec_M2.size();
			vpush = vec_M2[randidx];
//			printf("In iter %d: Choose %d, from M2\n\n",iter, vpush);
			push_with_tabu(vpush, vec_M2.size());
		}
//		printf("");
//		print_cur_sol();
		if (cur_X_sol->vnum - cur_Y_sol->vnum < -2 || cur_X_sol->vnum - cur_Y_sol->vnum>2){
			recover_balance();
		}else if (vpush == -1){
			stop = 1;
		}
//		printf("up %d, down %d\n", cur_X_sol->vnum, cur_Y_sol->vnum);
		if (SG_Balanced_Cur_Size() > local_best){
			local_best = SG_Balanced_Cur_Size();
			if (SG_Balanced_Cur_Size() > best_balanced_size)
				record_current();
		}
		if (iter >= max_iter || (float)(clock() - fs_start_time) / CLOCKS_PER_SEC > max_seconds){
			stop = 1;
			break;
		}
		iter++;
	}
ls_stop:
	result_ls->best_balanced_size = best_balanced_size;
	memcpy(result_ls->best_clique, best_clique, sizeof(int) * (best_up_size + best_down_size));
	result_ls->best_up_size = best_up_size;
	result_ls->best_down_size = best_down_size;
	result_ls->best_end_time = best_end_time;
	result_ls->best_iter = best_iter;

	*total_iter = iter;
	*total_sec = (float)(clock() - fs_start_time) / CLOCKS_PER_SEC;
	//result_ls->total_iter = total_iter;
	release_data();
	return;
}

