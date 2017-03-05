//============================================================================
// Name        : mbbp_greedy.cpp
// Author      : JOEY
// Version     :
// Copyright   : GPL
// Description : This is an implementation of A. Al-Yamani et al's algorithm.
// 				https://doi.org/10.1109/TCSI.2007.907875
//============================================================================
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <queue>
#include <set>
#include <list>

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "utils.hpp"
using namespace std;

#define DEBUG_INFO 0

#define PNUM_PARTIE 2
struct Edge {
	int v1;
	int v2;
	Edge(int v1, int v2){
		this->v1 = v1;
		this->v2 = v2;
	}
};

/**************************Section PARAMETERS******************************************/
char para_file[1024] =	"/home/zhou/benchmarks/MBBP/rand_test/rand_100_0.3.clq";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/youtube-groupmemberships/out.youtube-groupmemberships"; //94,238 + 30,087
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/escorts/out.escorts";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/brunson_southern-women/out.brunson_southern-women_southern-women";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/gottron-reuters/out.gottron-reuters";
//char para_file[1024] = "g14_example.clq";
//char para_file[1024] = "t1_example.clq";
int para_max_seconds = 10800;
int para_seed = 1234567;// Dummy variable, reserved for our python script
/**********************Section Graph DATA****************************************/

/*Current graph*/
int g_num_vertex[PNUM_PARTIE];
int g_num_edge;
int g_sum_ve;
int* g_v_edges_count; //v_edge[0][i] is vertex i's number of adjacent edges in the current grpah
int** g_v_adj_vertex;	//v_adj_vertex[i][j] is vertex i's jth adjacet vertex
int g_min_deg;

/*The complementent graph*/
int cmpg_num_edge;
int* cmpg_cnt_edges;
int** cmpg_v_adj_vtx;

#define LIdx2Part(idx) ((idx)<g_num_vertex[0]? 0 : 1)
#define LIdx2PartIdx(idx) ((idx)<g_num_vertex[0]? (idx) : (idx)-g_num_vertex[0])
#define PIdx2LIdx(part,idx) ((part)*g_num_vertex[0]+(idx))
#define MAX_VAL 9999999

/********************Section Algorithm DATA************************************/
clock_t start_time;
clock_t total_time;
clock_t best_time;

long long node_cnt;
int is_optima;
long long best_node_id;
int max_allowed_seconds;
int best_half_size;
int act_part;

int* sorted_act_set;
vector<int> best_solution;

/**********************Declare hopcrof-karp algorithm *******************************/
//extern void hk_init(int n1, int n2);
//extern int hk_maxMatching();
//extern void hk_addEdge(int u, int v);
int asc_idx_cmp(const void *a, const void *b) {
	return *(int*)a - *(int*)b;
}

int load_konect_instance(char *filename){
	ifstream infile(filename);
	if (infile == NULL){
		fprintf(stderr,"Can not open file %s\n", filename);
		return 0;
	}
	char line[1000];
	int vid1, vid2;

	/*ignore all the comment lines*/

	g_num_vertex[0] = 0;
	g_num_vertex[1] = 0;

	vector<Edge> *alledges= new vector<Edge>();
	/*Read the original graph, record the maximum id of vertices
	 * as the maximum number of the vertices*/
	while (infile.getline(line,  1000)){
		//cout << id1 << " " << id2 << endl;
		if (line[0] == '%')
			continue;
		else
			sscanf(line, "%d %d", &vid1, &vid2);
		//update the largest vertex number
		if (vid1 > g_num_vertex[0])
			g_num_vertex[0] = vid1;
		if (vid2 > g_num_vertex[1]){
			g_num_vertex[1] = vid2;
		}
		alledges->push_back(Edge(vid1,vid2));
	}
	//printf("Vertices %d %d, Edges: %d\n", g_num_vertex[0], g_num_vertex[1], (int)alledges->size());

	/*Allocate space for  current graph*/
	g_sum_ve = g_num_vertex[0] + g_num_vertex[1];

	/*estimate adjacent edge numbers for each nodes*/
	int *estimate_edge_cnt = new int[g_sum_ve];
	memset(estimate_edge_cnt, 0 , sizeof(int) * g_sum_ve);
	/*count the number of adjacent edges for each vertex*/
	for (vector<Edge>::iterator itr = alledges->begin();
			itr != alledges->end();
			itr++){
		vid1 = ((Edge)(*itr)).v1 - 1;
		vid2 = ((Edge)(*itr)).v2 - 1;
		//cout << id1 <<" "<< id2 << endl;
		int idx1 = PIdx2LIdx(0, vid1);
		int idx2 = PIdx2LIdx(1, vid2);
		estimate_edge_cnt[idx1]++;
		estimate_edge_cnt[idx2]++;
	}
	/*allocate space for the adjacent vertex list of each vertex*/
	g_v_edges_count = new int[g_sum_ve];
	g_v_adj_vertex = new int*[g_sum_ve];
	for (int idx = 0; idx < g_sum_ve; idx++){
		if (estimate_edge_cnt[idx] > 0)
			g_v_adj_vertex[idx] = new int[estimate_edge_cnt[idx]];
		else
			g_v_adj_vertex[idx] = (int*)0; //for secure
	}

	/*Build up the adjacent vertex list, remove the duplicated edges*/
	memset(g_v_edges_count, 0, sizeof(int) * g_sum_ve);
	g_num_edge = 0;
	for (vector<Edge>::iterator itr = alledges->begin();
				itr != alledges->end();
				itr++){
			vid1 = ((Edge)(*itr)).v1 - 1;
			vid2 = ((Edge)(*itr)).v2 - 1;
			//cout << id1 <<" "<< id2 << endl;
			int idx1 = PIdx2LIdx(0, vid1);
			int idx2 = PIdx2LIdx(1, vid2);
			int* rt_find1 = find(g_v_adj_vertex[idx1],
					g_v_adj_vertex[idx1] + g_v_edges_count[idx1], idx2);
			if (rt_find1 == g_v_adj_vertex[idx1] + g_v_edges_count[idx1]){
				g_v_adj_vertex[idx1][g_v_edges_count[idx1]++] = idx2;
				g_num_edge ++;

			}
			int* rt_find2 = find(g_v_adj_vertex[idx2],
					g_v_adj_vertex[idx2]+g_v_edges_count[idx2], idx1);
			if (rt_find2 == g_v_adj_vertex[idx2]+g_v_edges_count[idx2]){
				g_v_adj_vertex[idx2][g_v_edges_count[idx2]++] = idx1;
			}
		}
	for (int i = 0; i < g_sum_ve; i++){
		qsort(g_v_adj_vertex[i], g_v_edges_count[i], sizeof(int), asc_idx_cmp);
	}

	delete[] estimate_edge_cnt;
	delete alledges;
	return 1;
}
/*Load random generalized instances*/
int load_rand_instance(char* filename) {
	ifstream infile(filename);
	char line[1024];
	char tmps1[1024];
	char tmps2[1024];
	vector<Edge*>* all_edges;

	if (infile == NULL) {
		fprintf(stderr, "Can not find file %s\n", filename);
		return 0;
	}
	/*Read the comments*/
	infile.getline(line, 1024);
	while (line[0] != 'p')
		infile.getline(line, 1024);
	sscanf(line, "%s %s %d %d %d", tmps1, tmps2, &g_num_vertex[0], &g_num_vertex[1],
			&g_num_edge);

	/*The current graph*/
	all_edges = new vector<Edge*>();
	g_sum_ve = g_num_vertex[0] + g_num_vertex[1];
	g_v_edges_count = new int[g_sum_ve];
	//g_bmat = vector<bitarray>();

	fill(g_v_edges_count, g_v_edges_count + g_sum_ve, 0);
	while (infile.getline(line, 1024)) {
		int v1, v2;
		if (strlen(line) == 0)
			continue;
		if (line[0] != 'e')
			fprintf(stderr, "Line format error %s", line);
		sscanf(line, "%s %d %d", tmps1, &v1, &v2);
		v1--, v2--;
		all_edges->push_back(new Edge(v1,v2));

		g_v_edges_count[v1]++;
		g_v_edges_count[v2 + g_num_vertex[0]]++;
	}
	assert(all_edges->size() == (int )g_num_edge);

	//crt_v_component_id = new int[crt_num_sum_ve];
	g_v_adj_vertex = new int*[g_sum_ve];
	g_min_deg = MAX_VAL;
	vector<int> *init_part = new vector<int>(g_sum_ve);
	for (int i = 0; i < g_sum_ve; i++) {
		init_part->at(i) = i;
		//crt_v_component_id[i] = 0;
		g_v_adj_vertex[i] = new int[g_v_edges_count[i]];
		if (g_v_edges_count[i] < g_min_deg)
			g_min_deg = g_v_edges_count[i];
	}
	//connect_components.push_back(init_part);
	int *adj_size = new int[g_sum_ve];
	memset(adj_size, 0, sizeof(int) * g_sum_ve);
	for (int i = 0; i < (int) all_edges->size(); i++) {
		Edge *pe = all_edges->at(i);
		int v1 = pe->v1, v2 = g_num_vertex[0] + pe->v2;
		g_v_adj_vertex[v1][adj_size[v1]++] = v2;
		g_v_adj_vertex[v2][adj_size[v2]++] = v1;
		delete pe;
	}
	//debug
//	for (int i = 0; i < num_sum_ve; i++){
//		assert(v_edges_count[i] == adj_size[i]);
//	}
	printf("Load random graph %s, with %d (%d %d) vertices, %d edges\n",
			basename(filename), g_sum_ve, g_num_vertex[0],
			g_num_vertex[1], all_edges->size());
	delete[] adj_size;
	delete all_edges;

	return 1;
}

void buildComplementGraph(){
	/**
	 * 	int cmpg_num_edge;
		int* cmpg_cnt_edges;
		int** cmp_v_adj_vtx;
	 */
	cmpg_num_edge = g_num_vertex[0] * g_num_vertex[1] - g_num_edge;
	cmpg_cnt_edges = new int[g_sum_ve];
	cmpg_v_adj_vtx = new int*[g_sum_ve];

	int *mark = new int[g_sum_ve];
	for (int idx = 0; idx < g_num_vertex[0]; idx++){
		cmpg_cnt_edges[idx] = g_num_vertex[1] - g_v_edges_count[idx];
		cmpg_v_adj_vtx[idx] = new int[cmpg_cnt_edges[idx]];
		int cnt = 0;
		memset(mark + g_num_vertex[0], 0, sizeof(int) * g_num_vertex[1]);
		for (int i = 0; i < g_v_edges_count[idx]; i++){
			mark[g_v_adj_vertex[idx][i]] = 1;
		}
		for (int ridx = g_num_vertex[0]; ridx < g_sum_ve; ridx++){
			if (mark[ridx] == 0){
				cmpg_v_adj_vtx[idx][cnt++] = ridx;
			}
		}
		assert(cnt == cmpg_cnt_edges[idx]);
	}
	for (int idx = g_num_vertex[0]; idx < g_sum_ve; idx++){
		cmpg_cnt_edges[idx] = g_num_vertex[0] - g_v_edges_count[idx];
		cmpg_v_adj_vtx[idx] = new int[cmpg_cnt_edges[idx]];
		int cnt = 0;
		memset(mark, 0, sizeof(int) * g_num_vertex[0]);
		for (int i = 0; i < g_v_edges_count[idx];i++){
			mark[g_v_adj_vertex[idx][i]] = 1;
		}
		for (int ridx = 0; ridx < g_num_vertex[0]; ridx++){
			if (mark[ridx] == 0)
				cmpg_v_adj_vtx[idx][cnt++] = ridx;
		}
		assert(cnt == cmpg_cnt_edges[idx]);
	}
	delete[] mark;
}

static void print_vector(vector<int>::iterator it_start,
		vector<int>::iterator it_end){
	int real_size = it_end - it_start;
	int max_print_size = min(20, real_size);
	int n = 0;
	printf(" [%d] : ", real_size);
	while (it_start != it_end && n < max_print_size){
		printf("%d ", *it_start);
		it_start++;
	}
	if (n < max_print_size){
		printf("\n");
	}else {
		printf("...\n");
	}
}
static void print_vector(vector<int>& vec){
//	int max_print_size = min(20, (int)vec.size());
	int max_print_size = vec.size();
	printf(" [%d] : ", (int)vec.size());
	for (int i = 0; i < max_print_size; i++){
			printf("%d ", vec[i]);
	}
	if (max_print_size == (int)(vec.size()))
		printf("\n");
	else
		printf("...\n");
}


void greedySearch(vector<int> &greedySol){
	RandAccessList *lstU = ral_init(g_sum_ve);
	RandAccessList *lstV = ral_init(g_sum_ve);
	int *deg = new int[g_sum_ve];
	memcpy(deg, cmpg_cnt_edges, sizeof(int) * g_sum_ve);
	greedySol.clear();
	for (int i = 0; i < g_num_vertex[0]; i++)
		ral_add(lstU, i);
	for (int i = g_num_vertex[0]; i < g_sum_ve; i++)
		ral_add(lstV, i);

	for (int vtx = 0; vtx < g_sum_ve; vtx++){
		if (deg[vtx] == 0){
			greedySol.push_back(vtx);
			if (vtx < g_num_vertex[0])
				ral_delete(lstU, vtx);
			else
				ral_delete(lstV, vtx);
		}
	}
	int part = 0;
	int *weU = new int[g_num_vertex[0]];
	int *weV = new int[g_num_vertex[1]];
	while (1){
		//find the minim degre from part V
		int mindeg = MAX_VAL;
		if (part == 0){
			for (int i = 0; i < lstV->vnum; i++){
				int v = lstV->vlist[i];
				if (deg[v] < mindeg){
					mindeg = deg[v];
				}
			}
			//calculate the connection of each vertex in U
			memset(weU, 0, sizeof(int) * g_num_vertex[0]);
			for (int i = 0; i < lstV->vnum; i++){
				int v = lstV->vlist[i];
				if (mindeg == deg[v]){
					for (int i = 0; i < cmpg_cnt_edges[v]; i++){
						int u = cmpg_v_adj_vtx[v][i];
						weU[u]++;
					}
				}
			}
			//find the vertex in u with max weight in U
			int maxwe = 0;
			int u = -1;
			for (int i = 0; i < lstU->vnum; i++){
				int vtx = lstU->vlist[i];
				if (weU[vtx] > maxwe){
					maxwe = weU[vtx];
					u = vtx;
				}
			}
//			printf("Iter %lld: mindeg %d, remove %d \n", node_cnt, mindeg, u);
			assert(u >= 0 && u < g_num_vertex[0]);
			//remove u
			ral_delete(lstU, u);
			//update deg
			for (int i = 0; i < cmpg_cnt_edges[u]; i++){
				int v = cmpg_v_adj_vtx[u][i];
				deg[v]--;
				if (deg[v] == 0){
					greedySol.push_back(v);
					ral_delete(lstV, v);
				}
			}
		}else{
			//find a vertex with min deg in U
			for (int i = 0; i < lstU->vnum; i++){
				int u = lstU->vlist[i];
				if (deg[u] < mindeg){
					mindeg = deg[u];
				}
			}
			//calculate connection of each vertex in v
			memset(weV, 0, sizeof(int) * g_num_vertex[1]);
			for (int i = 0; i < lstU->vnum; i++){
				int u = lstU->vlist[i];
				if (mindeg == deg[u]){
					for (int i = 0; i < cmpg_cnt_edges[u]; i++){
						int v = cmpg_v_adj_vtx[u][i];
						weV[v - g_num_vertex[0]]++;
					}
				}
			}
			//find the vertex in v with max weight
			int maxwe = 0;
			int v = -1;
			for (int i = 0; i < lstV->vnum; i++){
				int vtx = lstV->vlist[i];
				if (weV[vtx - g_num_vertex[0]] > maxwe){
					maxwe = weV[vtx - g_num_vertex[0]];
					v = vtx;
				}
			}
			assert(v >= g_num_vertex[0] && v < g_sum_ve);
//			printf("Iter %d: mindeg %d, remove %d \n", node_cnt, mindeg, v);
			//remove u
			ral_delete(lstV, v);
			//update deg
			for (int i = 0; i < cmpg_cnt_edges[v]; i++){
				int u = cmpg_v_adj_vtx[v][i];
				deg[u]--;
				if (deg[u] == 0){
					greedySol.push_back(u);
					ral_delete(lstU, u);
				}
			}
		}
		if (lstV->vnum == 0 || lstU->vnum == 0)
			break;
		part = 1 - part;
		node_cnt++;
	}
	ral_release(lstU);
	ral_release(lstV);
	delete[] deg;
	delete[] weU;
	delete[] weV;
}


void read_parameters(int argc, char** argv) {
	for (int i = 1; i < argc; i += 2) {
		if (argv[i][0] != '-' || argv[i][2] != 0) {
			printf("arg format error\n");
			exit(0);
		} else if (argv[i][1] == 'f') {
			strncpy(para_file, argv[i + 1], 1000);
		} else if (argv[i][1] == 't') {
			para_max_seconds = atoi(argv[i + 1]);
		} else if (argv[i][1] == 's'){
			para_seed = atoi(argv[i+1]);
		}
	}
}

int verify_biclique(int half_size, vector<int>& clique){
	int *p1_vtx = new int[half_size];
	int np1 = 0;
	int *p2_vtx = new int[half_size];
	int np2 = 0;

	for (int i = 0; i < half_size * 2; i++){
		if (clique[i] < g_num_vertex[0]){
			p1_vtx[np1++] = clique[i];
		}else{
			p2_vtx[np2++] = clique[i];
		}
	}
	if (np1 != np2) return 0;
	for (int i = 0; i < np1; i++){
		int vp1 = p1_vtx[i];
		for (int j = 0; j < np2; j++){
			int vp2 = p2_vtx[j];
			int *pfind = find(g_v_adj_vertex[vp1],
					g_v_adj_vertex[vp1] + g_v_edges_count[vp1],
					vp2
					);
			if (pfind == g_v_adj_vertex[vp1] + g_v_edges_count[vp1]){
				return 0;
			}
		}
	}
	return 1;
}

void report_result(){
//	int test_half_size = 12;
//	int test_sol[24] = {2965, 3601, 10876, 17892, 21136, 29804, 29823, 34558, 35126, 36024, 61080, 68004, 94324, 94586, 96346, 96484, 96526, 96606, 97819, 98239, 98293, 98297, 98317, 98349,};
//	vector<int> test_best_sol(test_sol, test_sol+24);

	int rt_verify = verify_biclique(best_half_size, best_solution);
	if (rt_verify == 0){
		fprintf(stderr, "ERROR: The final solution is not feasible \n");
		return;
	}

	printf("$seed=%d\n", para_seed);
	printf("@solver=mwcp\n");
	printf("@para_file=%s\n", para_file);
	printf("@para_time=%d\n", para_max_seconds);
	printf("#orgvn1=%d\n",g_num_vertex[0]);
	printf("#orgvn2=%d\n",g_num_vertex[1]);
	printf("#orgnum=%d\n",g_sum_ve);
	printf("#orgenum=%d\n",g_num_edge);
	printf("#objv=%d\n", best_half_size);
	printf("#isopt=%d\n",is_optima);
	printf("#besttime=%.3f\n", (float)(best_time - start_time) / CLOCKS_PER_SEC);
	printf("#bestiter=%ld\n", best_node_id);
	printf("#totaltime=%.3f\n", (float)(total_time - start_time) /CLOCKS_PER_SEC);
	printf("#totaliter=%lld\n", node_cnt);
	//printf("#solutiontype = %d\n", stype);
	printf("Size of best solution %d \n", best_half_size);
	for (int i = 0; i < best_half_size; i++){
//		printf("%d ", Idx2Vertex(best_solution[i])+1);
		printf("%d ", best_solution[i]);
	}
	printf("\n");
	for (int i = best_half_size; i < best_solution.size(); i++){
//		printf("%d ", Idx2Vertex(best_solution[i])+1);
		printf("%d ", best_solution[i]);
	}
	printf("\n");
}

const char* file_suffix(char* filename){
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

void print_current_graph(){
	//print test
	cout << "Graph Size:" << g_sum_ve
			<< " up side:" << g_num_vertex[0]
			<< " down side:" << g_num_vertex[1] << endl;
	//omit the details if there are more than 50 vertices
	if (g_sum_ve > 50)
		return ;
	cout << "PARTIE 0 :" << endl;
	for (int vid = 0;  vid < g_num_vertex[0]; vid++){
		printf("Adjacent vertex of %d (size %d):",vid, g_v_edges_count[vid]);
		for (int j = 0; j < g_v_edges_count[vid]; j++){
			printf("%d ",g_v_adj_vertex[vid][j]) ;
		}
		cout << endl;
	}
	cout << "PARTIE 1 :" << endl;
	for (int vid = g_num_vertex[0];  vid < g_sum_ve; vid++){
		printf("Adjacent vertex of %d (size %d):", vid, g_v_edges_count[vid]);
		for (int j = 0; j < g_v_edges_count[vid]; j++){
			printf("%d ",g_v_adj_vertex[vid][j]) ;
		}
		cout << endl;
	}
}


int main(int argc, char** argv) {
	read_parameters(argc, argv);
	/*Load graph according to the suffix and prefix of filename*/
	int rt_load = 0;
	const char* fileext = file_suffix(para_file);
	if (0 == strcmp(fileext, "clq")){
		rt_load = load_rand_instance(para_file);
	}else{
		//prefix: the first 4 letters
		const char* bfname = basename(para_file);
		char fprefix[4];
		strncpy(fprefix, bfname, 4);
		if (0 == strncmp(fprefix, "out.", 4))
			rt_load = load_konect_instance(para_file);
	}
	if (rt_load != 1){
		fprintf(stderr, "failed in loading graph %s\n",para_file);
		exit(-1);
	}
	print_current_graph();
	buildComplementGraph();

	start_time = clock();
	total_time = 0;
	node_cnt = 0;

	vector<int> greedySol;
	greedySearch(greedySol);
	//reformat the best solution
	int usize = 0, vsize = 0;
	for (int i = 0; i < greedySol.size(); i++){
		if (greedySol[i] < g_num_vertex[0]) 	usize++;
		else	vsize++;
	}
	best_half_size = min(usize, vsize);
	best_solution.clear();
	int cnt = 0, i = 0;;
	while (cnt < best_half_size){
		if (greedySol[i] < g_num_vertex[0]){
			best_solution.push_back(greedySol[i]);
			cnt++;
		}
		i++;
	}
	cnt = 0;
	i = 0;
	while (cnt < best_half_size){
		if (greedySol[i] >= g_num_vertex[0]){
			best_solution.push_back(greedySol[i]);
			cnt++;
		}
		i++;
	}
	//set true initialy, if the algorithm stops before enumerate all the nodes, if will be reset as 0
	is_optima = 0;
	best_node_id = -1;
/**************************Precompute***********************************************/
	/***********OUTPUT Result*******************/
	total_time = clock();
	report_result();
	return 0;
}
