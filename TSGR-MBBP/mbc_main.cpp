//============================================================================
// Name        : mbbp_big.cpp
// Author      : JOEY
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "mbbp.hpp"

//using namespace std;

//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/actor-movie/out.actor-movie";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/brunson_revolution/out.brunson_revolution_revolution";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/dbpedia-producer/out.dbpedia-producer";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/real_life/livejournal-groupmemberships/out.livejournal-groupmemberships";
char para_file[1024] = "/home/zhou/benchmarks/MBBP/graphA/GraphU_250_0.05_1.clq";
//char para_file[1024] = "/home/zhou/benchmarks/MBBP/graphB/G_pbt_100_0.90_1.clq";
int para_optimum = 9999;
int para_seed = 1234567;
int para_max_seconds = 10000; //

/********config of heuristic**************/
int para_max_ls_iter = 10000;
/********config of exact algorithm**************/
float para_tt = 0.8;
int para_enum_max_vertex = 500;
int para_enum_max_seconds = 60;

/*************************Graph*************************************************/
typedef struct ST_Edge{
	int v1;
	int v2;
}Edge;
/**
 * To access the adjacent vertices of vertex i in upside
 * id of vertex i is : v_adj[0][i]
 * degree of vertex is: v_degree[i] v_degree[i] = v_adj_pos[0][i+1] - v_adj_pos[0][i]
 * for (i = adj_pos[0][i]; i < degree[i]; i++)
 * 		v_adj = v_adj_edges[0][i];
 */
static vector<Edge*>* orgs_edges;
static int org_v_num[PNUM_PARTIE];
static int org_e_num;
static int* org_v_edge_count;
static int** org_v_adj_vertex;

/*Current graph*/
int crt_num_vertex[PNUM_PARTIE];
int crt_num_sum_ve;
int* crt_v_org_id; // v_id[i] is the original id of vertex i
int* crt_v_edges_count; //v_edge[0][i] is vertex i's number of adjacent edges in the current grpah
int** crt_v_adj_vertex;	//v_adj_vertex[i][j] is vertex i's jth adjacet vertex
int crt_min_deg;

//static int cur_reduce_deg = 0;
#define Idx2Part(idx) ((idx)<crt_num_vertex[0]? 0 : 1)
#define Idx2Vertex(idx) ((idx)<crt_num_vertex[0]? (idx) : (idx)-crt_num_vertex[0])
#define PV2Idx(part,idx) ((part)*crt_num_vertex[0]+(idx));

//vector<int> vec_comp_size;
static int *crt_v_component_id = NULL;
static vector<vector<int>* > connect_components;
/*current sub graph*/
/*TODO: for graph with density less than 50%, we process it by the
 * complementary graph. */
int sg_comp_id = -1;
int sg_is_complement = 0;
int sg_num_vertex[PNUM_PARTIE];
int sg_num_sum_ve;
int *sg_v_org_id = NULL;	//this vector records the correpsonding id in the current graph(NOT the original graph)
int *sg_v_edge_count = NULL;
int **sg_v_adj_vertex = NULL;

static RT_LS *best_ever_rt = NULL;
int end_iter;
clock_t end_time;
int rt_is_opt = 0;

static const int Const_Sol_Max_Size = 2048;
static clock_t start_time;
//enum Sol_Type{REDUCE_OPT, BEST_EVER, BEST_KNOWN};
//static Sol_Type stype;

/*Load random generalized instances*/
int load_clq_instance(char* filename){
	ifstream infile(filename);
	char line[1024];
	char tmps1[1024];
	char tmps2[1024];

	if (infile == NULL){
		fprintf(stderr,"Can not find file %s\n", filename);
		return 0;
	}
	/*Read the comments*/
	infile.getline(line,  1024);
	while (line[0] != 'p')	infile.getline(line,1024);
	sscanf(line, "%s %s %d %d %d", tmps1, tmps2, &org_v_num[0],
			&org_v_num[1], &org_e_num);

	/*The current graph is equal to the original graph initially*/
	crt_num_vertex[0] = org_v_num[0];
	crt_num_vertex[1] = org_v_num[1];
	crt_num_sum_ve = crt_num_vertex[0] + crt_num_vertex[1];

	crt_v_org_id = new int[crt_num_sum_ve];
	for (int i = 0; i < crt_num_sum_ve; i++){
		crt_v_org_id[i] = i;
	}

	/*read edges*/
	orgs_edges = new vector<Edge*>();
	crt_v_edges_count = new int[crt_num_sum_ve];
	while (infile.getline(line, 1024)){
		int v1,v2;
		if (strlen(line) == 0)
			continue;
		if (line[0] != 'e')
			fprintf(stderr, "Line format error %s", line);
		sscanf(line, "%s %d %d", tmps1, &v1, &v2);
		v1--,v2--;
		Edge *pe = new Edge;
		pe->v1 = v1, pe->v2 = v2;
		orgs_edges->push_back(pe);

		crt_v_edges_count[v1]++;
		crt_v_edges_count[v2 + crt_num_vertex[0]]++;
	}
	assert(orgs_edges->size() == (int)org_e_num);

	crt_v_component_id = new int[crt_num_sum_ve];
	crt_v_adj_vertex = new int*[crt_num_sum_ve];
	crt_min_deg = MAX_VAL;
	vector<int> *init_part = new vector<int>(crt_num_sum_ve);
	for (int i = 0; i < crt_num_sum_ve; i++){
		init_part->at(i) = i;
		crt_v_component_id[i] = 0;
		crt_v_adj_vertex[i] = new int[crt_v_edges_count[i]];
		if (crt_v_edges_count[i] < crt_min_deg)
			crt_min_deg = crt_v_edges_count[i];
	}
	connect_components.push_back(init_part);

	int *adj_size = new int[crt_num_sum_ve];
	memset(adj_size, 0, sizeof(int) * crt_num_sum_ve);
	for (int i = 0; i < (int)orgs_edges->size() ;i++){
		Edge *pe = orgs_edges->at(i);
		int v1 = pe->v1, v2 = crt_num_vertex[0] + pe->v2;
		crt_v_adj_vertex[v1][adj_size[v1]++] = v2;
		crt_v_adj_vertex[v2][adj_size[v2]++] = v1;
	}
	//debug
//	for (int i = 0; i < num_sum_ve; i++){
//		assert(v_edges_count[i] == adj_size[i]);
//	}
	printf("Load clq graph %s, with %d (%d %d) vertices, %d edges\n",
			basename(filename), crt_num_sum_ve, crt_num_vertex[0], crt_num_vertex[1],
			orgs_edges->size());
	delete[] adj_size;

	/*reserve the orignial graph for verification*/
	org_v_edge_count = new int[crt_num_sum_ve];
	memcpy(org_v_edge_count, crt_v_edges_count, sizeof(int) * crt_num_sum_ve);
	org_v_adj_vertex = new int*[crt_num_sum_ve];
	for (int i = 0; i < org_v_num[0] + org_v_num[1]; i++){
		org_v_adj_vertex[i] = new int[org_v_edge_count[i]];
		memcpy(org_v_adj_vertex[i], crt_v_adj_vertex[i], sizeof(int) * org_v_edge_count[i]);
	}
	return 1;
}
int load_konect_instance(char* filename){
	ifstream infile(filename);
	if (infile == NULL){
		fprintf(stderr,"Can not open file %s\n", filename);
		return 0;
	}
	char line[1000];
	int vid1, vid2;

	/*ignore all the comment lines*/

	org_v_num[0] = 0;
	org_v_num[1] = 0;

	orgs_edges = new vector<Edge*>();
	/*Read the original graph, record the maximum id of vertices
	 * as the maximum number of the vertices*/
	while (infile.getline(line,  1000)){
		//cout << id1 << " " << id2 << endl;
		if (line[0] == '%')
			continue;
		else
			sscanf(line, "%d %d", &vid1, &vid2);
		if (vid1 > org_v_num[0])
			org_v_num[0] = vid1;
		if (vid2 > org_v_num[1]){
			org_v_num[1] = vid2;
		}
		Edge *pe = new Edge;
		pe->v1 = vid1;
		pe->v2 = vid2;
		orgs_edges->push_back(pe);
	}
	printf("Vertices %d %d, Edges: %d\n", org_v_num[0], org_v_num[1], (int)orgs_edges->size());

	/*Allocate space for  current graph*/
	crt_num_vertex[0] = org_v_num[0];
	crt_num_vertex[1] = org_v_num[1];
	crt_num_sum_ve = crt_num_vertex[0] + crt_num_vertex[1];
	crt_v_org_id = new int[crt_num_sum_ve];

	for (int vi = 0; vi < crt_num_sum_ve; vi++){
		crt_v_org_id[vi] = Idx2Vertex(vi)+1;
	}
	/*estimate adjacent edge numbers for each nodes*/
	int *estimate_edge_cnt = new int[crt_num_sum_ve];
	memset(estimate_edge_cnt, 0 , sizeof(int) * crt_num_sum_ve);
	/*count the number of adjacent edges for each vertex*/
	for (vector<Edge*>::iterator itr = orgs_edges->begin();
			itr != orgs_edges->end();
			itr++){
		vid1 = ((Edge*)(*itr))->v1 - 1;
		vid2 = ((Edge*)(*itr))->v2 - 1;
		//cout << id1 <<" "<< id2 << endl;
		int idx1 = PV2Idx(0, vid1);
		int idx2 = PV2Idx(1, vid2);
		estimate_edge_cnt[idx1]++;
		estimate_edge_cnt[idx2]++;
	}

	/*allocate space for the adjacent vertex list of each vertex*/
	crt_v_edges_count = new int[crt_num_sum_ve];
	crt_v_adj_vertex = new int*[crt_num_sum_ve];
	for (int idx = 0; idx < crt_num_sum_ve; idx++){
		if (estimate_edge_cnt[idx] > 0)
			crt_v_adj_vertex[idx] = new int[estimate_edge_cnt[idx]];
		else
			crt_v_adj_vertex[idx] = (int*)0; //for secure
//		if (estimate_edge_cnt[idx] < crt_min_deg)
//			crt_min_deg = crt_v_edges_count[idx];
	}

	/*Build up the adjacent vertex list*/
	memset(crt_v_edges_count, 0, sizeof(int) * crt_num_sum_ve);
	for (vector<Edge*>::iterator itr = orgs_edges->begin();
				itr != orgs_edges->end();
				itr++){
			vid1 = ((Edge*)(*itr))->v1 - 1;
			vid2 = ((Edge*)(*itr))->v2 - 1;
			//cout << id1 <<" "<< id2 << endl;
			int idx1 = PV2Idx(0, vid1);
			int idx2 = PV2Idx(1, vid2);
			int* rt_find1 = find(crt_v_adj_vertex[idx1],
					crt_v_adj_vertex[idx1] + crt_v_edges_count[idx1], idx2);
			if (rt_find1 == crt_v_adj_vertex[idx1] + crt_v_edges_count[idx1]){
				crt_v_adj_vertex[idx1][crt_v_edges_count[idx1]++] = idx2;
			}
			int* rt_find2 = find(crt_v_adj_vertex[idx2],
					crt_v_adj_vertex[idx2]+crt_v_edges_count[idx2], idx1);
			if (rt_find2 == crt_v_adj_vertex[idx2]+crt_v_edges_count[idx2]){
				crt_v_adj_vertex[idx2][crt_v_edges_count[idx2]++] = idx1;
			}
		}
	delete[] estimate_edge_cnt;

	crt_min_deg = MAX_VAL;
	crt_v_component_id = new int[crt_num_sum_ve];
	vector<int> *vec_init_part = new vector<int>();
	for (int v = 0; v < crt_num_sum_ve; v++){
		/*Suppose all the vertices are in the same connected graph*/
		vec_init_part->push_back(v);
		crt_v_component_id[v] = 0;
		/*update min degree*/
		if (crt_v_edges_count[v] <= crt_min_deg){
			crt_min_deg = crt_v_edges_count[v];
		}
	}
	connect_components.push_back(vec_init_part);

	/*reserve the orignial graph for verification*/
	org_e_num = orgs_edges->size();
	org_v_edge_count = new int[crt_num_sum_ve];
	memcpy(org_v_edge_count, crt_v_edges_count, sizeof(int) * crt_num_sum_ve);
	org_v_adj_vertex = new int*[crt_num_sum_ve];
	for (int i = 0; i < org_v_num[0] + org_v_num[1]; i++){
		org_v_adj_vertex[i] = new int[org_v_edge_count[i]];
		memcpy(org_v_adj_vertex[i], crt_v_adj_vertex[i], sizeof(int) * org_v_edge_count[i]);
	}
	return 1;
}


void print_current_graph(){
	//print test
	cout << "Graph Size:" << crt_num_sum_ve
			<< " up side:" << crt_num_vertex[0]
			<< " down side:" << crt_num_vertex[1] << endl;
	if (crt_num_sum_ve > 50)
		return ;
	cout << "PARTIE 0 :" << endl;
	for (int id = 0;  id < crt_num_vertex[0]; id++){
		int pi = Idx2Part(id);
		int vi = Idx2Vertex(id);
		printf("Adjacent vertex of %d(%d) (size %d):",id, vi, crt_v_edges_count[id]);
		for (int j = 0; j < crt_v_edges_count[id]; j++){
			printf("%d(%d) ",crt_v_adj_vertex[id][j], Idx2Vertex(crt_v_adj_vertex[id][j])) ;
		}
		cout << endl;
	}
	cout << "PARTIE 1 :" << endl;
	for (int id = crt_num_vertex[0];  id < crt_num_sum_ve; id++){
		int pi = Idx2Part(id);
		int vi = Idx2Vertex(id);
		printf("Adjacent vertex of %d(%d) (size %d):", id, vi, crt_v_edges_count[id]);
		for (int j = 0; j < crt_v_edges_count[id]; j++){
			printf("%d(%d) ",crt_v_adj_vertex[id][j], Idx2Vertex(crt_v_adj_vertex[id][j])) ;
		}
		cout << endl;
	}
}

void rebuild_rest_graph(vector<int> &rm_vertices){
	int part_rm_cnt[2] = {0,0};
	/*Mark the vertices to be removed*/
	int *rmflag = new int[crt_num_sum_ve];
	/*The left adjacent edges after reduce*/
	//int *reduced_edge_count = new int[num_sum_ve];
	//memcpy(reduced_edge_count, v_edges_count, sizeof(int) * num_sum_ve);
	memset(rmflag, 0, sizeof(int) * crt_num_sum_ve);
	for (int i = 0; i < rm_vertices.size(); i++){
		int v = rm_vertices[i];
		rmflag[v] = 1;
		if (v < crt_num_vertex[0])
			part_rm_cnt[0]++;
		else
			part_rm_cnt[1]++;
	}
	int num_rest_ve = crt_num_sum_ve - rm_vertices.size();
	int *new_id = new int[crt_num_sum_ve]; // resin
	int *org_id = new int[num_rest_ve];
	int count = 0;	//count the rest vertices
	//resign id to the rest vertices
	for (int idx = 0; idx < crt_num_sum_ve; idx++){
		if (!rmflag[idx]){
			new_id[idx] = count;
			org_id[count] = crt_v_org_id[idx]; //the id in original graph, not the current graph
			count++;
		}
	}
	int n_edges= 0;
	int **new_adj_tbl = new int*[num_rest_ve];
	int *new_edge_count = new int[num_rest_ve];
	crt_min_deg = MAX_VAL;	//recaclculate minimum degree
	for (int idx_prev = 0; idx_prev < crt_num_sum_ve; idx_prev++){
		if (rmflag[idx_prev]){
			delete[] crt_v_adj_vertex[idx_prev];
			continue;
		}
		int idx_new = new_id[idx_prev];
		vector<int> rest_edges;
		for (int i = 0; i < crt_v_edges_count[idx_prev]; i++){
			int vi_adj = crt_v_adj_vertex[idx_prev][i];
			if (!rmflag[vi_adj]){
				rest_edges.push_back(new_id[vi_adj]);
				n_edges++;
			}
		}
		new_edge_count[idx_new] = rest_edges.size();
		new_adj_tbl[idx_new] = new int[rest_edges.size()];
		copy(rest_edges.begin(), rest_edges.end(), new_adj_tbl[idx_new]);
		if (new_edge_count[idx_new] < crt_min_deg){
			crt_min_deg = new_edge_count[idx_new];
		}
//		for (int i = 0; i < rest_edges.size(); i++){
//			new_adj_tbl[idx_new][i] = rest_edges[i];
//		}
		delete[] crt_v_adj_vertex[idx_prev];
	}

	/*assign to the new graph*/
	delete[] crt_v_adj_vertex;
	crt_v_adj_vertex = new_adj_tbl;
	delete[] crt_v_edges_count;
	crt_v_edges_count = new_edge_count;

	crt_num_sum_ve = num_rest_ve;
	crt_num_vertex[0] = crt_num_vertex[0] - part_rm_cnt[0];
	crt_num_vertex[1] = crt_num_vertex[1] - part_rm_cnt[1];

	delete[] crt_v_org_id;
	crt_v_org_id = org_id;

	delete[] new_id;
	delete[] rmflag;

}

/**
 * remove all the vertex with degree equal to or less than reduce_deg
 */
	/*The left adjacent edges after reduce*/
void find_peel_vertices(int reduce_deg, vector<int>& vtx_2_peel){
	int *reduced_edge_count = new int[crt_num_sum_ve];
	int *mark = new int[crt_num_sum_ve];
	memset(mark, 0, sizeof(int) * crt_num_sum_ve);
	memcpy(reduced_edge_count, crt_v_edges_count, sizeof(int) * crt_num_sum_ve);
	queue<int> rm_que;
	if (reduce_deg >= crt_min_deg){
		for (int idx = 0; idx < crt_num_sum_ve; idx++){
			if (reduced_edge_count[idx] <= reduce_deg){
				rm_que.push(idx);
				mark[idx] = 1;
			}
		}
		while (!rm_que.empty()){
			int idx = rm_que.front();
			rm_que.pop();
			vtx_2_peel.push_back(idx);
			for (int i = 0; i < crt_v_edges_count[idx]; i++){
				int adjv = crt_v_adj_vertex[idx][i];
				reduced_edge_count[adjv]--;
				if (!mark[adjv] && reduced_edge_count[adjv] <= reduce_deg){
					mark[adjv] = 1;
					rm_que.push(adjv);
				}
			}
		}
	}
	delete[] reduced_edge_count;
	delete[] mark;
	//make all the vertices with degree less than lb by BFS search
	//debug
//	printf("remove %d %d \n", part_rm_cnt[0], part_rm_cnt[1]);
//	for (int i = 0; i < num_sum_ve; i++){
//		printf("%d:%d  ", Idx2Vertex(i), rmflag[i]);
//		if (i == num_vertex[0])
//			printf("\n");
//	}
//	printf("\n");
}



/**
 *	find all the connect components, mark them in v_component_id
 */
void partition_connect_commponent(){
	queue<int> q;
	int comp_id = 0;

	if (crt_v_component_id == NULL){
		crt_v_component_id = new int[crt_num_sum_ve];
	}
	for (int idx = 0; idx < crt_num_sum_ve; idx++)
		crt_v_component_id[idx] = -1;
	/*remove all the exisiting paritions*/
	for (int i = 0; i < (int)connect_components.size(); i++){
		delete connect_components[i];
	}
	connect_components.clear();

	for (int vi = 0; vi < crt_num_sum_ve; vi++){
		if (crt_v_component_id[vi] == -1){
			vector<int>* vec_cur_comp = new vector<int>();
			vec_cur_comp->push_back(vi);
			crt_v_component_id[vi] = comp_id;
			q.push(vi);
			while (!q.empty()){
				int v_cur = q.front();
				q.pop();
				for (int i = 0; i < crt_v_edges_count[v_cur]; i++){
					int vadj = crt_v_adj_vertex[v_cur][i];
					if (crt_v_component_id[vadj] == -1){
						q.push(vadj);
						vec_cur_comp->push_back(vadj);
						crt_v_component_id[vadj] = comp_id;
					}
				}
			}
			comp_id++;
			connect_components.push_back(vec_cur_comp);
			printf("component %d size %d \n", comp_id, vec_cur_comp->size());
		}
	}
}

/**
 * Construct the current subgraph from the connect component comp_id
 * the graph is mantained by
 *	int sg_is_complement = 0;
 *	int sg_num_vertex[PNUM_PARTIE];
 *	int *sg_v_org_id;
 *  int *sg_v_edge_count;
 *	int **sg_v_adj_vertex;
 */
void reconstruct_subgraph(int comp_id){
	int *new_id = new int[crt_num_sum_ve];

	assert(comp_id < (int)connect_components.size());
	/*reallocate memeory*/
	if (sg_v_edge_count != NULL){
		delete[] sg_v_org_id;
		delete[] sg_v_edge_count;
		for (int i = 0; i < sg_num_sum_ve; i++)
			delete[] sg_v_adj_vertex[i];
		delete[] sg_v_adj_vertex;
		sg_v_adj_vertex = NULL;
		sg_v_edge_count = NULL;
	}
	sg_comp_id = comp_id;
	sg_is_complement = 0;
	sg_num_sum_ve = connect_components[comp_id]->size();
	sg_num_vertex[0] = sg_num_vertex[1] = 0;
	sg_v_org_id = new int[sg_num_sum_ve];
	sg_v_edge_count = new int[sg_num_sum_ve];
	sg_v_adj_vertex = new int*[sg_num_sum_ve];

	/*renumber the vertices, rearrange the vertices by their ids so that
	 * the vertices of the up side locate before that of the down side*/
	sort(connect_components[comp_id]->begin(), connect_components[comp_id]->end());
	for (int idx = 0; idx < (int)(connect_components[comp_id]->size()); idx++){
		int org_vtx = connect_components[comp_id]->at(idx);
		new_id[org_vtx] = idx;
	}

	for (int idx = 0; idx < (int)(connect_components[comp_id]->size()); idx++){
		int org_vtx = connect_components[comp_id]->at(idx);
		sg_num_vertex[Idx2Part(org_vtx)]++;
		sg_v_org_id[idx] = org_vtx;
		sg_v_edge_count[idx] = 0;
		vector<int> vec_tmp_adj_vtx;
		for (int adj_idx = 0; adj_idx < crt_v_edges_count[org_vtx]; adj_idx++){
			int org_adj_vtx = crt_v_adj_vertex[org_vtx][adj_idx];
			/*The adjacent vertices must be in the same sugraph of vidx, but
			 * we keep such judge in case of new partition methods*/
			if (crt_v_component_id[org_adj_vtx] == comp_id){
				sg_v_edge_count[idx]++;
				vec_tmp_adj_vtx.push_back(new_id[org_adj_vtx]);
			}
		}
		sg_v_adj_vertex[idx] = new int[sg_v_edge_count[idx]];
		//copy
		copy(vec_tmp_adj_vtx.begin(), vec_tmp_adj_vtx.end(), sg_v_adj_vertex[idx]);
//		for (int i = 0; i < sg_v_edge_count[idx]; i++)
//			sg_v_adj_vertex[idx][i] = vec_tmp_adj_vtx[i];
	}
	delete[] new_id;
}

void print_subgraph(){
	printf("sub_graph %d, total vtx: %d, up %d down %d \n",
			sg_comp_id, sg_num_sum_ve, sg_num_vertex[0], sg_num_vertex[1]);
	printf("is complete %d\n", sg_is_complement);
	if (crt_num_sum_ve > 50)
		return;
	for (int idx = 0; idx < crt_num_sum_ve; idx++){
		printf("Vertex %d (%d) part %d deg (%d): ", idx, sg_v_org_id[idx], Idx2Part(sg_v_org_id[idx]),
				sg_v_edge_count[idx]);
		for (int adj_idx = 0; adj_idx < sg_v_edge_count[idx]; adj_idx++){
			printf(" %d(%d) ", crt_v_adj_vertex[idx][adj_idx],
					sg_v_org_id[crt_v_adj_vertex[idx][adj_idx]]);
		}
		printf("\n");
	}
}

void initial_data(){
	srand(para_seed);
	start_time = clock();
	end_time = 0;
	end_iter = 0;
	rt_is_opt = 0;
	best_ever_rt = new RT_LS;
	best_ever_rt->best_clique = new int[Const_Sol_Max_Size];
	best_ever_rt->best_balanced_size = 0;
	best_ever_rt->best_up_size = 0;
	best_ever_rt->best_down_size = 0;
	best_ever_rt->best_end_time = 0;
}
void copy_RT_LS(RT_LS *dest, RT_LS *src){
	assert(dest->best_clique != NULL);
	memcpy(dest->best_clique, src->best_clique,
			sizeof(int) * (src->best_down_size + src->best_up_size));
	dest->best_down_size = src->best_down_size;
	dest->best_up_size = src->best_up_size;
	dest->best_balanced_size = src->best_balanced_size;
	dest->best_iter = src->best_iter;
	//dest->start_time = src->start_time;
	//dest->total_iter = src->total_iter;
	dest->best_end_time = src->best_end_time;
}

void exact_search(){
	int total_sec;
	int total_iter;
	int is_opt = 0;
	RT_LS *rt_ls = new RT_LS;
	rt_ls->best_clique = new int[Const_Sol_Max_Size];

	initial_data();
	reconstruct_subgraph(0);
	enum_search(0, para_max_seconds, rt_ls, &total_iter, &total_sec, &is_opt);
	copy_RT_LS(best_ever_rt, rt_ls);
	rt_is_opt = is_opt;

	for (int i = 0; i < (int)connect_components.size(); i++)
		delete connect_components[i];
	delete[] rt_ls->best_clique;
	delete rt_ls;
	end_time = clock();
	end_iter = total_iter;
	return;
}


void reduce_and_search() {
	vector<int> vec_sol;
	int end = 0;
	RT_LS *rt_ls = new RT_LS;
	rt_ls->best_clique = new int[Const_Sol_Max_Size];

	initial_data();
	while (!end){
		int ls_iter;
		int ls_sec;
		//local search
		feasible_local_search(para_max_ls_iter, para_max_seconds, rt_ls, &ls_iter, &ls_sec);
		end_iter += ls_iter;
//		printf("A new solution %d after %d iter\n",rt_ls->best_balanced_size, ls_iter);
		if (rt_ls->best_balanced_size > best_ever_rt->best_balanced_size){
			copy_RT_LS(best_ever_rt, rt_ls);
		}

		/*peel procedure*/
		while (best_ever_rt->best_balanced_size >= crt_min_deg){
			vector<int> vtx_2_remove;
			//peel graph
			find_peel_vertices(best_ever_rt->best_balanced_size, vtx_2_remove);
			printf("Peel %d vertices, left %d \n", vtx_2_remove.size(), crt_num_sum_ve - vtx_2_remove.size());
			if (vtx_2_remove.size() == 0)
				break;
			rebuild_rest_graph(vtx_2_remove);
			//partition graph to connect subgaphs
			partition_connect_commponent();
			printf("Partition into %d subgraph\n", connect_components.size());
			vtx_2_remove.clear();
			for (int com_id = 0; com_id < (int)(connect_components.size()); com_id++){
				if ((int)(connect_components[com_id]->size()) <= para_enum_max_vertex){
					reconstruct_subgraph(com_id);
					int enum_sec;
					int enum_iter;
					int is_opt;
					//exhaust search
					enum_search(rt_ls->best_balanced_size, 10, rt_ls, &enum_iter, &enum_sec, &is_opt);
					end_iter += enum_iter;
					if (rt_ls->best_balanced_size > best_ever_rt->best_balanced_size){
						printf("Eact find a better solution %d in sg %d\n",
								rt_ls->best_balanced_size, com_id);
						copy_RT_LS(best_ever_rt, rt_ls);
					}
					if (is_opt){
					//remove these vertices
						vtx_2_remove.insert(vtx_2_remove.end(), connect_components[com_id]->begin(),
												connect_components[com_id]->end());
					}
				}
				if ((float)(clock() - start_time) / CLOCKS_PER_SEC >= para_max_seconds){
					goto stop;
				}
			}
			if (vtx_2_remove.size() > 0)
				printf("Prune %d vertices by enumerate \n", vtx_2_remove.size());
				rebuild_rest_graph(vtx_2_remove);
		}
		if (crt_num_sum_ve <= best_ever_rt->best_balanced_size * 2){
			rt_is_opt = 1;
			end = 1;
		}
		if ((float)(clock() - start_time) / CLOCKS_PER_SEC > para_max_seconds){
			end = 1;
		}
	}
stop:
	for (int i = 0; i < (int)connect_components.size(); i++)
		delete connect_components[i];
	delete[] rt_ls->best_clique;
	delete rt_ls;
	end_time = clock();
	return;
}

void show_usage(){
	printf("-f <filename> -s <seed> -o <best know result> -t <maximum seconds>\n");
}

void read_parameters(int argc ,char** argv){
	for (int i = 1; i < argc; i+=2){
		if (argv[i][0] != '-' || argv[i][2] != 0){
			show_usage();
			exit(0);
		}else if (argv[i][1] == 'f'){
			strncpy(para_file, argv[i+1],1000);
		}else if (argv[i][1] == 's'){
			para_seed = atoi(argv[i+1]);
		}else if (argv[i][1] == 't'){
			para_max_seconds = atoi(argv[i+1]);
		}else if (argv[i][1] == 'o'){
			para_optimum = atoi(argv[i+1]);
		}else if (argv[i][1] == 'l'){
			para_max_ls_iter = atoi(argv[i+1]);
		}else if (argv[i][1] == 'a'){
			para_tt = atof(argv[i+1]);
		}
	}
}


void report_result(){
	/** $seed=158462
		@solver=splex
		@param1=
		@para2=
		#iter= 123545
		#objv=45
		#bestiter= 1235
	 */

	/*format print */
	printf("$seed=%d\n", para_seed);
	printf("@solver=mwcp\n");
	printf("@para_file=%s\n", para_file);
	printf("@para_time=%d\n", para_max_seconds);
	printf("@para_T=%.2f\n", para_tt);
	printf("#orgvn1=%d\n",org_v_num[0]);
	printf("#orgvn2=%d\n",org_v_num[1]);
	printf("#orgnum=%d\n",org_v_num[0]+org_v_num[1]);
	printf("#orgenum=%d\n",(int)orgs_edges->size());
	printf("#redvnum=%d\n",crt_num_sum_ve);
	printf("#redvn1=%d\n",crt_num_vertex[0]);
	printf("#redvn2=%d\n",crt_num_vertex[1]);
	printf("#objv=%d\n", best_ever_rt->best_balanced_size);
	printf("#isopt=%d\n",rt_is_opt);
	printf("#besttime=%.3f\n", (float)(best_ever_rt->best_end_time - start_time) / CLOCKS_PER_SEC);
	printf("#bestiter=%d\n", best_ever_rt->best_iter);
	printf("#totaltime=%.3f\n", (float)(end_time - start_time) /CLOCKS_PER_SEC);
	printf("#totaliter=%d\n", end_iter);
	//printf("#solutiontype = %d\n", stype);
	sort(best_ever_rt->best_clique, best_ever_rt->best_clique +
			best_ever_rt->best_up_size + best_ever_rt->best_down_size);
	for (int i = 0; i < best_ever_rt->best_up_size + best_ever_rt->best_down_size; i++){
		printf("%d ",best_ever_rt->best_clique[i]);
	}
	printf("\n");
}

const char* file_suffix(char* filename){
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}
int main(int argc, char **argv) {
	int rt_load = 0;

	read_parameters(argc, argv);

	const char* fileext = file_suffix(para_file);
	if (0 == strcmp(fileext, "clq")){
		rt_load = load_clq_instance(para_file);
	}else{
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
	//debug
	//print_current_graph();
	/*The reduce and search algorithm*/
	reduce_and_search();

	//exact_search();

	report_result();

	return 0;
}
