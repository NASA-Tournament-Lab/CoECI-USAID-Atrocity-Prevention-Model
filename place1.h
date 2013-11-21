/*
 * Copyright 2012-2013 TopCoder, Inc.

 *
 * This code was developed under U.S. government contract NNH10CD71C. 

 *
 * Licensed under the Apache License, Version 2.0 (the "License");

 * You may not use this file except in compliance with the License.

 * You may obtain a copy of the License at:
 *     http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software

 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

 * See the License for the specific language governing permissions and

 * limitations under the License.
 */

#include "base.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <iterator>
#include <algorithm>
#include <queue>
#include <functional>
#include <sstream>
#include <complex>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <sys/time.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <string.h>

#ifdef AST
#include <assert.h>
#endif

using namespace std;

namespace place1 {

#define FOR(i,a,b) for(int i=(a);i<(b);++i)
#define REP(i,a) FOR(i,0,a)
#define ZERO(m) memset(m,0,sizeof(m));
#define DB(a) cerr << #a << ":" << (a) << endl;
#define SIZE(a) ((int) a.size())

#ifdef LOCAL
	ofstream out_info("out_info.txt");
#else
	ostream& out_info = cerr;
#endif

double time_begin;
double cur_time;
int calledgettime;

void starttime()
{
	unsigned long long time;
	__asm__ volatile ("rdtsc" : "=A" (time));
#ifdef LOCAL
	time_begin = time / 3.595e9; 
#else
	time_begin = time / 3.595e9;
#endif
}

double gettime()
{
	++calledgettime;
	unsigned long long time;
	__asm__ volatile ("rdtsc" : "=A" (time));
#ifdef LOCAL
	cur_time = time / 3.595e9 - time_begin; 
	return cur_time;
#else
	cur_time = time / 3.595e9 - time_begin;
	return cur_time;
#endif
}

#define REG_SIZE 3671
#define COU_SIZE 254
#define FIRST_DAY 11284
#define LAST_DAY 17644
#define WHOLE_PERIOD (LAST_DAY - FIRST_DAY + 31)

vector<int> inner_regs[COU_SIZE];
int in_country[REG_SIZE];
double reg_longi[REG_SIZE];
double reg_lati[REG_SIZE];


int glb_reg_atro[REG_SIZE][WHOLE_PERIOD];
int glb_cou_atro[COU_SIZE][WHOLE_PERIOD];
vector<int> glb_reg_history[REG_SIZE];
vector<int> glb_cou_history[COU_SIZE];

int rct_reg_atro[REG_SIZE]; //correct answer at day now - 30

#define REGION_RATIO_SAMPLE 1000
#define REG_CHECK_DAY_LEN 18
int reg_check_day[REG_CHECK_DAY_LEN + 1] = {3, 7, 14, 21, 28, 42, 63, 84, 168, 252, 364, 532, 728, 992, 1456, 2184, 2912, 10000, REGION_RATIO_SAMPLE};
#define COU_CHECK_DAY_LEN 18
int cou_check_day[COU_CHECK_DAY_LEN + 2] = {3, 7, 14, 21, 28, 42, 63, 84, 168, 252, 364, 532, 728, 992, 1456, 2184, 2912, 10000, REGION_RATIO_SAMPLE, 35};

int fix_day_len[2] = {30, 0};
int reg_past_atro[2][REG_SIZE][REG_CHECK_DAY_LEN + 1];
int cou_past_atro[2][COU_SIZE + 1][COU_CHECK_DAY_LEN + 2];

#define EVENT_TYPE_CNT 31
#define SOC_BUF_SIZE 450
// bigger than max check day + 30.
#define REG_EVENT_CHECK_LEN 2
int reg_eve_check_day[REG_EVENT_CHECK_LEN] = {63, 364};//don't forget to update SOC_BUF_SIZE
#define COU_EVENT_CHECK_LEN 2
int cou_eve_check_day[COU_EVENT_CHECK_LEN] = {63, 364};
int reg_past_soc[2][REG_SIZE][REG_EVENT_CHECK_LEN][EVENT_TYPE_CNT];
int cou_past_soc[2][COU_SIZE + 1][COU_EVENT_CHECK_LEN][EVENT_TYPE_CNT];
int soc_buf_ptr = 0;
int reg_soc_buf[REG_SIZE][SOC_BUF_SIZE][EVENT_TYPE_CNT];
int cou_soc_buf[COU_SIZE + 1][SOC_BUF_SIZE][EVENT_TYPE_CNT];
#define START_COLLECT_SOC_STA 14000

//soc sta only:
#define SOC_CANDI_CNT (4 * EVENT_TYPE_CNT)
double soc_tot_r1[SOC_CANDI_CNT][2];
double soc_tot_r2[SOC_CANDI_CNT][2];
int soc_sta_def_times[SOC_CANDI_CNT][2];
int soc_sta_undef_times[SOC_CANDI_CNT][2];

#define THRE_LEN 8
double threshold[THRE_LEN] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 100};
//double thre_ratio[THRE_LEN] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double thre_ratio[THRE_LEN] = {1.0, 1.0, 1.0, 1.0, 1.1, 1.2, 1.2, 1.2};

#define SCATER_FEA_CNT 7
#define SOC_FEA_CNT 31
#define DT_LIB_SIZE (REG_SIZE * 270)
#define DT_PRO_SIZE (REG_CHECK_DAY_LEN + COU_CHECK_DAY_LEN + SCATER_FEA_CNT + SOC_FEA_CNT)
double dt_lib_pro[DT_LIB_SIZE][DT_PRO_SIZE];
double dt_prd_pro[REG_SIZE][DT_PRO_SIZE];
int dt_lib_ans[DT_LIB_SIZE];
int dt_lib_pos = 0;

int prop_weight[DT_PRO_SIZE]
={0,286,0,405,0,162,0,286,286,0,286,0,286,102,286,141,286,436,
219,401,0,0,401,194,0,107,401,273,0,0,401,0,0,0,274,0,
0,0,0,396,396,487,315,
221,0,0,48,40,0,0,0,143,144,0,0,48,99,315,0,349,0,315,0,0,315,0,0,272,0,0,0,0,174,0};

#define CONFIG_COUNT 4
double config_v[CONFIG_COUNT] ={209260,168,0.0015,0.3};
 
// dt sta
int used_pos_cnt[DT_LIB_SIZE];
int pick_pro_cnt[DT_PRO_SIZE];
int best_pro_cnt[DT_PRO_SIZE];

unsigned long long now_rand = 1;

void set_rand_seed(unsigned long long seed)
{
	now_rand = seed;
}

unsigned long long math_rand()
{
	now_rand = ((now_rand * 6364136223846793005ULL + 1442695040888963407ULL) >> 1);
	return now_rand;
}

int get_conf_id(double v)
{
	int conf_id = 0;
	for (int k = 0; k + 1 < THRE_LEN; ++k)
	{
		if (v > threshold[k])
		{
			++conf_id;
		}
	}
	return conf_id;
}

struct Node
{
	int type;
	int pid;
	double vle; // vle for leaf, split vle for mid.
	int left_id;
	int right_id;
	vector<int> temp_dids;
};

struct DT
{
	int construct_day;
	vector<Node> nodes;
	double zero_log2(double vle)
    {
        if (vle < 1e-8)
        {
			return 0;
		}
		return log2(vle);
    }
    double get_entropy(int a0, int a1, int b0, int b1)
    {
        double an = a0 + a1;
        double bn = b0 + b1;
        double ent = - a0 * zero_log2(a0 / an) 
				- a1 * zero_log2(a1 / an)
				- b0 * zero_log2(b0 / bn)
				- b1 * zero_log2(b1 / bn);
        return ent;
    }
	double get_gini(int a0, int a1, int b0, int b1)
	{
		double an = a0 + a1;
        double bn = b0 + b1;
		double ga = 1 - (a0 / an) * (a0 / an) - (a1 / an) * (a1 / an);
		double gb = 1 - (b0 / bn) * (b0 / bn) - (b1 / bn) * (b1 / bn);
		return ga * an + gb * bn;
	}
	double get_conf_score(int a0, int a1, int b0, int b1)
	{
		double an = a0 + a1;
        double bn = b0 + b1;
		double pa = a1 / an;
		double pb = b1 / bn;
		double scorea = a1 * (pa - pa * pa / 2) - a0 * pa * pa / 2;
		double scoreb = b1 * (pb - pb * pb / 2) - b0 * pb * pb / 2;
		return 1e+7 - (scorea + scoreb);
	}
	double get_avg(vector<double>& svles)
	{
		int sz = SIZE(svles);
		double sum = 0;
		REP (i, sz)
		{
			sum += svles[i];
		}
		return sum / sz;
	}
	double get_value(vector<double>& prop)
	{
		int now_node = 0;
		while (nodes[now_node].type == 1)
		{
			double pv = prop[nodes[now_node].pid];
			double cv = nodes[now_node].vle;
			if (pv < cv)
			{
				now_node = nodes[now_node].left_id;
			}
			else
			{
				now_node = nodes[now_node].right_id;
			}
		}
		return nodes[now_node].vle;
	}
	void build_dt(int day_seed)
	{
		set_rand_seed(day_seed);
		nodes.reserve(2000); //ensure not extend space.
		nodes.push_back(Node());
		Node& node0 = nodes[0];
		int data_size = (int) config_v[0];
		node0.temp_dids.reserve(data_size);
		REP (i, data_size)
		{
			int pos = math_rand() % DT_LIB_SIZE;
			used_pos_cnt[pos]++;
			node0.temp_dids.push_back(pos);
		}
		int head = 0;
		vector<int> prop_pot;
		prop_pot.reserve(20000);
		REP (i, DT_PRO_SIZE)
		{
			REP (j, prop_weight[i])
			{
				prop_pot.push_back(i);
			}
		}
		while(head < nodes.size())
		{
			Node& node = nodes[head];
			if (SIZE(node.temp_dids) == 0)
			{
				cerr << "WARNING DT CONSTRUCT" << endl;
			}
			if (SIZE(node.temp_dids) < (int) config_v[1])
			{
				node.type = 0;
				vector<double> svles;
				svles.reserve(SIZE(node.temp_dids));
				REP (i, SIZE(node.temp_dids))
				{
					svles.push_back(dt_lib_ans[node.temp_dids[i]]);
				}
				node.vle = get_avg(svles);
			}
			else
			{
				vector<int> use_prp;
				int choose_size = 1;
				use_prp.reserve(choose_size);
				REP (i, choose_size)
				{
					int pos = math_rand() % SIZE(prop_pot);
					pos = prop_pot[pos];
					pick_pro_cnt[pos]++;
					use_prp.push_back(pos);
				}
				int best_cri_id = -1;
				double best_vle = -1;
				double best_eva = 1e+300;
				int a0 = 0;
				int a1 = 0;
				int b0 = 0;
				int b1 = 0;
				int best_a = 0;
				int best_b = 0;
				REP (i, SIZE(node.temp_dids))
				{
					if (dt_lib_ans[node.temp_dids[i]] > 0)
					{
						b1++;
					}
					else
					{
						b0++;
					}
				}
				int bak_a0 = a0;
				int bak_a1 = a1;
				int bak_b0 = b0;
				int bak_b1 = b1;
				if (b0 <= 1 || b1 <= 1) // may return even at bigger.
				{
					node.type = 0;
					vector<double> svles;
					svles.reserve(SIZE(node.temp_dids));
					REP (i, SIZE(node.temp_dids))
					{
						svles.push_back(dt_lib_ans[node.temp_dids[i]]);
					}
					node.vle = get_avg(svles);
				}
				else
				{
					//cout << head << " " << b1 << " " << b0 << endl;
					REP (i, SIZE(use_prp))
					{
						a0 = bak_a0;
						a1 = bak_a1;
						b0 = bak_b0;
						b1 = bak_b1;
						int cri_id = use_prp[i];
						vector<pair<double, double> > it;
						it.reserve(SIZE(node.temp_dids));
						REP (j, SIZE(node.temp_dids))
						{
							int id = node.temp_dids[j];
							it.push_back(make_pair(dt_lib_pro[id][cri_id], dt_lib_ans[id]));
						}
						sort(it.begin(), it.end());
						REP (j, SIZE(it) - 1)
						{
							if (it[j].second > 0)
							{
								a1++;
								b1--;
							}
							else
							{
								a0++;
								b0--;
							}
							if (fabs(it[j].first - it[j + 1].first) < 1e-7)
							{
								continue;
							}
							if (it[j].second == it[j + 1].second)
							{
								continue;
							}
							double eva = get_entropy(a0, a1, b0, b1);
							if (eva < best_eva)
							{
								best_eva = eva;
								best_cri_id = cri_id;
								best_vle = (it[j].first + it[j + 1].first) / 2;
								best_a = a0 + a1;
								best_b = b0 + b1;
							}
						}
					}
					if (best_eva < 1e+100)
					{
						node.type = 1;
						node.pid = best_cri_id;
						best_pro_cnt[best_cri_id]++;
						node.vle = best_vle;
						node.left_id = SIZE(nodes);
						node.right_id = SIZE(nodes) + 1;
							
						Node leftnode;
						Node rightnode;
						nodes.push_back(leftnode);
						nodes.push_back(rightnode);
							
						vector<int>& left_ids = nodes[node.left_id].temp_dids;
						vector<int>& right_ids = nodes[node.right_id].temp_dids;
						left_ids.reserve(best_a);
						right_ids.reserve(best_b);
												
						REP (i, SIZE(node.temp_dids))
						{
							if (dt_lib_pro[node.temp_dids[i]][node.pid] < node.vle)
							{
								left_ids.push_back(node.temp_dids[i]);
							}
							else
							{
								right_ids.push_back(node.temp_dids[i]);
							}
						}
					}
					else 
					{
						node.type = 0;
						vector<double> svles;
						svles.reserve(SIZE(node.temp_dids));
						REP (i, SIZE(node.temp_dids))
						{
							svles.push_back(dt_lib_ans[node.temp_dids[i]]);
						}
						node.vle = get_avg(svles);
					}
				}
			}
			vector<int> empty_vec;
			nodes[head].temp_dids.swap(empty_vec);
			head++;
		}
		#ifdef TEST
			//out_info << "Build tree with index " << day_seed << " : " << nodes.size() << endl;
		#endif
	}
};

#define START_DT_DAY 15000
#define DT_DAY_INTERVAL 5
#define MAX_DT_SIZE 1500

int now_dt_size = 0;
DT all_dt[MAX_DT_SIZE];

char tch1[100];
char tch2[100];

int read_soc_emp_int(string& str, int& pos)
{
	if (str[pos] == '_')
	{
		pos += 2;
		return -1;
	}
	int num = 0;
	while(pos < SIZE(str) && str[pos] != ' ')
	{
		num = num * 10 + (str[pos] - '0');
		++pos;
	}
	++pos;
	return num;
}

void read_soc_skip_str(string& str, int& pos)
{
	while(pos < SIZE(str) && str[pos] != ' ')
	{
		++pos;
	}
	++pos;
}

char read_soc_char(string& str, int& pos)
{
	char ch = str[pos];
	pos += 2;
	return ch;
}

int read_soc_nonemp_int(string& str, int& pos)
{
	int num = 0;
	while(pos < SIZE(str) && str[pos] != ' ')
	{
		num = num * 10 + (str[pos] - '0');
		++pos;
	}
	++pos;
	return num;
}

void setup_env()
{
	REP (i, COU_SIZE)
	{
		inner_regs[i].clear();
	}
	ZERO(in_country);
	ZERO(reg_longi);
	ZERO(reg_lati);
	
	ZERO(glb_reg_atro);
	ZERO(glb_cou_atro);
	REP (i, REG_SIZE)
	{
		glb_reg_history[i].clear();
	}
	REP (i, COU_SIZE)
	{
		glb_cou_history[i].clear();
	}
	ZERO(rct_reg_atro);
	
	ZERO(reg_past_atro);
	ZERO(cou_past_atro);
	
	ZERO(reg_past_soc);
	ZERO(cou_past_soc);
	soc_buf_ptr = 0;
	ZERO(reg_soc_buf);
	ZERO(cou_soc_buf);
	
	ZERO(soc_tot_r1);
	ZERO(soc_tot_r2);
	ZERO(soc_sta_def_times);
	ZERO(soc_sta_undef_times);
	
	ZERO(dt_lib_pro);
	ZERO(dt_prd_pro);
	ZERO(dt_lib_ans);
	dt_lib_pos = 0;
	
	ZERO(used_pos_cnt);
	ZERO(pick_pro_cnt);
	ZERO(best_pro_cnt);
	
	now_rand = 1;
		
	now_dt_size = 0;
	REP (i, MAX_DT_SIZE)
	{
		all_dt[i].construct_day = 0;
		all_dt[i].nodes.clear();
	}
}

struct MassAtrocityPredictor : public BaseMassAtrocityPredictor
{
	void generate_soc_feature(int rid, vector<double>& soc_va, vector<double>& soc_vb,
		vector<double>& soc_vc, vector<double>& soc_vd)
	{
		
		soc_va.clear();
		soc_va.reserve(SOC_CANDI_CNT);
		soc_vb.clear();
		soc_vb.reserve(SOC_CANDI_CNT);
		soc_vc.clear();
		soc_vc.reserve(SOC_CANDI_CNT);
		soc_vd.clear();
		soc_vd.reserve(SOC_CANDI_CNT);
		int cid = in_country[rid];
		
		// given type : absolute value in reg -> absolute value for world avg. 
		REP (i, EVENT_TYPE_CNT)
		{
			soc_va.push_back(reg_past_soc[0][rid][1][i]);
			soc_vb.push_back(1.0);
			soc_vc.push_back(cou_past_soc[0][COU_SIZE][1][i] * 1.0 / REG_SIZE);
			soc_vd.push_back(1.0);
		}
		
		// given type : type/all in p1 -> that in p2 
		REP (i, EVENT_TYPE_CNT)
		{
			soc_va.push_back(reg_past_soc[0][rid][0][i]);
			soc_vb.push_back(reg_past_soc[0][rid][0][0]);
			soc_vc.push_back(reg_past_soc[0][rid][1][i]);
			soc_vd.push_back(reg_past_soc[0][rid][1][0]);
		}
		
		// given type : reg/glb avg in p1 -> that in p2.
		REP (i, EVENT_TYPE_CNT)
		{
			soc_va.push_back(reg_past_soc[0][rid][0][i]);
			soc_vb.push_back(cou_past_soc[0][COU_SIZE][0][i] * 1.0 / REG_SIZE);
			soc_vc.push_back(reg_past_soc[0][rid][1][i]);
			soc_vd.push_back(cou_past_soc[0][COU_SIZE][1][i] * 1.0 / REG_SIZE);
		}
		
		// given type : reg/cou avg in p1 -> that in p2 
		REP (i, EVENT_TYPE_CNT)
		{
			soc_va.push_back(reg_past_soc[0][rid][0][i]);
			soc_vb.push_back(cou_past_soc[0][cid][0][i] * 1.0 / SIZE(inner_regs[cid]));
			soc_vc.push_back(reg_past_soc[0][rid][1][i]);
			soc_vd.push_back(cou_past_soc[0][cid][1][i] * 1.0 / SIZE(inner_regs[cid]));
		}
	}
		
	void update_library(int dayID)
	{
		vector<double> soc_va, soc_vb, soc_vc, soc_vd;
		vector<int> focus_type[2];
		REP (i, EVENT_TYPE_CNT)
		{
			focus_type[0].push_back(i);
		}
		/*
		focus_type[0].push_back(0);
		//atype begins from 11.
		focus_type[0].push_back(14);
		focus_type[0].push_back(18);
		focus_type[0].push_back(21);
		*/
		
		//focus_type[1].push_back(0);
		
		REP (i, REG_SIZE)
		{
			int cid = in_country[i];
			dt_lib_ans[dt_lib_pos] = 0;
			if (rct_reg_atro[i] > 0)
			{
				dt_lib_ans[dt_lib_pos] = 1;
			}
			
			if (dayID >= START_COLLECT_SOC_STA)
			{
				generate_soc_feature(i, soc_va, soc_vb, soc_vc, soc_vd);
				int t_fact = dt_lib_ans[dt_lib_pos];
				REP (j, SOC_CANDI_CNT)
				{
					if (soc_vb[j] == 0 || soc_vd[j] == 0)
					{
						soc_sta_undef_times[j][t_fact]++;
					}
					else
					{
						soc_tot_r1[j][t_fact] += soc_va[j] / soc_vb[j];
						soc_tot_r2[j][t_fact] += soc_vc[j] / soc_vd[j];
						soc_sta_def_times[j][t_fact]++;
					}
				}		
				/*
				//soc sta only:
				#define SOC_CANDI_CNT 3
				double soc_tot_r1[SOC_CANDI_CNT][2];
				double soc_tot_r2[SOC_CANDI_CNT][2];
				
				int soc_sta_def_times[SOC_CANDI_CNT][2];
				int soc_sta_undef_times[SOC_CANDI_CNT][2];
				*/
			}
						
			int column_pos = 0;
			REP (j, REG_CHECK_DAY_LEN)
			{
				dt_lib_pro[dt_lib_pos][column_pos] = reg_past_atro[0][i][j];
				dt_prd_pro[i][column_pos] = reg_past_atro[1][i][j];
				++column_pos;
			}
			
			double region_ratio1[2];
			double region_ratio2[2];
			double region_ratio[2];
			REP (j, 2)
			{
				region_ratio1[j] = 0;
				if (cou_past_atro[j][cid][COU_CHECK_DAY_LEN] != 0)
				{
					region_ratio1[j] = reg_past_atro[j][i][REG_CHECK_DAY_LEN] * 1.0 / cou_past_atro[j][cid][COU_CHECK_DAY_LEN];
				}
				region_ratio2[j] = 1.0 / SIZE(inner_regs[cid]);
				region_ratio[j] = region_ratio1[j] * 0.7 + region_ratio2[j] * 0.3;
			}
			
			REP (j, COU_CHECK_DAY_LEN)
			{
				dt_lib_pro[dt_lib_pos][column_pos] = cou_past_atro[0][cid][j] * 1.0 * region_ratio[0];
				dt_prd_pro[i][column_pos] = cou_past_atro[1][cid][j] * 1.0 * region_ratio[1];
				++column_pos;
			}
			{
				dt_lib_pro[dt_lib_pos][column_pos] = SIZE(inner_regs[cid]);
				dt_prd_pro[i][column_pos] = SIZE(inner_regs[cid]);
				++column_pos;
			}
			{
				dt_lib_pro[dt_lib_pos][column_pos] = region_ratio1[0];
				dt_prd_pro[i][column_pos] = region_ratio1[1];
				++column_pos;
			}
			{
				dt_lib_pro[dt_lib_pos][column_pos] = cou_past_atro[0][COU_SIZE][COU_CHECK_DAY_LEN + 1];
				dt_prd_pro[i][column_pos] = cou_past_atro[1][COU_SIZE][COU_CHECK_DAY_LEN + 1];
				++column_pos;
			}
			
			{
				int pos = SIZE(glb_reg_history[i]) - 1;
				int td = dayID - FIRST_DAY - 30;
				while (pos >= 0 && glb_reg_history[i][pos] > td)
				{
					--pos;
				}
				for (int k = 0; k < 1; ++k)
				{
					int lastt1 = 10000;
					if (pos >= k)
					{
						lastt1 = td - glb_reg_history[i][pos - k];
					}
					int lastt2 = 10000;
					if (SIZE(glb_reg_history[i]) > k)
					{
						lastt2 = dayID - FIRST_DAY - glb_reg_history[i][SIZE(glb_reg_history[i]) - 1 - k];
					}
					dt_lib_pro[dt_lib_pos][column_pos] = lastt1;
					dt_prd_pro[i][column_pos] = lastt2;
					++column_pos;
				}
			}
			{
				int pos = SIZE(glb_cou_history[cid]) - 1;
				int td = dayID - FIRST_DAY - 30;
				while (pos >= 0 && glb_cou_history[cid][pos] > td)
				{
					--pos;
				}
				for (int k = 0; k < 1; ++k)
				{
					int lastt1 = 10000;
					if (pos >= k)
					{
						lastt1 = td - glb_cou_history[cid][pos - k];
					}
					int lastt2 = 10000;
					if (SIZE(glb_cou_history[cid]) > k)
					{
						lastt2 = dayID - FIRST_DAY - glb_cou_history[cid][SIZE(glb_cou_history[cid]) - 1 - k];
					}
					dt_lib_pro[dt_lib_pos][column_pos] = lastt1;
					dt_prd_pro[i][column_pos] = lastt2;
					++column_pos;
				}
			}
			
			{
				dt_lib_pro[dt_lib_pos][column_pos] = reg_longi[i];
				dt_prd_pro[i][column_pos] = reg_longi[i];
				++column_pos;
				dt_lib_pro[dt_lib_pos][column_pos] = reg_lati[i];
				dt_prd_pro[i][column_pos] = reg_lati[i];
				++column_pos;
			}
			
			
			
			//SOC PART
			REP (k, SIZE(focus_type[0]))
			{
				double d1 = 1;
				if (cou_past_soc[0][COU_SIZE][0][focus_type[0][k]] != 0)
				{
					d1 = reg_past_soc[0][i][0][focus_type[0][k]] / (cou_past_soc[0][COU_SIZE][0][focus_type[0][k]] * 1.0 / REG_SIZE);
				}
				double d2 = 1;
				if (cou_past_soc[1][COU_SIZE][0][focus_type[0][k]] != 0)
				{
					d2 = reg_past_soc[1][i][0][focus_type[0][k]] / (cou_past_soc[1][COU_SIZE][0][focus_type[0][k]] * 1.0 / REG_SIZE);
				}
				dt_lib_pro[dt_lib_pos][column_pos] = d1;
				dt_prd_pro[i][column_pos] = d2;
				++column_pos;
			}
			REP (k, SIZE(focus_type[1]))
			{
				double d1 = 1;
				if (cou_past_soc[0][COU_SIZE][0][focus_type[1][k]] != 0)
				{
					d1 = cou_past_soc[0][cid][0][focus_type[1][k]] / (cou_past_soc[0][COU_SIZE][0][focus_type[1][k]] * 1.0 / COU_SIZE);
				}
				double d2 = 1;
				if (cou_past_soc[1][COU_SIZE][0][focus_type[1][k]] != 0)
				{
					d2 = cou_past_soc[1][cid][0][focus_type[1][k]] / (cou_past_soc[1][COU_SIZE][0][focus_type[1][k]] * 1.0 / COU_SIZE);
				}
				dt_lib_pro[dt_lib_pos][column_pos] = d1;
				dt_prd_pro[i][column_pos] = d2;
				++column_pos;
			}
			
			++dt_lib_pos;
			if (dt_lib_pos == DT_LIB_SIZE)
			{
				dt_lib_pos = 0;
			}
		}
	}
	
	void read_geo_data(vector<string>& data)
	{
		int cid;
		int rid;
		for (int i = 0; i < SIZE(data);)
		{
			sscanf(data[i].c_str(), "%d %s", &cid, tch1);
			++i;
			sscanf(data[i].c_str(), "%d %s", &rid, tch1);
			++i;
			inner_regs[cid].push_back(rid);
			in_country[rid] = cid;
			/*
			#ifdef TEST
				out_info << "Read geo: " << cid << " " << rid << endl;
			#endif
			*/
			bool first_outer = false;
			while (i < SIZE(data) && (data[i] == "outer" || data[i] == "inner"))
			{
				if (!first_outer && data[i] == "outer")
				{
					first_outer = true;
					++i;
					sscanf(data[i].c_str(), "%lf,%lf", &reg_longi[rid], &reg_lati[rid]);
					++i;
				}
				else
				{
					i += 2;
				}
			}
		}
	}
	void read_soc_data(int dayID, vector<string>& data)
	{
		/*
		int read_soc_emp_int(string& str, int& pos)
		void read_soc_skip_str(string& str, int& pos)
		char read_soc_char(string& str, int& pos)
		int read_soc_nonemp_int(string& str, int& pos)
		*/
		int day_index = dayID - FIRST_DAY;
		int cur_pos = 0;
		int p1_id, p1_prc, p1_cid, p1_rid;
		int p2_id, p2_prc, p2_cid, p2_rid;
		char ach;
		int atype;
		int a_prc, a_cid, a_rid;
		char imp_ch;
		int imp_vle;
		int med_cov, med_sen;
		vector<bool> event_flag(EVENT_TYPE_CNT, false);
		for (int i = 0; i < SIZE(data); ++i)
		{
			cur_pos = 0;
			string& str = data[i];
			p1_id = read_soc_emp_int(str, cur_pos);
			p1_prc = read_soc_emp_int(str, cur_pos);
			read_soc_skip_str(str, cur_pos);
			p1_cid = read_soc_emp_int(str, cur_pos);
			p1_rid = read_soc_emp_int(str, cur_pos);
			
			p2_id = read_soc_emp_int(str, cur_pos);
			p2_prc = read_soc_emp_int(str, cur_pos);
			read_soc_skip_str(str, cur_pos);
			p2_cid = read_soc_emp_int(str, cur_pos);
			p2_rid = read_soc_emp_int(str, cur_pos);
			
			ach = read_soc_char(str, cur_pos);
			atype = (ach - 'a');
			#ifdef AST
				assert((atype >= 0) && (atype <= 19));
			#endif
			
			a_prc = read_soc_emp_int(str, cur_pos);
			read_soc_skip_str(str, cur_pos);
			a_cid = read_soc_emp_int(str, cur_pos);
			a_rid = read_soc_emp_int(str, cur_pos);
			
			imp_ch = read_soc_char(str, cur_pos);
			#ifdef AST
				assert((imp_ch == 't') || (imp_ch == 'f'));
			#endif
			imp_vle = ((imp_ch == 't') ? 1 : 0);
			
			med_cov = read_soc_nonemp_int(str, cur_pos);
			#ifdef AST
				assert(med_cov >= 1 && med_cov <= 100);
			#endif
			med_sen = read_soc_nonemp_int(str, cur_pos);
			#ifdef AST
				assert(med_sen >= 0 && med_sen <= 50);
			#endif
			
			int peo_cnt = 0;
			if (p1_id >= 0)
			{
				++peo_cnt;
			}
			if (p2_id >= 0)
			{
				++peo_cnt;
			}
			event_flag[0] = true;
			event_flag[1] = (imp_vle != 0);
			event_flag[2] = (!event_flag[1]);
			event_flag[3] = (peo_cnt == 1);
			event_flag[4] = (peo_cnt == 2);
			event_flag[5] = (a_cid >= 0 && ((p1_cid >= 0 && p1_cid != a_cid) || (p2_cid >= 0 && p2_cid != a_cid)));
			event_flag[6] = (!event_flag[5]);
			event_flag[7] = (a_rid >= 0 && ((p1_rid >= 0 && p1_rid != a_rid) || (p2_rid >= 0 && p2_rid != a_rid)));
			event_flag[8] = (!event_flag[7]);
			event_flag[9] = (med_cov > 15);
			event_flag[10] = (!event_flag[9]);
			REP (j, 20)
			{
				event_flag[11 + j] = (atype == j);
			}
			REP (j, EVENT_TYPE_CNT)
			{
				if (event_flag[j])
				{
					if (a_rid >= 0)
					{
						reg_soc_buf[a_rid][soc_buf_ptr][j]++;
					}
					if (a_cid >= 0)
					{
						cou_soc_buf[a_cid][soc_buf_ptr][j]++;
					}
					cou_soc_buf[COU_SIZE][soc_buf_ptr][j]++;
				}
			}
					
			/*
			#ifdef AST
			out_info << p1_id << " " << p1_prc << " " << p1_cid << " " << p1_rid
			  << " " << p2_id << " " << p2_prc << " " << p2_cid << " " << p2_rid
			  << " " << atype << " " << a_prc << " " << a_cid << " " << a_rid
			  << " " << imp_vle << " " << med_cov << " " << med_sen << endl;
			#endif
			*/
		}
		REP (ii, 2)
		{
			REP (i, REG_EVENT_CHECK_LEN)
			{
				int period = reg_eve_check_day[i];
				if (day_index >= fix_day_len[ii])
				{
					int target_ptr = soc_buf_ptr - fix_day_len[ii];
					if (target_ptr < 0)
					{
						target_ptr += SOC_BUF_SIZE;
					}
					REP (j, REG_SIZE)
					{
						REP (k, EVENT_TYPE_CNT)
						{
							reg_past_soc[ii][j][i][k] += reg_soc_buf[j][target_ptr][k];
						}
					}
				}
				if (day_index >= fix_day_len[ii] + period)
				{
					int target_ptr = soc_buf_ptr - fix_day_len[ii] - period;
					if (target_ptr < 0)
					{
						target_ptr += SOC_BUF_SIZE;
					}
					REP (j, REG_SIZE)
					{
						REP (k, EVENT_TYPE_CNT)
						{
							reg_past_soc[ii][j][i][k] -= reg_soc_buf[j][target_ptr][k];
						}
					}
				}
			}
			REP (i, COU_EVENT_CHECK_LEN)
			{
				int period = cou_eve_check_day[i];
				if (day_index >= fix_day_len[ii])
				{
					int target_ptr = soc_buf_ptr - fix_day_len[ii];
					if (target_ptr < 0)
					{
						target_ptr += SOC_BUF_SIZE;
					}
					REP (j, COU_SIZE + 1)
					{
						REP (k, EVENT_TYPE_CNT)
						{
							cou_past_soc[ii][j][i][k] += cou_soc_buf[j][target_ptr][k];
						}
					}
				}
				if (day_index >= fix_day_len[ii] + period)
				{
					int target_ptr = soc_buf_ptr - fix_day_len[ii] - period;
					if (target_ptr < 0)
					{
						target_ptr += SOC_BUF_SIZE;
					}
					REP (j, COU_SIZE + 1)
					{
						REP (k, EVENT_TYPE_CNT)
						{
							cou_past_soc[ii][j][i][k] -= cou_soc_buf[j][target_ptr][k];
						}
					}
				}
			}
		}
		soc_buf_ptr++;
		if (soc_buf_ptr == SOC_BUF_SIZE)
		{
			soc_buf_ptr = 0;
		}
		//clear old data.
		REP (i, REG_SIZE)
		{
			REP (j, EVENT_TYPE_CNT)
			{
				reg_soc_buf[i][soc_buf_ptr][j] = 0;
			}
		}
		REP (i, COU_SIZE + 1)
		{
			REP (j, EVENT_TYPE_CNT)
			{
				cou_soc_buf[i][soc_buf_ptr][j] = 0;
			}
		}
		/*
		cout << endl;
		cout << cou_past_soc[1][COU_SIZE][0][0] << " " << cou_past_soc[1][COU_SIZE][0][1] << " " << cou_past_soc[1][COU_SIZE][0][2] << endl;
		cout << cou_past_soc[0][COU_SIZE][0][0] << " " << cou_past_soc[0][COU_SIZE][0][1] << " " << cou_past_soc[0][COU_SIZE][0][2] << endl;
		*/
	}
	
	void read_atro_data(int dayID, vector<string>& data)
	{
		int day_index = dayID - FIRST_DAY;
		int cid;
		int rid;
		for (int i = 0; i < SIZE(data); ++i)
		{
			sscanf(data[i].c_str(), "%s %s %d %d", tch1, tch2, &cid, &rid);
			glb_reg_atro[rid][day_index]++;
			glb_cou_atro[cid][day_index]++;
			glb_reg_history[rid].push_back(day_index);
			glb_cou_history[cid].push_back(day_index);
		}
		
		//dy update
		REP (i, REG_SIZE)
		{
			rct_reg_atro[i] += glb_reg_atro[i][day_index];
			if (day_index >= 30)
			{
				rct_reg_atro[i] -= glb_reg_atro[i][day_index - 30];
			}
		}
		
		REP (ii, 2)
		{
			REP (i, REG_CHECK_DAY_LEN + 1)
			{
				int period = reg_check_day[i];
				if (day_index >= fix_day_len[ii])
				{
					REP (j, REG_SIZE)
					{
						reg_past_atro[ii][j][i] += glb_reg_atro[j][day_index - fix_day_len[ii]];
					}
				}
				if (day_index >= fix_day_len[ii] + period)
				{
					REP (j, REG_SIZE)
					{
						reg_past_atro[ii][j][i] -= glb_reg_atro[j][day_index - fix_day_len[ii] - period];
					}
				}
			}
			REP (i, COU_CHECK_DAY_LEN + 2)
			{
				int period = cou_check_day[i];
				if (day_index >= fix_day_len[ii])
				{
					REP (j, COU_SIZE)
					{
						cou_past_atro[ii][j][i] += glb_cou_atro[j][day_index - fix_day_len[ii]];
						cou_past_atro[ii][COU_SIZE][i] += glb_cou_atro[j][day_index - fix_day_len[ii]];
					}
				}
				if (day_index >= fix_day_len[ii] + period)
				{
					REP (j, COU_SIZE)
					{
						cou_past_atro[ii][j][i] -= glb_cou_atro[j][day_index - fix_day_len[ii] - period];
						cou_past_atro[ii][COU_SIZE][i] -= glb_cou_atro[j][day_index - fix_day_len[ii] - period];
					}
				}
			}
		}
		update_library(dayID);
		
		#ifndef SIMPLE_VER
		{
			// build dt
			if (dayID >= START_DT_DAY && dayID % DT_DAY_INTERVAL == 0)
			{
				all_dt[now_dt_size].build_dt(dayID);
				all_dt[now_dt_size].construct_day = dayID;
				now_dt_size++;
			}
		}
		#endif
	}
	int receiveData(int dataSourceId, int dayID, vector<string> data)
	{
		if (dataSourceId == 0)
		{
			read_atro_data(dayID, data);
		}
		else if (dataSourceId == 1)
		{
			read_soc_data(dayID, data);
		}
		else if (dataSourceId == 2)
		{
			setup_env();
			read_geo_data(data);
		}
		return 0;
	}
	vector<double> predictAtrocities(int dayID)
	{
		int last_dt_idx = now_dt_size - 1;
		int first_dt_idx = max(last_dt_idx - 9999, 0);
		vector<double> res(REG_SIZE, 0);
		vector<double> props(DT_PRO_SIZE, 0);
		REP (i, REG_SIZE)
		{
			REP (j, DT_PRO_SIZE)
			{
				props[j] = dt_prd_pro[i][j];
			}
			double sum_wt = 0;
			for (int j = first_dt_idx; j <= last_dt_idx; ++j)
			{
				double tv = all_dt[j].get_value(props);
				int from_now = dayID - all_dt[j].construct_day;
				double wt = max(1.0 - from_now * config_v[2], config_v[3]);
				sum_wt += wt;
				res[i] += tv * wt;
			}
			res[i] /= sum_wt;
			res[i] *= thre_ratio[get_conf_id(res[i])];
			res[i] = max(res[i], 0.0);
			res[i] = min(res[i], 1.0);
		}
		return res;
	}
};

}