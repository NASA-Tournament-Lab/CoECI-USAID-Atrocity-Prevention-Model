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

#include <vector>
#include <valarray>
#include <list>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <stack>
#include <algorithm>
#include <sstream>
#include <bitset>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <complex>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef WIN32
	#include <windows.h>
	#include <time.h>
	#include <conio.h>
	#pragma warning (disable:4996)
	#define Sleep(x) 	Sleep(x)
	#define popen	_popen
	#define pclose	_pclose
	#define getpid() 0
#else
	#define getch()
	#include <unistd.h>
#endif


using namespace std;

namespace place5 {

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define EPSILON 0.00000000000001



#define COEF_NX				0
#define ADD_NX				1
#define AVG_REGION			2
#define DAY_AVG_REGION		3
#define MIN_CONFIANCE		4

#define COEF_NY				5
#define ADD_NY				6
#define AVG_COUNTRY			7
#define DAY_AVG_COUNTRY		8

#define COEF_PROBA_1		9

#define COEF_RET			10

#define COEF_NZ				11
#define ADD_NZ				12
#define COEF_PROBA_2		13
#define DAY_NZ				14

#define AVG_REGION_2		15
#define DAY_AVG_REGION_2	16

#define NB_COEF				17

#define MAX_NB_PROBA_TYPE	3
#define MAX_NB_PROBA		100

double proba_tab[MAX_NB_PROBA_TYPE][MAX_NB_PROBA];

// Obtained with training data :
double COEF[NB_COEF] = {
5.95760655251397,3.9630680431128,0.9962580703125,9278.41321108638,9.21988499749627,2.05419334625942,21.826013621842,0.998836458864393,64.1139416992664,1,0.87,60.878849821296,100,0.344239511426576,214.568233390546,1,213.766133761389};
int nb_next[MAX_NB_PROBA_TYPE][MAX_NB_PROBA] = {
{12549,0,0,1864,3529,2396,1571,1558,1728,2318,1104,903,1078,738,772,795,548,555,485,349,511,423,309,213,184,211,230,280,233,222,184,163,146,161,139,149,102,83,85,65,59,44,36,46,50,43,57,62,58,53,36,39,38,29,21,16,27,26,20,25,17,19,23,20,14,19,22,23,23,23,23,16,11,18,16,9,7,10,13,4,2,3,5,14,12,11,12,11,13,9,6,6,9,12,13,8,8,10,13,628},
{15765,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7876,0,4962,0,2713,0,1900,6,1023,12,827,15,624,7,440,29,396,16,248,6,350,11,210,64,275,0,298,0,359,14,234,47,186,0,48,0,32,0,82,0,35,0,131,0,148,3,90,15,150,10,109,22,112,10,96,5,73,12,55,36,76,43,100,42,111,12,121,60,63,6,12,0,20,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,43,467,40343}};

int nb_trig[MAX_NB_PROBA_TYPE][MAX_NB_PROBA] = {
{14888596,0,0,325919,350570,103581,70336,55001,46627,40877,16047,12254,10281,8709,7348,5810,4359,3882,3202,2725,2447,1956,1655,1342,1151,1059,873,815,738,680,645,591,510,456,389,366,322,285,269,212,186,163,146,138,121,112,106,106,99,98,76,72,66,58,54,42,51,47,40,45,36,35,42,37,37,38,39,36,32,36,35,28,21,25,24,17,12,13,17,7,6,6,9,17,15,15,17,17,19,13,12,12,15,18,19,14,14,15,16,645},
{14505202,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,880356,0,291771,0,121483,0,71678,147,33827,24,19373,30,10800,53,7713,90,5656,76,2379,15,3252,28,2074,189,2316,6,1852,0,2122,101,1729,168,1752,0,198,0,264,0,438,0,354,0,960,0,828,18,546,54,1050,114,648,228,846,36,378,18,252,36,270,108,306,126,324,144,378,36,450,180,234,18,36,0,54,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,125,215,328,4103,15971421}};


void split(const std::string &s, char delim, std::vector<std::string> &output)
{
	typedef string::const_iterator iter;
	iter beg;
	bool in_token = false;
	for( string::const_iterator it = s.begin(), end = s.end();it != end; ++it )
	{
		if( *it==delim )
		{
			if( in_token )
			{
				output.push_back(std::vector<std::string>::value_type(beg, it));
				in_token = false;
			}
		}
		else if( !in_token )
		{
			beg = it;
			in_token = true;
		}
	}
	if( in_token )
		output.push_back(std::vector<std::string>::value_type(beg, s.end()));
}

#define MAX_COUNTRY	254
#define MAX_REGION	3671

map<int,int> country_att[MAX_COUNTRY]; // array of map nb atoricities of one day by country
map<int,int> region_att[MAX_REGION];   // array of map nb atoricities of one day by region
map<int,int> region_lastatt[MAX_REGION];

map<int,double> region_avgatt[MAX_REGION];
map<int,double> region_avgatt2[MAX_REGION];
map<int,double> country_avgatt[MAX_COUNTRY];

int country_of_region[MAX_REGION];


int firstday=-1;


void reset_probas()
{
	for(int i=0;i<MAX_NB_PROBA_TYPE;i++)
	{
		for(int j=0;j<MAX_NB_PROBA;j++)
		{
			nb_trig[i][j] = 0;
			nb_next[i][j] = 0;
		}
	}
}

inline int getNTabProba(int type, int day, int region, int country)
{
	if(type==0)
	{
		if(region_avgatt[region][day]>0.0)
			return max(0,min((int)(region_avgatt[region][day]*COEF[COEF_NX] + COEF[ADD_NX]), MAX_NB_PROBA-1));
		return 0;
	}
	if(type==1)
	{
		if(country_avgatt[country][day]>0.0)
			return max(0,min((int)(country_avgatt[country][day]*COEF[COEF_NY] + COEF[ADD_NY]), MAX_NB_PROBA-1));
		return 0;
	}
	if(type==2)
	{
		int d2 = day-(int)COEF[DAY_NZ];
		if(d2<0)
			d2=0;
		double p = (double)(region_avgatt2[region][day]-region_avgatt2[region][d2])/(double)(day-d2);
		return max(0,min((int)(p*COEF[COEF_NZ] + COEF[ADD_NZ]), MAX_NB_PROBA-1));
	}

	return 0;
}

bool extractProbas = false;

class MassAtrocityPredictor : public BaseMassAtrocityPredictor
{
public:
	int receiveData(int dataSourceId, int dayID, vector<string> data)
	{
		if(firstday<0)
		{
			firstday = dayID;
			//reset_probas();
		}

		if(dataSourceId==0)
		{	
			for(vector<string>::iterator it = data.begin(); it != data.end(); ++it)
			{
				vector<string> v;
				v.reserve(4);
				split(*it,' ',v);
				int country = atoi(v[2].c_str());
				int region  = atoi(v[3].c_str());
				country_att[country][dayID]++;
				region_att[region][dayID]++;

				region_lastatt[region][dayID]=dayID;

				double avg = 1;
				for(int d=dayID;d<=dayID+COEF[DAY_AVG_REGION];d++)
				{
					region_avgatt[region][d] += avg;
					avg*=COEF[AVG_REGION];
				}

				avg = 1;
				for(int d=dayID;d<=dayID+COEF[DAY_AVG_COUNTRY];d++)
				{
					country_avgatt[country][d] += avg;
					avg*=COEF[AVG_COUNTRY];
				}

				avg = 1;
				for(int d=dayID;d<=dayID+COEF[DAY_AVG_REGION_2];d++)
				{
					region_avgatt2[region][d] += avg;
					avg*=COEF[AVG_REGION_2];
				}
			}
			
			for(int region=0;region<MAX_REGION;region++)
			{
				if(region_lastatt[region][dayID]==0)
					region_lastatt[region][dayID]=region_lastatt[region][dayID-1];
			}
			
			if(dayID>firstday)
			{
				// save samples for probas
				if(dayID>firstday+30)
				{
					for(int region=0;region<MAX_REGION;region++)
					{
						int country = country_of_region[region];
						int day = dayID-30;
						
						int	nbx = getNTabProba(0, day, region, country);
						int	nby = getNTabProba(1, day, region, country);
						int nbz = getNTabProba(2, day, region, country);
						
						if(region_lastatt[region].find(day+30)!=region_lastatt[region].end() && region_lastatt[region][day+30]>=day+1)
						{
							nb_next[0][nbx]++;
							nb_next[1][nby]++;
							nb_next[2][nbz]++;
						}
						
						nb_trig[0][nbx]++;
						nb_trig[1][nby]++;
						nb_trig[2][nbz]++;
					}
				}

				// clears map to save memory
				if(dayID>firstday+120)
				{
					for(int region=0;region<MAX_REGION;region++)
					{
						region_avgatt[region].erase(dayID-120);
						region_avgatt2[region].erase(dayID-((int)COEF[DAY_NZ]+30));
						region_att[region].erase(dayID-120);
						region_lastatt[region].erase(dayID-120);
					}
					for(int country=0;country<MAX_COUNTRY;country++)
					{
						country_avgatt[country].erase(dayID-120);
						country_att[country].erase(dayID-120);
					}
				}
			}
		}
		else if(dataSourceId==1)
		{
			
		}
		else if(dataSourceId==2)
		{
			for(vector<string>::iterator it = data.begin(); it != data.end(); ++it)
			{
				if( strstr(it->c_str(),"inner") ||
					strstr(it->c_str(),"outer") ||
					strchr(it->c_str(),',') )
					continue;

				int country = atoi(it->c_str());
				++it;
				int region = atoi(it->c_str());
				country_of_region[region] = country;
			}
		}

		return 0;
	}

	vector<double> predictAtrocities(int dayID)
	{
		vector<double> ret(MAX_REGION);

		for(int region=0;region<MAX_REGION;region++)
		{
			int country = country_of_region[region];

			// p1
			double p1;

			int	nbx = getNTabProba(0, dayID, region, country);

			while((nb_trig[0][nbx]<=0.0 || nb_trig[0][nbx]<COEF[MIN_CONFIANCE]) && nbx>=0)
				nbx--;

			if(nbx<0)
			{
				p1 = 0.0;
			}
			else
			{
				p1 = (double)nb_next[0][nbx]/(double)nb_trig[0][nbx];
			}

			// p2
			double p2;

			int	nby = getNTabProba(1, dayID, region, country);

			while((nb_trig[1][nby]<=0.0 || nb_trig[1][nby]<COEF[MIN_CONFIANCE]) && nby>=0)
				nby--;

			if(nby<0)
			{
				p2 = 0.0;
			}
			else
			{
				p2 = (double)nb_next[1][nby]/(double)nb_trig[1][nby];
			}

			// p3
			double p3;

			int	nbz = getNTabProba(2, dayID, region, country);

			while((nb_trig[2][nbz]<=0.0 || nb_trig[2][nbz]<COEF[MIN_CONFIANCE]) && nbz>=0)
				nbz--;

			if(nbz<0)
			{
				p3 = 0.0;
			}
			else
			{
				p3 = (double)nb_next[2][nbz]/(double)nb_trig[2][nbz];
			}

			double p = p1*COEF[COEF_PROBA_1] + p2*(1.0*COEF[COEF_PROBA_2]) + p3*(1.0-(COEF[COEF_PROBA_1]+COEF[COEF_PROBA_2]));
			ret[region] = max(0.0,min(p*COEF[COEF_RET],1.0));
		}

		return ret;
	}
};


}