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

// File Name: code_v0_0.cpp
// Author: ***
// Created Time: Fri Aug 23 21:53:21 2013

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <bitset>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <cstring>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cmath>
//#define LOCAL 
using namespace std;
namespace place3 {
template<class T> string toString(T n) { ostringstream ost; ost<<n; ost.flush(); return ost.str(); }
struct SRoll
{//a struct to preserve a queue of qsize, then calculate the sum this queue in a rolling way
  int *q;
  int num;
  int count;//how many 1
  int top;
  int qsize;
  SRoll()
  {
    qsize = 0;
    q = NULL;
    num = count = top = 0;
  }
  SRoll(int _qsize)
  {
    qsize = _qsize;
    q = new int [qsize];
    num = count = top = 0;
  }
  SRoll(const struct SRoll &that)
  {
    num = that.num;
    count = that.num;
    top = that.top;
    qsize = that.qsize;
    if (qsize == 0) {
      q = NULL;
    } else {
      q = new int[qsize];
      memcpy(q, that.q, qsize * sizeof(int));
    }
  }
  ~SRoll() { if (q != NULL) delete[]q; }
  int getSize() { return num; }
  bool isFull() { return (num == qsize); }
  int getCount() { return count; }
  void print() 
  {
    for (int i = (top - num + qsize) % qsize, c = 0; c < num;
        i = (i + 1 + qsize) % qsize, c++) {
      cout<<q[i]<<' ';
    }
    cout<<endl;
    //cin.get();
  }
  void push(int mark)
  {
    if (num == qsize) {
      count -= q[top];
      q[top] = mark;
      count += mark;
      top = (top + 1) % qsize;
    } else {
      count += mark;
      q[top] = mark;
      top = (top + 1) % qsize;
      num++;
    }
  }
};
class CPredictor2
{
 //just look back for kWindowSize day, and calculate the probability, because, this will overestimate the probability
 //use a factor to adjust it
 //it use the same algo as CPredictor3, difference is that using a decay factor to calculate the probability, weight more for recent days
  static const int kRegionNum = 3671;
  static const int kCountryNum = 254;
  static const int kReserveDay = 30;
  static const int kFeatureNum = 9;
  //static const int kWindowSize0 = 100;
  //static const int kWindowSize1 = 500;
  double factor;
  int kWindowSize1;
  double ow;
  double w;
  int region_num;
  deque<vector<int > > hist_sociopolitical_data;
  vector<int> last30day; 
  vector<int> last_di;
  double prob1[kRegionNum];
  deque<vector<int> > hist_atrocities;
  deque<vector<int> > hist_windows0;//100
  deque<vector<int> > hist_windows1;//500
  vector<SRoll> queue_atrocity;
public:
  CPredictor2(double f, int win_size):factor(f), kWindowSize1(win_size)
  {
    memset(prob1, 0, sizeof(double) * kRegionNum);
    w = 0.998;
    ow = 1;
    for (int i = 0;i < kWindowSize1;++i) {
      ow *= w;
    }
    queue_atrocity = vector<SRoll>(kRegionNum, SRoll(30));
  }
  ~CPredictor2()
  {
  }
  vector<int> readAtrocities_data(vector <string> &data)
  {
    vector<int> tmp_region;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      tmp_region.push_back(rid);
    }
    return tmp_region;
  }
  int ID(int a, int b)
  {
    return (a + 1) * (1 + kCountryNum) + (b + 1);
  }
  void getLast30Day(vector<int> &pid)
  {
    deque<vector<int > >::iterator it;
    for (it = hist_sociopolitical_data.begin();
        it != hist_sociopolitical_data.end();++it) {
      pid.insert(pid.end(), (*it).begin(), (*it).end() );
    }
  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      vector<int> rid = readAtrocities_data(data);
      vector<int> tmp_aid(kRegionNum, 0);
      for (int i = 0;i < rid.size();++i) {
        tmp_aid[rid[i]] += 1;
      }
      for (int i = 0;i < kRegionNum;++i) {
        queue_atrocity[i].push(tmp_aid[i]);
      }
      hist_atrocities.push_back(tmp_aid);
      if (hist_atrocities.size() > kReserveDay) {
        hist_atrocities.pop_front();
      }
      vector<int> count(kRegionNum, 0);
      deque<vector<int> >::iterator it;
      for (it = hist_atrocities.begin();
          it != hist_atrocities.end();++it) {
        for (int j = 0;j < (*it).size();++j) {
          if ((*it)[j]) {
            count[j] = 1;
          }
        }
      }
      hist_windows1.push_back(count);
      for (int j = 0;j < (kRegionNum);++j) {
        prob1[j] = w * prob1[j] +  (double)count[j];
      }
      if (hist_windows1.size() > kWindowSize1) {
        vector<int> &old = hist_windows1.front();
        for (int j = 0;j < kRegionNum;++j) {
          prob1[j] -= ow * (double)old[j];
        }
        hist_windows1.pop_front();
      }
    } else if(dataSourceId == 1) {//1
    } else if(dataSourceId == 2) {//0
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    double wn = (1 - ow) / (1 - w);
    vector<double> ans(region_num, 0.0);
    vector<pair<int, int > > prv;
    for (int i = 0;i < kRegionNum;++i) {
      int num = queue_atrocity[i].getCount();
      prv.push_back(make_pair(num, i));
    }
    sort(prv.begin(), prv.end(), greater<pair<int, int> > ());
    vector<double> weight(region_num, 0.0);
    for (int i = 0;i < kRegionNum;++i) {
      int num = queue_atrocity[i].getCount();
      //if (num > 10) num = 10;
      weight[i] = tanh((double)(num) / 3.0) + 1.0;
    }

    for (int i = 0;i < kRegionNum;++i) {
      double ans1 =  (double) prob1[i] / (double) wn;
      if (ans1 <= 0) ans1 = 0.0;
      if (ans1 >= 1) {
        ans1 = 1;
      }
      //assert(ans1 <= 1 && ans1 >= 0);
      ans[i] = ans1 * factor;//multiple by a factor
      ans[i] = ans[i] * weight[i];
      if (ans[i] >= 1) ans[i] = 1;
      if (ans[i] <= 0) ans[i] = 0;
    }
    return ans;
  }
};
class CPredictor3
{//suppose last kWindowSize days probability is the same, same algo as CPredictor2, without decay.
  static const int kRegionNum = 3671;
  static const int kCounryNum = 254;
  static const int kReserveDay = 30;
  static const int kFeatureNum = 9;
  static const int kWindowSize0 = 100;
  double factor;
  int kWindowSize1;
  //static const int kWindowSize1 = 500;
  int region_num;
  deque<vector<int > > hist_sociopolitical_data;
  vector<int> last30day; 
  vector<int> last_di;
  int prob0[kRegionNum];
  int prob1[kRegionNum];
  deque<vector<int> > hist_atrocities;
  deque<vector<int> > hist_windows0;//100
  deque<vector<int> > hist_windows1;//500
public:
  CPredictor3(double f, int win_size1):
    factor(f), kWindowSize1(win_size1)
  {
    memset(prob0, 0, sizeof(int) * kRegionNum);
    memset(prob1, 0, sizeof(int) * kRegionNum);
  }
  ~CPredictor3()
  {
  }
  vector<int> readAtrocities_data(vector <string> &data)
  {
    vector<int> tmp_region;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      tmp_region.push_back(rid);
    }
    return tmp_region;
  }
  int ID(int a, int b)
  {
    return (a + 1) * (1 + kCounryNum) + (b + 1);
  }
  void getLast30Day(vector<int> &pid)
  {
    deque<vector<int > >::iterator it;
    for (it = hist_sociopolitical_data.begin();
        it != hist_sociopolitical_data.end();++it) {
      pid.insert(pid.end(), (*it).begin(), (*it).end() );
    }
  }

  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      vector<int> rid = readAtrocities_data(data);
      vector<int> tmp_aid(kRegionNum, 0);
      for (int i = 0;i < rid.size();++i) {
        tmp_aid[rid[i]] = 1;
      }
      hist_atrocities.push_back(tmp_aid);
      if (hist_atrocities.size() > kReserveDay) {
        hist_atrocities.pop_front();
      }
      vector<int> count(kRegionNum, 0);
      deque<vector<int> >::iterator it;
      for (it = hist_atrocities.begin();
          it != hist_atrocities.end();++it) {
        for (int j = 0;j < (*it).size();++j) {
          if ((*it)[j]) {
            count[j] = 1;
          }
        }
      }
      hist_windows0.push_back(count);
      for (int j = 0;j < (kRegionNum);++j) {
        prob0[j] += count[j];
      }
      hist_windows1.push_back(count);
      for (int j = 0;j < (kRegionNum);++j) {
        prob1[j] += count[j];
      }
      if (hist_windows0.size() > kWindowSize0) {
        vector<int> &old = hist_windows0.front();
        for (int j = 0;j < kRegionNum;++j) {
          prob0[j] -= old[j];
        }
        hist_windows0.pop_front();
      }
      if (hist_windows1.size() > kWindowSize1) {
        vector<int> &old = hist_windows1.front();
        for (int j = 0;j < kRegionNum;++j) {
          prob1[j] -= old[j];
        }
        hist_windows1.pop_front();
      }
    } else if(dataSourceId == 1) {//1
    } else if(dataSourceId == 2) {//0
      //readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    for (int i = 0;i < kRegionNum;++i) {
      double ans0 =  (double) prob0[i] / (double) hist_windows0.size();
      double ans1 =  (double) prob1[i] / (double) hist_windows1.size();
      ans[i] = (0.1 * ans0 + 0.9 * ans1) * factor;
    }
    return ans;
  }
};

class CPredictor4
{//calculate the probability in a country for last kWindowSize days, and use it to estimate probability in each region
 //it is useful when the country is big.
  static const int kRegionNum = 3671;
  static const int kCountryNum = 254;
  static const int kReserveDay = 30;
  static const int kFeatureNum = 9;
  //static const int kWindowSize1 = 500;
  double factor;
  int kWindowSize1;
  double ow;
  double w;
  int region_num;
  int country_num;
  double prob1[kCountryNum];
  deque<vector<int> > hist_atrocities;
  deque<vector<int> > hist_windows1;//500
  vector<SRoll> queue_atrocity;
  vector<int> rig_count;
public:
  CPredictor4(double f, int win_size):factor(f), kWindowSize1(win_size)
  {
    memset(prob1, 0, sizeof(double) * kCountryNum);
    w = 0.998;
    ow = 1;
    for (int i = 0;i < kWindowSize1;++i) {
      ow *= w;
    }
    queue_atrocity = vector<SRoll>(kCountryNum, SRoll(30));
    rig_count = vector<int>(kCountryNum, 0);
  }
  ~CPredictor4()
  {
  }
  vector<int> readAtrocities_data(vector <string> &data)
  {
    vector<int> tmp_region;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      tmp_region.push_back(rid);
    }
    return tmp_region;
  }
  map<int, int> cr_map;//country region map
  void readRegion(vector<string> &data)
  {
    vector<int> reg_start_point;
    for (int i = 0;i < data.size();) {
      reg_start_point.push_back(i);
      istringstream sin1(data[i]);
      istringstream sin2(data[i + 1]);
      int cid, rid;
      sin1>>cid;//assert(cid >= 0 && cid < kCountryNum);
      sin2>>rid;//assert(rid >= 0 && rid < kRegionNum);
      cr_map[rid] = cid;
      rig_count[cid] += 1;
      i += 2;
      while (i < data.size() && (data[i][0] == 'o' || data[i][0] == 'i') ) {
        i += 2;
      }
    }
  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      vector<int> rid = readAtrocities_data(data);
      vector<int> tmp_aid(kCountryNum, 0);
      for (int i = 0;i < rid.size();++i) {
        tmp_aid[cr_map[rid[i]]] += 1;
      }
      for (int i = 0;i < kCountryNum;++i) {
        queue_atrocity[i].push(tmp_aid[i]);
      }
      hist_atrocities.push_back(tmp_aid);
      if (hist_atrocities.size() > kReserveDay) {
        hist_atrocities.pop_front();
      }
      vector<int> count(kCountryNum, 0);
      deque<vector<int> >::iterator it;
      for (it = hist_atrocities.begin();
          it != hist_atrocities.end();++it) {
        for (int j = 0;j < (*it).size();++j) {
          if ((*it)[j]) {
            count[j] = 1;
          }
        }
      }
      hist_windows1.push_back(count);
      for (int j = 0;j < (kCountryNum);++j) {
        prob1[j] = w * prob1[j] +  (double)count[j];
      }
      if (hist_windows1.size() > kWindowSize1) {
        vector<int> &old = hist_windows1.front();
        for (int j = 0;j < kCountryNum;++j) {
          prob1[j] -= ow * (double)old[j];
        }
        hist_windows1.pop_front();
      }
    } else if(dataSourceId == 1) {//1
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    double wn = (1 - ow) / (1 - w);
    vector<double> ans(region_num, 0.0);

    for (int i = 0;i < kRegionNum;++i) {
      int cid = cr_map[i];
      double div = (rig_count[cid] > 0)?(double)(rig_count[cid]):1.0;
      double ans1 =  (double) prob1[cid] / (double) wn / div;
      if (ans1 <= 0) ans1 = 0.0;
      if (ans1 >= 1) {
        ans1 = 1;
      }
      ans[i] = ans1 * factor;//multiple by a factor
      if (ans[i] >= 1) ans[i] = 1;
      if (ans[i] <= 0) ans[i] = 0;
    }
    return ans;
  }
};
class CPredictor5
{//conditioned on whether atrocity happened for the last kBackDay days, in a region,
 //adjust kBackDay to tune the model
  static const int kRegionNum = 3671;
  static const int kWindowSize = 700;
  static const int kEventType = 20;
  static const int kReserveDay = 30;//di + 1,di + 30
  int kBackDay;
  vector<SRoll> queue_atrocity_big;
  vector<SRoll> queue_atrocity_small;
  vector<SRoll> queue_atrocity;//this is for prediction
  deque<vector<bool> > qab[2];
  deque<vector<bool> > qa[2];
  int pab[kRegionNum][2];
  int pa[kRegionNum][2];
  vector<vector<bool> > cur_sociopolitical;
  vector<bool> cur_atrocities;
  double factor;
public:
  CPredictor5(int backday, double f)
  {
    kBackDay = backday;
    factor = f;
    for (int i = 0;i < kRegionNum;++i) {
      pab[i][0] = pa[i][0] = 0;
      pab[i][1] = pa[i][1] = 0;
    }
    queue_atrocity_big = vector<SRoll>(kRegionNum, SRoll(kReserveDay + kBackDay));
    queue_atrocity_small = vector<SRoll>(kRegionNum, SRoll(kReserveDay));
    queue_atrocity = vector<SRoll>(kRegionNum, SRoll(kBackDay));

  }
  ~CPredictor5()
  {
  }
  void readRegion(vector<string> &data)
  {
  }
  vector<bool> readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
    }

    vector<bool> tmp_aid(kRegionNum, false);
    for (int i = 0;i < region_id.size();++i) {
      tmp_aid[region_id[i]] = true;
    }
    return tmp_aid;
  }
  void updateRollingArray()
  {
    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity[i].push(cur_atrocities[i]);
      queue_atrocity_big[i].push(cur_atrocities[i]);
      queue_atrocity_small[i].push(cur_atrocities[i]);
    }
    if (queue_atrocity_big[0].isFull()) {//if is full
      vector<bool> tmp_atrocity(kRegionNum, false);
      for (int i = 0;i < kRegionNum;++i) {
        if (queue_atrocity[i].getCount() > 0) {
          tmp_atrocity[i] = true;
        }
      }
      for (int j = 0;j <= 1;++j) {
        vector<bool> tmp_qab(kRegionNum, false);
        vector<bool> tmp_qa(kRegionNum, false);
        for (int i = 0;i < kRegionNum;++i) {
          int event_num = queue_atrocity_big[i].getCount() - 
            queue_atrocity_small[i].getCount();
          int type = (event_num > 0)?1:0;
          if (type == j) {
            tmp_qa[i] = true;
            if (tmp_atrocity[i]) {
              tmp_qab[i] = true;
            }
          }
        }
        for (int i = 0;i < kRegionNum;++i) {
          if (tmp_qab[i]) {
            pab[i][j] += 1;
          }
          if (tmp_qa[i]) {
            pa[i][j] += 1;
          }
        }
        qab[j].push_back(tmp_qab);
        qa[j].push_back(tmp_qa);
        if (qab[j].size() > kWindowSize) {
          vector<bool> &tv = qab[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pab[i][j] -= 1;
            }
          }
          qab[j].pop_front();
        }
        if (qa[j].size() > kWindowSize) {
          vector<bool> &tv = qa[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pa[i][j] -= 1;
            }
          }
          qa[j].pop_front();
        }
      }
    }
  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      cur_atrocities.clear();
      cur_atrocities = readAtrocities_data(dayID, data);
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      //cur_sociopolitical.clear();
      //cur_sociopolitical = readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }

  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);

    for (int i = 0;i < region_num;++i) {
      //double prob = 0.0;
      //double n = 0;
      int count = queue_atrocity[i].getCount();
      //double allp = ((pa[i][0] + pa[i][1]) == 0)?0.0:(double) (pab[i][0] + pab[i][1]) / (double)(pa[i][0] + pa[i][1]);
      if (count == 0) {
        if (pa[i][0] > 0) {
          ans[i] = (double)pab[i][0] / (double)pa[i][0];
        }
      } else {
        double w = tanh((double)(count) / 13.0);
        //cout<<w * 2 + 1.0<<' ';cin.get();
        if (pa[i][1] > 0) {
          ans[i] = (double)pab[i][1] / (double)pa[i][1];// * (w * 9 + 1.0);
        } 
        ans[i] = ans[i] + w * (1 - ans[i]);
      }
      ans[i] *= factor;//adjust by a factor
    }
    return ans;
  }
};
class CPredictor6
{//same as CPredictor5, but use decay factor to weight more on recent information, 
  static const int kRegionNum = 3671;
  static const int kWindowSize = 700;
  static const int kEventType = 20;
  static const int kReserveDay = 30;//di + 1,di + 30
  int kBackDay;
  double w;
  double ow;
  vector<SRoll> queue_atrocity_big;
  vector<SRoll> queue_atrocity_small;
  vector<SRoll> queue_atrocity;//this is for prediction
  deque<vector<bool> > qab[2];
  deque<vector<bool> > qa[2];
  double pab[kRegionNum][2];
  double pa[kRegionNum][2];
  vector<vector<bool> > cur_sociopolitical;
  vector<bool> cur_atrocities;
  double factor;
public:
  CPredictor6(int backday, double f):factor(f)
  {
    kBackDay = backday;
    for (int i = 0;i < kRegionNum;++i) {
      pab[i][0] = pa[i][0] = 0;
      pab[i][1] = pa[i][1] = 0;
    }
    queue_atrocity_big = vector<SRoll>(kRegionNum, SRoll(kReserveDay + kBackDay));
    queue_atrocity_small = vector<SRoll>(kRegionNum, SRoll(kReserveDay));
    queue_atrocity = vector<SRoll>(kRegionNum, SRoll(kBackDay));
    w = 0.998;
    ow = 1;
    for (int i = 0;i < kWindowSize;++i) {
      ow *= w;
    }
  }
  void readRegion(vector<string> &data)
  {
  }
  vector<bool> readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
    }

    vector<bool> tmp_aid(kRegionNum, false);
    for (int i = 0;i < region_id.size();++i) {
      tmp_aid[region_id[i]] = true;
    }
    return tmp_aid;
  }
  void updateRollingArray()
  {
    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity[i].push(cur_atrocities[i]);
      queue_atrocity_big[i].push(cur_atrocities[i]);
      queue_atrocity_small[i].push(cur_atrocities[i]);
    }
    if (queue_atrocity_big[0].isFull()) {//if is full
      vector<bool> tmp_atrocity(kRegionNum, false);
      for (int i = 0;i < kRegionNum;++i) {
        //if (queue_atrocity[i].getCount() > 0) {
        if (queue_atrocity_small[i].getCount() > 0) {
          tmp_atrocity[i] = true;
        }
      }
      for (int j = 0;j <= 1;++j) {
        vector<bool> tmp_qab(kRegionNum, false);
        vector<bool> tmp_qa(kRegionNum, false);
        for (int i = 0;i < kRegionNum;++i) {
          int event_num = queue_atrocity_big[i].getCount() - 
            queue_atrocity_small[i].getCount();
          int type = (event_num > 0)?1:0;
          if (type == j) {
            tmp_qa[i] = true;
            if (tmp_atrocity[i]) {
              tmp_qab[i] = true;
            }
          }
        }
        for (int i = 0;i < kRegionNum;++i) {
          if (tmp_qab[i]) {
            //pab[i][j] += 1;
            pab[i][j] = w * pab[i][j] + 1.0;
          } else {
            pab[i][j] = w * pab[i][j];
          }
          if (tmp_qa[i]) {
            //pa[i][j] += 1;
            pa[i][j] = w * pa[i][j] + 1.0;
          } else {
            pa[i][j] = w * pa[i][j];
          }
        }
        qab[j].push_back(tmp_qab);
        qa[j].push_back(tmp_qa);
        if (qab[j].size() > kWindowSize) {
          vector<bool> &tv = qab[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              //pab[i][j] -= 1;
              pab[i][j] = pab[i][j] - ow;
            }
          }
          qab[j].pop_front();
        }
        if (qa[j].size() > kWindowSize) {
          vector<bool> &tv = qa[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              //pa[i][j] -= 1;
              pa[i][j] = pa[i][j] - ow;
            }
          }
          qa[j].pop_front();
        }
      }
    }
  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      cur_atrocities.clear();
      cur_atrocities = readAtrocities_data(dayID, data);
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      //cur_sociopolitical.clear();
      //cur_sociopolitical = readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    for (int i = 0;i < region_num;++i) {
      //double prob = 0.0;
      //double n = 0;
      int count = queue_atrocity[i].getCount();
      //double allp = ((pa[i][0] + pa[i][1]) == 0)?0.0:(double) (pab[i][0] + pab[i][1]) / (double)(pa[i][0] + pa[i][1]);
      if (count == 0) {
        if (pa[i][0] > 0) {
          ans[i] = (double)pab[i][0] / (double)pa[i][0];
        }
      } else {
        double w = tanh((double)(count) / 13.0);
        //cout<<w * 2 + 1.0<<' ';cin.get();
        if (pa[i][1] > 0) {
          ans[i] = (double)pab[i][1] / (double)pa[i][1];// * (w * 9 + 1.0);
        } 
        ans[i] = ans[i] + w * (1 - ans[i]);
      }
      ans[i] *= factor;
    }
    return ans;
  }
};
class CPredictor7
{//condition on whether atrocity happen in a country for pass kBackDay
  static const int kRegionNum = 3671;
  static const int kCounryNum = 254;
  static const int kWindowSize = 500;
  static const int kEventType = 20;
  static const int kReserveDay = 30;//di + 1,di + 30
  int kBackDay;
  int kLen;
  vector<SRoll> queue_atrocity_big;//country
  vector<SRoll> queue_atrocity_small;//country
  vector<SRoll> queue_atrocity;//this is for prediction
  vector<SRoll> queue_atrocity_len;//

  deque<vector<bool> > qab[2];
  deque<vector<bool> > qa[2];
  int pab[kRegionNum][2];
  int pa[kRegionNum][2];
  vector<bool> cur_country_atrocities;//in country
  vector<bool> cur_region_atrocities;//in region
  double factor;
public:
  CPredictor7(int backday, int len, double f)
  {
    kBackDay = backday;
    kLen = len;
    factor = f;
    for (int i = 0;i < kRegionNum;++i) {
      pab[i][0] = pa[i][0] = 0;
      pab[i][1] = pa[i][1] = 0;
    }
    queue_atrocity_big = vector<SRoll>(kCounryNum, SRoll(kReserveDay + kBackDay));
    queue_atrocity_small = vector<SRoll>(kCounryNum, SRoll(kReserveDay));
    queue_atrocity = vector<SRoll>(kCounryNum, SRoll(kBackDay));//for prediction
    queue_atrocity_len = vector<SRoll>(kRegionNum, SRoll(kLen));//for prediction

  }
  ~CPredictor7()
  {
  }
  map<int, int> cr_map;//country region map
  void readRegion(vector<string> &data)
  {
    vector<int> reg_start_point;
    for (int i = 0;i < data.size();) {
      reg_start_point.push_back(i);
      //cout<<reg_start_point.size() - 1<<' '<<data[i + 1]<<endl;cin.get();
      istringstream sin1(data[i]);
      istringstream sin2(data[i + 1]);
      int cid, rid;
      sin1>>cid;//assert(cid >= 0 && cid < kCounryNum);
      sin2>>rid;//assert(rid >= 0 && rid < kRegionNum);
      cr_map[rid] = cid;
      i += 2;
      while (i < data.size() && (data[i][0] == 'o' || data[i][0] == 'i') ) {
        i += 2;
      }
    }
  }
  void readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    vector<int> country_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
      country_id.push_back(cid);
    }
    vector<bool> tmp_rid(kRegionNum, false);
    cur_region_atrocities.clear();
    cur_country_atrocities.clear();
    cur_region_atrocities = vector<bool>(kRegionNum, false); 
    cur_country_atrocities = vector<bool>(kCounryNum, false); 
    for (int i = 0;i < region_id.size();++i) {
      cur_region_atrocities[region_id[i]] = true;
    }
    for (int i = 0;i < country_id.size();++i) {
      cur_country_atrocities[country_id[i]] = true;
    }
  }
  void updateRollingArray()
  {
    for (int i = 0;i < kCounryNum;++i) {
      queue_atrocity[i].push(cur_country_atrocities[i]);
      queue_atrocity_big[i].push(cur_country_atrocities[i]);
      queue_atrocity_small[i].push(cur_country_atrocities[i]);
    }

    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity_len[i].push(cur_region_atrocities[i]);
    }
    if (queue_atrocity_big[0].isFull()) {//if is full
      vector<bool> tmp_atrocity(kRegionNum, false);
      for (int i = 0;i < kRegionNum;++i) {
        if (queue_atrocity_len[i].getCount() > 0) {
          tmp_atrocity[i] = true;
        }
      }
      for (int j = 0;j <= 1;++j) {
        vector<bool> tmp_qab(kRegionNum, false);
        vector<bool> tmp_qa(kRegionNum, false);
        for (int i = 0;i < kRegionNum;++i) {
          int cid = cr_map[i];
          int event_num = queue_atrocity_big[cid].getCount() - 
            queue_atrocity_small[cid].getCount();
          int type = (event_num > 0)?1:0;
          if (type == j) {
            tmp_qa[i] = true;
            if (tmp_atrocity[i]) {
              tmp_qab[i] = true;
            }
          }
        }
        for (int i = 0;i < kRegionNum;++i) {
          if (tmp_qab[i]) {
            pab[i][j] += 1;
          }
          if (tmp_qa[i]) {
            pa[i][j] += 1;
          }
        }
        qab[j].push_back(tmp_qab);
        qa[j].push_back(tmp_qa);
        if (qab[j].size() > kWindowSize) {
          vector<bool> &tv = qab[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pab[i][j] -= 1;
            }
          }
          qab[j].pop_front();
        }
        if (qa[j].size() > kWindowSize) {
          vector<bool> &tv = qa[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pa[i][j] -= 1;
            }
          }
          qa[j].pop_front();
        }
      }
    }

  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      readAtrocities_data(dayID, data);
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      //cur_sociopolitical.clear();
      //cur_sociopolitical = readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    for (int i = 0;i < region_num;++i) {
      //double prob = 0.0;
      //double n = 0;
      int cid = cr_map[i];
      int count = queue_atrocity[cid].getCount();
      double allp = ((pa[i][0] + pa[i][1]) == 0)?0.0:(double) (pab[i][0] + pab[i][1]) / (double)(pa[i][0] + pa[i][1]);
      if (count == 0) {
        if (pa[i][0] > 10) {
          ans[i] = (double)pab[i][0] / (double)pa[i][0];
        } else {
          ans[i] = allp;
        }
      } else {
        double w = tanh((double)(count) / 13.0);
        if (pa[i][1] > 10) {
          ans[i] = (double)pab[i][1] / (double)pa[i][1];
        } else {
          ans[i] = allp;
        }
        ans[i] = ans[i] + w * (1 - ans[i]);
      }
      ans[i] *= factor;
    }
    return ans;
  }
};
class CPredictor8
{//condition on whether atrocity happen in a country for pass kBackDay, use a decay factor to weight recent information more
  static const int kRegionNum = 3671;
  static const int kCounryNum = 254;
  static const int kWindowSize = 500;
  static const int kEventType = 20;
  static const int kReserveDay = 30;//di + 1,di + 30
  int kBackDay;
  int kLen;
  double ow;
  double w;
  vector<SRoll> queue_atrocity_big;//country
  vector<SRoll> queue_atrocity_small;//country
  vector<SRoll> queue_atrocity;//this is for prediction
  vector<SRoll> queue_atrocity_len;//

  deque<vector<bool> > qab[2];
  deque<vector<bool> > qa[2];
  double pab[kRegionNum][2];
  double pa[kRegionNum][2];
  vector<bool> cur_country_atrocities;//in country
  vector<bool> cur_region_atrocities;//in region
  double factor;
public:
  CPredictor8(int backday, int len, double f)
  {
    w = 0.999;
    //w = 1.0;
    ow = 1;
    for (int i = 0;i < kWindowSize;++i) {
      ow *= w;
    }
    kBackDay = backday;
    kLen = len;
    factor = f;
    for (int i = 0;i < kRegionNum;++i) {
      pab[i][0] = pa[i][0] = 0.0;
      pab[i][1] = pa[i][1] = 0.0;
    }
    queue_atrocity_big = vector<SRoll>(kCounryNum, SRoll(kReserveDay + kBackDay));
    queue_atrocity_small = vector<SRoll>(kCounryNum, SRoll(kReserveDay));
    queue_atrocity = vector<SRoll>(kCounryNum, SRoll(kBackDay));//for prediction
    queue_atrocity_len = vector<SRoll>(kRegionNum, SRoll(kLen));//for prediction

  }
  ~CPredictor8()
  {
  }
  map<int, int> cr_map;//country region map
  void readRegion(vector<string> &data)
  {
    vector<int> reg_start_point;
    for (int i = 0;i < data.size();) {
      reg_start_point.push_back(i);
      //cout<<reg_start_point.size() - 1<<' '<<data[i + 1]<<endl;cin.get();
      istringstream sin1(data[i]);
      istringstream sin2(data[i + 1]);
      int cid, rid;
      sin1>>cid;//assert(cid >= 0 && cid < kCounryNum);
      sin2>>rid;//assert(rid >= 0 && rid < kRegionNum);
      cr_map[rid] = cid;
      i += 2;
      while (i < data.size() && (data[i][0] == 'o' || data[i][0] == 'i') ) {
        i += 2;
      }
    }
  }
  void readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    vector<int> country_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
      country_id.push_back(cid);
    }
    vector<bool> tmp_rid(kRegionNum, false);
    cur_region_atrocities.clear();
    cur_country_atrocities.clear();
    cur_region_atrocities = vector<bool>(kRegionNum, false); 
    cur_country_atrocities = vector<bool>(kCounryNum, false); 
    for (int i = 0;i < region_id.size();++i) {
      cur_region_atrocities[region_id[i]] = true;
    }
    for (int i = 0;i < country_id.size();++i) {
      cur_country_atrocities[country_id[i]] = true;
    }
  }
  void updateRollingArray()
  {
    for (int i = 0;i < kCounryNum;++i) {
      queue_atrocity[i].push(cur_country_atrocities[i]);
      queue_atrocity_big[i].push(cur_country_atrocities[i]);
      queue_atrocity_small[i].push(cur_country_atrocities[i]);
    }

    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity_len[i].push(cur_region_atrocities[i]);
    }
    if (queue_atrocity_big[0].isFull()) {//if is full
      vector<bool> tmp_atrocity(kRegionNum, false);
      for (int i = 0;i < kRegionNum;++i) {
        if (queue_atrocity_len[i].getCount() > 0) {
          tmp_atrocity[i] = true;
        }
      }
      for (int j = 0;j <= 1;++j) {
        vector<bool> tmp_qab(kRegionNum, false);
        vector<bool> tmp_qa(kRegionNum, false);
        for (int i = 0;i < kRegionNum;++i) {
          int cid = cr_map[i];
          int event_num = queue_atrocity_big[cid].getCount() - 
            queue_atrocity_small[cid].getCount();
          int type = (event_num > 0)?1:0;
          if (type == j) {
            tmp_qa[i] = true;
            if (tmp_atrocity[i]) {
              tmp_qab[i] = true;
            }
          }
        }
        for (int i = 0;i < kRegionNum;++i) {
          pab[i][j] = w * pab[i][j] + ((tmp_qab[i])?1.0:0.0);
          pa[i][j] = w * pa[i][j] + ((tmp_qa[i])?1.0:0.0);
          //if (tmp_qab[i]) {
          //  pab[i][j] += 1;
          //}
          //if (tmp_qa[i]) {
          //  pa[i][j] += 1;
          //}
        }
        qab[j].push_back(tmp_qab);
        qa[j].push_back(tmp_qa);
        if (qab[j].size() > kWindowSize) {
          vector<bool> &tv = qab[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              //pab[i][j] -= 1;
              //cout<<ow<<' ';
              pab[i][j] -= (ow * 1.0);
            }
          }
          qab[j].pop_front();
        }
        if (qa[j].size() > kWindowSize) {
          vector<bool> &tv = qa[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pa[i][j] -= ow * 1.0;
              //pa[i][j] -= 1;
            }
          }
          qa[j].pop_front();
        }
      }
    }

  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      readAtrocities_data(dayID, data);
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      //cur_sociopolitical.clear();
      //cur_sociopolitical = readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    for (int i = 0;i < region_num;++i) {
      //double prob = 0.0;
      //double n = 0;
      int cid = cr_map[i];
      int count = queue_atrocity[cid].getCount();
      double allp = ((pa[i][0] + pa[i][1]) == 0)?0.0:(double) (pab[i][0] + pab[i][1]) / (double)(pa[i][0] + pa[i][1]);
      if (count == 0) {
        if (pa[i][0] > 10) {
          ans[i] = (double)pab[i][0] / (double)pa[i][0];
        } else {
          ans[i] = allp;
        }
      } else {
        double w = tanh((double)(count) / 13.0);
        if (pa[i][1] > 10) {
          ans[i] = (double)pab[i][1] / (double)pa[i][1];
        } else {
          ans[i] = allp;
        }
        ans[i] = ans[i] + w * (1 - ans[i]);
      }
      ans[i] *= factor;
    }
    return ans;
  }
};
class CPredictor9
{//condition on which organization do something in previous day
  static const int kRegionNum = 3671;
  static const int kWindowSize = 500;
  static const int kEventType = 20;
  static const int kReserveDay = 30;//di + 1,di + 30
  static const int kBackDay = 10;
  static const int kBigDay = 50;//
  static const int kSize1 = 10;//
  static const int kSize2 = 20;//30 - 180
  static const int kSize3 = 30;//30 - 180
  static const int kSize4 = 40;//30 - 180
  static const int kSize5 = 50;//30 - 180
  static const int kSize6 = 60;//30 - 180
  static const int kMax = 30000;//suppose the bigget id is kMax - 1
  deque<vector<map<int, int> > > hist_sociopolitical_group;//for calculate statistic
  deque<vector<map<int, int> > > cur_sociopolitical_group;//for prediction
  vector<map<int, int> > old_sociopolitical_group;//kReserveDay's ago's information
  vector<bool> cur_atrocities;
  vector<SRoll> queue_atrocity;
  vector<SRoll> queue_big_atrocity1;//kReserveDay + kBigDay
  vector<SRoll> queue_predict_atrocity1;//kBigDay
  vector<SRoll> queue_predict_atrocity2;//kBigDay
  vector<SRoll> queue_predict_atrocity3;//kBigDay
  vector<SRoll> queue_predict_atrocity4;//kBigDay
  vector<SRoll> queue_predict_atrocity5;//kBigDay
  vector<SRoll> queue_predict_atrocity6;//kBigDay
  vector<double> c1;
  vector<double> c2;

public:
  CPredictor9()
  {
    c1 = vector<double> (kMax, 0);
    c2 = vector<double> (kMax, 0);
    queue_atrocity = vector<SRoll>(kRegionNum, SRoll(kReserveDay));
    queue_big_atrocity1 = vector<SRoll>(kRegionNum, SRoll(kReserveDay + kBigDay));
    queue_predict_atrocity1 = vector<SRoll>(kRegionNum, SRoll(kSize1));
    queue_predict_atrocity2 = vector<SRoll>(kRegionNum, SRoll(kSize2));
    queue_predict_atrocity3 = vector<SRoll>(kRegionNum, SRoll(kSize3));
    queue_predict_atrocity4 = vector<SRoll>(kRegionNum, SRoll(kSize4));
    queue_predict_atrocity5 = vector<SRoll>(kRegionNum, SRoll(kSize5));
    queue_predict_atrocity6 = vector<SRoll>(kRegionNum, SRoll(kSize6));
  }

  void readSociopolitical_data(vector <string> &data)
  {
    vector<map<int, int> > event(kRegionNum, map<int, int>() );
    for (size_t i = 0;i < data.size();++i) {
      int start = 0;
      int end = 0;
      int count = 0;
      int len = data[i].size();
      int rid = -1;//region id
      int pid1 = -1;
      int pid2 = -1;
      //int find = 0;
      for (int j = 0;j < len;++j) {
        if (data[i][j] == ' ' || j == len - 1) {
          count++;//
          if (count == 1) {
            if (data[i][0] == '_') continue;
            int t= 0;
            for (int s = start;s < end;++s) {
              t= t * 10 + data[i][s] - '0';
            }
            pid1 = t;
          } else if(count == 6) {
            if (data[i][start] == '_') continue;
            int t= 0;
            for (int s = start;s < end;++s) {
              t= t * 10 + data[i][s] - '0';
            }
            pid2 = t;
          } else if (count == 15) {
            if (data[i][start] == '_') continue;
            int t= 0;
            for (int s = start;s < end;++s) {
              t= t * 10 + data[i][s] - '0';
            }
            rid = t;
          } else if (count == 11) {
            //if (data[i][start] == 'a' || data[i][start] == 'h' ||data[i][start] == 'n') find = 1;//h very bad, i 3 j5 k2 lm4 n13
          }
          start = j + 1;
          end = j + 1;
        } else {
          end ++;
        }
      }
      if (rid >= 0) {
        if (pid1 >= 0) {
          if (event[rid].count(pid1) == 0) {
            event[rid][pid1] = 1;
          } else {
            event[rid][pid1] = 1;
          }
          if (pid2 >= 0) {
            if (event[rid].count(pid2) == 0) {
              event[rid][pid2] = 1;
            } else {
              event[rid][pid2] = 1;
            }
          }
        }
      }
    }
    old_sociopolitical_group.clear();
    if (hist_sociopolitical_group.size() >= kReserveDay) {
      old_sociopolitical_group = hist_sociopolitical_group.front();
      hist_sociopolitical_group.pop_front();
      hist_sociopolitical_group.push_back(event);
    } else {
      hist_sociopolitical_group.push_back(event);
    }
    if (cur_sociopolitical_group.size() >= kBackDay) {
      cur_sociopolitical_group.pop_front();
      cur_sociopolitical_group.push_back(event);
    } else {
      cur_sociopolitical_group.push_back(event);
    }
  }
  vector<bool> readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
    }
    vector<bool> tmp_aid(kRegionNum, false);
    for (int i = 0;i < region_id.size();++i) {
      tmp_aid[region_id[i]] = true;
    }
    return tmp_aid;
  }
  //void readRegion(vector<string> &data)
  //{
  //}
  void updateRollingArray()
  {
    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity[i].push(cur_atrocities[i]);
      queue_big_atrocity1[i].push(cur_atrocities[i]);
      queue_predict_atrocity1[i].push(cur_atrocities[i]);
      queue_predict_atrocity2[i].push(cur_atrocities[i]);
      queue_predict_atrocity3[i].push(cur_atrocities[i]);
      queue_predict_atrocity4[i].push(cur_atrocities[i]);
      queue_predict_atrocity5[i].push(cur_atrocities[i]);
      queue_predict_atrocity6[i].push(cur_atrocities[i]);
    }
    if (old_sociopolitical_group.size() > 0) {
      vector<bool> mark(kRegionNum, false);
      for (int i = 0;i < kRegionNum;++i) {
        if (queue_atrocity[i].getCount() > 0) {
          mark[i] = true;
        }
      }
      for (int i = 0;i < kRegionNum;++i) {
        map<int, int>::iterator it;
        int oldnum = queue_big_atrocity1[i].getCount() - queue_atrocity[i].getCount();
        if (oldnum > 0) continue;
        for (it = old_sociopolitical_group[i].begin();
            it != old_sociopolitical_group[i].end();++it) {
          int gid = (*it).first;//group id;
          int num = (*it).second;
          if (gid >= kMax || gid < 0) continue;
          c1[gid] += (double) num;
          if (mark[i]) {
            c2[gid] += (double) num;
          }
        }
      }
    }
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    vector<pair<double, int> > prob;prob.reserve(5000);
    for (int k = 0;k < kMax;++k) {
      if (c1[k] == 0) continue;
      double tp = (double)c2[k] / (double)c1[k];
      //if (tp <= 0.2 || c1[k] > 40) continue;
      if (tp <= 0.0002 || c1[k] <= 10) continue;
      prob.push_back(make_pair(tp, k));
    }

    int topnum = 1000;
    if (topnum > prob.size()) topnum = prob.size();
    partial_sort(prob.begin(), prob.begin() + topnum, prob.end(), greater<pair<double, int> > ());
    map<int, int> bad_group;
    for (int i = 0;i < topnum;++i) {
      bad_group[prob[i].second] = 1;
    }
    for (int i = 0;i < region_num;++i) {
      double prob = 0.0;
      double n = 0;
      double h = 1.0;
      if (queue_predict_atrocity6[i].getCount() == 0) {
        h = 1.0;
      } else if(queue_predict_atrocity5[i].getCount() == 0) {
        h = 0.9;
      } else if(queue_predict_atrocity4[i].getCount() == 0) {
        h = 0.8;
      } else if(queue_predict_atrocity3[i].getCount() == 0) {
        h = 0.7;
      } else if(queue_predict_atrocity2[i].getCount() == 0) {
        h = 0.6;
      } else if(queue_predict_atrocity1[i].getCount() == 0) {
        h = 0.5;
      } else {
        continue;
      }
      deque<vector<map<int, int> > >::iterator it;
      double w = 1.0;
      for (it = cur_sociopolitical_group.begin();
          it != cur_sociopolitical_group.end();++it) {
        map<int, int>::iterator ip;
        for (ip = (*it)[i].begin();ip != (*it)[i].end();++ip) {
          int gid = (*ip).first;
          //int num = (*ip).second;
          if (bad_group.count(gid) > 0) {//is bad group;
            double tmp_prob = 0.0;
            if (c1[gid] > 0) {
              tmp_prob = (double)c2[gid] / (double)c1[gid];
            }
            double tw = (c1[gid] > 100)?2.0:1.0;
            prob += w * tmp_prob * tw;;
            n += w * tw;
          }
        }
        w *= 0.8;
      }
      if (n > 0) ans[i] = prob / n / (double)kBackDay * h;
      //if (n > 0) ans[i] = prob / n / (double)kBackDay;
    }
    return ans;

  }

  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      cur_atrocities.clear();
      cur_atrocities = readAtrocities_data(dayID, data);
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      //readRegion(data);
    } 
    return 0;
  }

};
class CPredictor10
{//condition on whether an atrocity happen in a area and its neighor in last few days
  static const int kRegionNum = 3671;
  //static const int kWindowSize = 500;
  static const int kEventType = 20;
  static const int kReserveDay = 30;//di + 1,di + 30
  //static const int kBackDay = 20;
  static const int kOffset = 90;
  int kBackDay;
  int kWindowSize;
  int kLen;
  vector<SRoll> queue_atrocity_big;
  vector<SRoll> queue_atrocity_small;
  vector<SRoll> queue_atrocity;//this is for prediction
  vector<SRoll> queue_atrocity_tmp;//this is for prediction
  deque<vector<bool> > qab[2];
  deque<vector<bool> > qa[2];
  vector<bool> cur_atrocities;


  int pab[kRegionNum][2];
  int pa[kRegionNum][2];
  vector<int> earth[2 * kOffset + 1][2 * kOffset + 1];
  vector<int> neighbor[kRegionNum];
  double factor;
public:
  CPredictor10(int kbd, int len, double f)
  {
    kBackDay = kbd;
    kLen = len;
    factor = f;
    kWindowSize = 500;
    for (int i = 0;i < kRegionNum;++i) {
      pab[i][0] = pa[i][0] = 0;
      pab[i][1] = pa[i][1] = 0;
    }
    queue_atrocity_big = vector<SRoll>(kRegionNum, SRoll(kReserveDay + kBackDay));
    queue_atrocity_small = vector<SRoll>(kRegionNum, SRoll(kReserveDay));
    queue_atrocity = vector<SRoll>(kRegionNum, SRoll(kBackDay));
    queue_atrocity_tmp = vector<SRoll>(kRegionNum, SRoll(kLen));
  }
  void readRegion(vector<string> &data)
  {
    vector<vector<int> > region_neighbor;
    int p;
    int rid = 0;;
    for (p = 0;p < data.size();) {
      //cout<<p<<' '<<endl;fflush(stdout);
      istringstream sin2(data[p + 1]);
      int tmp_rid;
      sin2>>tmp_rid;//assert(rid == tmp_rid);
      p++;
      p++;
      while(p < data.size() && (data[p][0] == 'o' || data[p][0] == 'i')) {
        p++;
        for (int i = 0;i < data[p].size();++i) {
          if(data[p][i] == ',') {
            data[p][i] = ' ';
          }
        }
        istringstream sin(data[p]);
        vector<double> edge_x;
        vector<double> edge_y;
        double tx, ty;
        //cout<<data[p]<<endl;
        while(sin>>tx>>ty) {
          edge_x.push_back(tx);
          edge_y.push_back(ty);
        }
        for (int i = 0;i < edge_x.size();++i) {
          int mx = floor(edge_x[i] / 2 + kOffset);//assert(mx >= 0 && mx < 2 * kOffset + 1);
          int my = floor(edge_y[i] + kOffset);//assert(my >= 0 && my < 2 * kOffset + 1);
          vector<int>::iterator it = find(earth[mx][my].begin(), earth[mx][my].end(), rid);
          if (it == earth[mx][my].end()) {
            earth[mx][my].push_back(rid);
          }
        }
        p++;
      }
      rid++;
    }
    int M = 2 * kOffset + 1;
    //assert(rid == kRegionNum);cout<<"pass number check"<<endl;
    int dx[] = {1, 0, -1, 0};
    int dy[] = {0, 1, 0, -1};
    for (int i = 0;i < M;++i) {
      for (int j = 0;j < M;++j) {
        for (int k = 0;k < 4;++k) {
          int cx = (i + dx[k] + M) % M;
          int cy = (j + dy[k] + M) % M;
          for (int s = 0;s < earth[cx][cy].size();++s) {
            for (int t = 0;t < earth[i][j].size();++t) {
              int tn = earth[cx][cy][s];
              int cn = earth[i][j][t];
              vector<int>::iterator it = find(neighbor[cn].begin(), neighbor[cn].end(), tn);
              if (it == neighbor[cn].end()) {
                neighbor[cn].push_back(tn);
              }
            }
          }
        }
      }
    }
  }
  vector<bool> readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
    }

    vector<bool> tmp_aid(kRegionNum, false);
    for (int i = 0;i < region_id.size();++i) {
      tmp_aid[region_id[i]] = true;
    }
    return tmp_aid;
  }
  void updateRollingArray()
  {
    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity[i].push(cur_atrocities[i]);
      queue_atrocity_tmp[i].push(cur_atrocities[i]);
      queue_atrocity_big[i].push(cur_atrocities[i]);
      queue_atrocity_small[i].push(cur_atrocities[i]);
    }
    if (queue_atrocity_big[0].isFull()) {//if is full
      vector<bool> tmp_atrocity(kRegionNum, false);//atrocity happen in a place
      vector<bool> neighbor_atrocity(kRegionNum, false);//in old day whether neighbor has atrocity;
      for (int i = 0;i < kRegionNum;++i) {
        if (queue_atrocity_tmp[i].getCount() > 0) {
          tmp_atrocity[i] = true;
          //for (int j = 0;j < neighbor[i].size();++j) {
          //  tmp_atrocity[neighbor[i]] = true;
          //}
        }
      }
      for (int i = 0;i < kRegionNum;++i) {
        int event_num = queue_atrocity_big[i].getCount() - 
          queue_atrocity_small[i].getCount();
        if (event_num > 0) {
          for (int j = 0;j < neighbor[i].size();++j) {
            neighbor_atrocity[neighbor[i][j]] = true;
          }
        }
      }
      for (int j = 0;j <= 1;++j) {
        vector<bool> tmp_qab(kRegionNum, false);
        vector<bool> tmp_qa(kRegionNum, false);
        for (int i = 0;i < kRegionNum;++i) {
          int type = (neighbor_atrocity[i])?1:0;
          if (type == j) {
            tmp_qa[i] = true;
            if (tmp_atrocity[i]) {
              tmp_qab[i] = true;
            }
          }
        }
        for (int i = 0;i < kRegionNum;++i) {
          if (tmp_qab[i]) {
            pab[i][j] += 1;
          }
          if (tmp_qa[i]) {
            pa[i][j] += 1;
          }
        }
        qab[j].push_back(tmp_qab);
        qa[j].push_back(tmp_qa);
        if (qab[j].size() > kWindowSize) {
          vector<bool> &tv = qab[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pab[i][j] -= 1;
            }
          }
          qab[j].pop_front();
        }
        if (qa[j].size() > kWindowSize) {
          vector<bool> &tv = qa[j].front();
          for (int i = 0;i < kRegionNum;++i) {
            if (tv[i]) {
              pa[i][j] -= 1;
            }
          }
          qa[j].pop_front();
        }
      }
    }

  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      cur_atrocities.clear();
      cur_atrocities = readAtrocities_data(dayID, data);
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      //readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    vector<int> neighbor_atrocity(kRegionNum, 0);//in old day whether neighbor has atrocity;
    for (int i = 0;i < kRegionNum;++i) {
      if (queue_atrocity[i].getCount() > 0) {
        neighbor_atrocity[i] += 1;
        for (int j = 0;j < neighbor[i].size();++j) {
          neighbor_atrocity[neighbor[i][j]] += 1;
        }
      }
    }
    for (int i = 0;i < region_num;++i) {
      //double prob = 0.0;
      //double n = 0;
      //int count = queue_atrocity[i].getCount();
      double allp = ((pa[i][0] + pa[i][1]) == 0)?0.0:(double) (pab[i][0] + pab[i][1]) / (double)(pa[i][0] + pa[i][1]);
      if (neighbor_atrocity[i] == 0) {
      //if (neighbor_atrocity[i]) {
        if (pa[i][0] > 0) {
          ans[i] = (double)pab[i][0] / (double)pa[i][0];
        } else {
          ans[i] = allp;
        }
      } else {
        double w = tanh((double)(neighbor_atrocity[i]) / 13.0);
        if (pa[i][1] > 0) {
          ans[i] = (double)pab[i][1] / (double)pa[i][1];
        } else {
          ans[i] = allp;
        }
        ans[i] = ans[i] + w * (1 - ans[i]);
      }
      ans[i] *= factor;
      //ans[i] /= 3.5;//4 is just a adjusted factor, need to tune
    }
    return ans;
  }

};


class CPredictor13
{//calculate these three things:1. which organizetion cause atrocity most frequently
 //                             2. which two countries' interaction causes atrocity most frequently
 //                             3. which two neighboring regions' interaction causes atrocity most frequently
 //                             4. one country do a special event 'a' - 't'
 //can get a positive result, but not very good yet, so do not use when combining, just put it here. 
  static const int kOffset = 90;
  static const int kCountryNum = 254;
  static const int kRegionNum = 3671;
  static const int kBackDay = 7;
  static const int kReserveDay = 30;//di + 1,di + 30
  static const int kMax = 70000;//suppose the bigget id is kMax - 1
  static const int kBigDay = 50;//
  static const int kType = 4;
  vector<double> c1[kType];//0 for bad organization, 1 for bad country relationship, 2 for bad neighbor relationship
  vector<double> c2[kType];
  deque<vector<map<int, int> > > hist_sociopolitical_group[kType];//for calculate statistic
  deque<vector<map<int, int> > > cur_sociopolitical_group[kType];//for prediction
  vector<map<int, int> > old_sociopolitical_group[kType];//kReserveDay's ago's information
  vector<int> earth[2 * kOffset + 1][2 * kOffset + 1];
  vector<int> neighbor[kRegionNum];//put each region's neighbor here
  vector<bool> cur_atrocities;
  vector<bool> old_atrocities;
  deque<vector<bool> > atro_queue;
  map<pair<int, int>, int> nei_mp;
  vector<int> latest;
  vector<double> weight;
  vector<SRoll> queue_atrocity;
  vector<SRoll> queue_big_atrocity1;//kReserveDay + kBigDay
public:
  CPredictor13()
  {
    for (int s = 0;s < kType;++s) {
      c1[s] = vector<double> (kMax, 0);
      c2[s] = vector<double> (kMax, 0);
    }
    latest = vector<int>(kRegionNum, -1);
    weight = vector<double>(kRegionNum, 1.0);
    queue_atrocity = vector<SRoll>(kRegionNum, SRoll(kReserveDay));
    queue_big_atrocity1 = vector<SRoll>(kRegionNum, SRoll(kReserveDay + kBigDay));
  }
  void readRegion(vector<string> &data)
  {
    vector<vector<int> > region_neighbor;
    int p;
    int rid = 0;;
    for (p = 0;p < data.size();) {
      istringstream sin2(data[p + 1]);
      int tmp_rid;
      sin2>>tmp_rid;//assert(rid == tmp_rid);
      p++;
      p++;
      while(p < data.size() && (data[p][0] == 'o' || data[p][0] == 'i')) {
        p++;
        for (int i = 0;i < data[p].size();++i) {
          if(data[p][i] == ',') {
            data[p][i] = ' ';
          }
        }
        istringstream sin(data[p]);
        vector<double> edge_x;
        vector<double> edge_y;
        double tx, ty;
        //cout<<data[p]<<endl;
        while(sin>>tx>>ty) {
          edge_x.push_back(tx);
          edge_y.push_back(ty);
        }
        for (int i = 0;i < edge_x.size();++i) {
          int mx = floor(edge_x[i] / 2 + kOffset);//assert(mx >= 0 && mx < 2 * kOffset + 1);
          int my = floor(edge_y[i] + kOffset);//assert(my >= 0 && my < 2 * kOffset + 1);
          vector<int>::iterator it = find(earth[mx][my].begin(), earth[mx][my].end(), rid);
          if (it == earth[mx][my].end()) {
            earth[mx][my].push_back(rid);
          }
        }
        p++;
      }
      rid++;
    }
    int M = 2 * kOffset + 1;
    int dx[] = {1, 0, -1, 0};
    int dy[] = {0, 1, 0, -1};
    for (int i = 0;i < M;++i) {
      for (int j = 0;j < M;++j) {
        for (int k = 0;k < 4;++k) {
          int cx = (i + dx[k] + M) % M;
          int cy = (j + dy[k] + M) % M;
          for (int s = 0;s < earth[cx][cy].size();++s) {
            for (int t = 0;t < earth[i][j].size();++t) {
              int tn = earth[cx][cy][s];
              int cn = earth[i][j][t];
              vector<int>::iterator it = find(neighbor[cn].begin(), neighbor[cn].end(), tn);
              if (it == neighbor[cn].end()) {
                neighbor[cn].push_back(tn);
              }
            }
          }
        }
      }
    }

    for (int i = 0;i < kRegionNum;++i) {
      for (int j = 0;j < neighbor[i].size();++j) {
        int rid1 = i;
        int rid2 = neighbor[i][j];
        if (rid1 == rid2) continue;
        if (rid2 > rid1) swap(rid2, rid1);
        if (nei_mp.count(make_pair(rid1, rid2)) == 0) {
          int pin = nei_mp.size();
          nei_mp[make_pair(rid1, rid2)] = pin;
        }
      }
    }
  }
  int getID(int id1, int id2)
  {
    return id1 * kCountryNum + id2;
  }
  void readSociopolitical_data(vector <string> &data)
  {
    vector<map<int, int> > event[kType];
    for (int i = 0;i < kType;++i) {
      event[i] = vector<map<int, int> >(kRegionNum, map<int, int>() );
    }
    //vector<map<int, int> > event1(kRegionNum, map<int, int>() );
    //vector<map<int, int> > event2(kRegionNum, map<int, int>() );
    for (size_t i = 0;i < data.size();++i) {
      int start = 0;
      int end = 0;
      int count = 0;
      int len = data[i].size();
      int rid = -1;//region id
      int pid1 = -1;
      int pid2 = -1;
      int rid1 = -1;
      int rid2 = -1;
      int cid1 = -1;
      int cid2 = -1;
      //int find = 0;
      int event_weight = 1;
      int event_id = -1;
      for (int j = 0;j < len;++j) {
        if (data[i][j] == ' ' || j == len - 1) {
          count++;//
          if (data[i][start] != '_') {// continue;
            if (count == 1) {
              pid1 = 0;
              for (int s = start;s < end;++s) pid1 = pid1 * 10 + data[i][s] - '0';
            } else if(count == 4) {
              cid1 = 0;
              for (int s = start;s < end;++s) cid1 = cid1 * 10 + data[i][s] - '0';
            } else if(count == 5) {
              rid1 = 0;
              for (int s = start;s < end;++s) rid1 = rid1 * 10 + data[i][s] - '0';
            } else if(count == 6) {
              pid2 = 0;
              for (int s = start;s < end;++s) pid2 = pid2 * 10 + data[i][s] - '0';
            } else if(count == 9) {
              cid2 = 0;
              for (int s = start;s < end;++s) cid2 = cid2 * 10 + data[i][s] - '0';
            } else if(count == 10) {
              rid2 = 0;
              for (int s = start;s < end;++s) rid2 = rid2 * 10 + data[i][s] - '0';
            } else if (count == 15) {
              rid = 0;
              for (int s = start;s < end;++s) rid = rid * 10 + data[i][s] - '0';
            } else if (count == 11) {

              event_id = data[i][start] - 'a';
              //if (data[i][start] == 't') event_weight = 2;
              //cover = 0;
              //for (int s = start;s < end;++s) cover = cover * 10 + data[i][s] - '0';
              //if (data[i][start] == 'a' || data[i][start] == 'h' ||data[i][start] == 'n') find = 1;//h very bad, i 3 j5 k2 lm4 n13
            }
          }
          start = j + 1;
          end = j + 1;
        } else {
          end ++;
        }
      }
      //if (cover > 50) event_weight = 2;
      //continue;
      if (rid < 0) continue;
      //for 0, different organization
      if (pid1 >= 0) {
        if (event[0][rid].count(pid1) == 0) event[0][rid][pid1] = event_weight;
        else event[0][rid][pid1] += event_weight;
      }
      if (pid2 >= 0) {
        if (event[0][rid].count(pid2) == 0) event[0][rid][pid2] = event_weight;
        else event[0][rid][pid2] = event_weight;
      }

      //for 1, different relationship between countries
      if (cid1 >= 0 && cid2 >= 0) {
        if (cid2 > cid1) swap(cid1, cid2); 
        int ID = getID(cid1, cid2);
        if (event[1][rid].count(ID) == 0) event[1][rid][ID] = event_weight;
        else event[1][rid][ID] += event_weight;
      }

      //for 2, bad neighbor relationship
      if (rid1 >= 0 && rid2 >= 0) {
        if (rid2 > rid1) swap(rid2, rid1);
        if (nei_mp.count(make_pair(rid1, rid2)) > 0) {
          int ID = nei_mp[make_pair(rid1, rid2)];
          if (event[2][rid].count(ID) == 0) event[2][rid][ID] = event_weight;
          else event[2][rid][ID] += event_weight;
        }
      }

      if (cid1 >= 0 && cid1 < kCountryNum && event_id >= 0) {
        int pin = cid1 * 20 + event_id;//
        if (event[3][rid].count(pin) == 0) event[3][rid][pin] = event_weight;
        else event[3][rid][pin] += event_weight;
      }
      if (cid2 >= 0 && cid2 < kCountryNum && event_id >= 0) {
        int pin = cid2 * 20 + event_id;//
        if (event[3][rid].count(pin) == 0) event[3][rid][pin] = event_weight;
        else event[3][rid][pin] += event_weight;
      }
    }
    //return;
    for (int i = 0;i < kType;++i) {
      old_sociopolitical_group[i].clear();
      if (hist_sociopolitical_group[i].size() >= kReserveDay) {
        old_sociopolitical_group [i]= hist_sociopolitical_group[i].front();
        hist_sociopolitical_group[i].pop_front();
        hist_sociopolitical_group[i].push_back(event[i]);
      } else {
        hist_sociopolitical_group[i].push_back(event[i]);
      }
      if (cur_sociopolitical_group[i].size() >= kBackDay) {
        cur_sociopolitical_group[i].pop_front();
        cur_sociopolitical_group[i].push_back(event[i]);
      } else {
        cur_sociopolitical_group[i].push_back(event[i]);
      }
    }
  }
  void updateRollingArray()
  {
    //return;
    for (int i = 0;i < kRegionNum;++i) {
      queue_atrocity[i].push(cur_atrocities[i]);
      queue_big_atrocity1[i].push(cur_atrocities[i]);
    }
    if (old_sociopolitical_group[0].size() > 0) {
      vector<bool> mark(kRegionNum, false);
      for (int i = 0;i < kRegionNum;++i) {
        if (queue_atrocity[i].getCount() > 0) {
          mark[i] = true;
        }
      }
      for (int s = 0;s < kType;++s) {
        for (int i = 0;i < kRegionNum;++i) {
          map<int, int>::iterator it;
          int oldnum = queue_big_atrocity1[i].getCount() - queue_atrocity[i].getCount();
          if (oldnum > 0) continue;
          for (it = old_sociopolitical_group[s][i].begin();
              it != old_sociopolitical_group[s][i].end();++it) {
            int gid = (*it).first;//group id;
            int num = (*it).second;
            if (gid >= kMax || gid < 0) continue;
            c1[s][gid] += (double) num * weight[i];//
            if (mark[i]) {
              c2[s][gid] += (double) num * weight[i];
            }
          }
        }
      }
    }
  }
  vector<bool> readAtrocities_data(int dayID, vector <string> &data)
  {
    vector<int> region_id;
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      //int did;
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      region_id.push_back(rid);
    }
    vector<bool> tmp_aid(kRegionNum, false);
    for (int i = 0;i < region_id.size();++i) {
      tmp_aid[region_id[i]] = true;
    }
    return tmp_aid;
  }
  int receiveData(int dataSourceId, int dayID, vector <string> &data)
  {
    if (dataSourceId == 0) {//2
      cur_atrocities.clear();
      cur_atrocities = readAtrocities_data(dayID, data);
      atro_queue.push_back(cur_atrocities);
      old_atrocities.clear();
      if (atro_queue.size() > kReserveDay) {
        old_atrocities = atro_queue.front();
        atro_queue.pop_front();
      }
      if (old_atrocities.size() > 0) {
        for (int i = 0;i < kRegionNum;++i) {
          if (old_atrocities[i]) {
            latest[i] = dayID - 30;
          }
        }
        for (int i = 0;i < kRegionNum;++i) {
          if (latest[i] > 0) {
            weight[i] = tanh(((double)(dayID - 30 - latest[i] + 10) / 180.0));
          }
        }
      }
      updateRollingArray();
    } else if(dataSourceId == 1) {//1
      readSociopolitical_data(data);
    } else if(dataSourceId == 2) {//0
      readRegion(data);
    } 
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {
    int region_num = kRegionNum;
    vector<double> ans(region_num, 0.0);
    vector<double> ans_count(region_num, 0.0);
    //return ans;
    for (int s = 0;s <= 3;++s) {
      vector<pair<double, int> > prob;prob.reserve(50000);
      for (int k = 0;k < kMax;++k) {
        if (c1[s][k] == 0) continue;
        double tp = (double)c2[s][k] / (double)c1[s][k];
        //if (tp <= 0.2 || c1[k] > 40) continue;
        if (tp <= 0.0002 || c1[s][k] <= 10) continue;
        prob.push_back(make_pair(tp, k));
      }

      int topnum = 1000;
      if (topnum > prob.size()) topnum = prob.size();
      partial_sort(prob.begin(), prob.begin() + topnum, prob.end(), greater<pair<double, int> > ());
      map<int, int> bad_group;
      for (int i = 0;i < topnum;++i) {
        bad_group[prob[i].second] = 1;
      }
      for (int i = 0;i < region_num;++i) {
        double prob = 0.0;
        double n = 0;
        deque<vector<map<int, int> > >::iterator it;
        double w = 1.0;
        for (it = cur_sociopolitical_group[s].begin();
            it != cur_sociopolitical_group[s].end();++it) {
          map<int, int>::iterator ip;
          for (ip = (*it)[i].begin();ip != (*it)[i].end();++ip) {
            int gid = (*ip).first;
            //int num = (*ip).second;
            if (bad_group.count(gid) > 0) {//is bad group;
              double tmp_prob = 0.0;
              if (c1[s][gid] > 0) {
                tmp_prob = (double)c2[s][gid] / (double)c1[s][gid];
              }
              prob += w * tmp_prob;
              n += w;
            }
          }
          w *= 0.8;
        }
        //if (n > 0) ans[i] = prob / n / (double)kBackDay * h;
        if (n > 0) {
          //if (prob / n / (double)kBackDay > ans[i]) {
          //  ans[i] = prob / n / (double)kBackDay;
          //}
          ans[i] += prob / n / (double)kBackDay * 2;
          ans_count[i] += 1.0;
        }
      }
    }
    for (int i = 0;i < kRegionNum;++i) {
      //if (ans_count[i] > 0) ans[i] /= ans_count[i];
      ans[i] /= 4.0;
    }
    return ans;

  }

  
};
class MassAtrocityPredictor : public BaseMassAtrocityPredictor
{
  static const int kRegionNum = 3671;
  static const int kCounryNum = 254;
  static const int kReserveDay = 30;
  static const int kFeatureNum = 9;
  vector<CPredictor2> obj2;
  vector<CPredictor3> obj3;
  vector<CPredictor4> obj4;
  vector<CPredictor5> obj5;
  vector<CPredictor6> obj6;
  vector<CPredictor7> obj7;
  vector<CPredictor8> obj8;
  CPredictor9 obj9;
  vector<CPredictor10> obj10;
  CPredictor13 obj13;
  vector<SRoll> queue_atrocity_reg;
public:

  MassAtrocityPredictor() 
  {
    queue_atrocity_reg = vector<SRoll>(kRegionNum, SRoll(500));

    obj2.push_back(CPredictor2(0.68, 700));
    obj2.push_back(CPredictor2(0.85, 600));
    obj2.push_back(CPredictor2(0.85, 500));

    //obj2.push_back(CPredictor2(0.50, 500));
    //obj2.push_back(CPredictor2(0.50, 500));

    obj3.push_back(CPredictor3(0.5, 700));
    obj3.push_back(CPredictor3(0.85, 600));
    obj3.push_back(CPredictor3(0.85, 500));

    //obj3.push_back(CPredictor3(0.5, 700));
    //obj3.push_back(CPredictor3(0.68, 600));

    obj4.push_back(CPredictor4(0.7, 500));
    obj4.push_back(CPredictor4(0.8, 600));
    obj4.push_back(CPredictor4(0.9, 500));

    obj5.push_back(CPredictor5(18, 1.0));
    obj5.push_back(CPredictor5(21, 1.0));
    obj5.push_back(CPredictor5(24, 1.0));
    obj5.push_back(CPredictor5(24, 1.0));
    obj5.push_back(CPredictor5(27, 1.0));
    obj5.push_back(CPredictor5(30, 1.0));

    obj6.push_back(CPredictor6(30, 0.68));
    obj6.push_back(CPredictor6(25, 0.68));
    obj6.push_back(CPredictor6(20, 0.68));
    obj6.push_back(CPredictor6(30, 0.9));
    obj6.push_back(CPredictor6(25, 0.9));
    obj6.push_back(CPredictor6(20, 0.9));

    obj7.push_back(CPredictor7(7, 15, 0.5));
    obj7.push_back(CPredictor7(15, 15, 0.9));
    //obj7.push_back(CPredictor7(15, 15, 0.8));
    obj8.push_back(CPredictor8(8, 15, 0.68));
    obj8.push_back(CPredictor8(15, 15, 0.8));
    //obj8.push_back(CPredictor8(18, 18, 0.5));
    obj10.push_back(CPredictor10(10, 30, 0.3));
    obj10.push_back(CPredictor10(7, 30, 0.4));
    obj10.push_back(CPredictor10(7, 30, 0.5));
  }
  //map<int, int> cr_map;//country region map
  //vector<int> rig_count;
  //void readRegion(vector<string> &data)
  //{
  //  rig_count = vector<int> (kCounryNum, 0);
  //  vector<int> reg_start_point;
  //  for (int i = 0;i < data.size();) {
  //    reg_start_point.push_back(i);
  //    istringstream sin1(data[i]);
  //    istringstream sin2(data[i + 1]);
  //    int cid, rid;
  //    sin1>>cid;//assert(cid >= 0 && cid < kCountryNum);
  //    sin2>>rid;//assert(rid >= 0 && rid < kRegionNum);
  //    cr_map[rid] = cid;
  //    rig_count[cid] += 1;
  //    i += 2;
  //    while (i < data.size() && (data[i][0] == 'o' || data[i][0] == 'i') ) {
  //      i += 2;
  //    }
  //  }
  //}
  void readAtrocities_data(vector <string> &data, vector<int> &reg_num)
  {
    vector<int> tmp_region;
    reg_num = vector<int> (kRegionNum, 0);
    for (int i = 0;i < data.size();++i) {
      istringstream sin(data[i]);
      int cid;
      int rid;
      string pos1, pos2;
      sin>>pos1>>pos2>>cid>>rid;
      reg_num[rid] += 1;
    }
  }
  //TUNE
  int receiveData(int dataSourceId, int dayID, vector <string> data)
  {
    if (dataSourceId == 0) {
      vector<int> reg_num;
      readAtrocities_data(data, reg_num);
      for (int i = 0;i < kRegionNum;++i) {
        queue_atrocity_reg[i].push(reg_num[i]);
      }
    //} else if(dataSourceId == 2) {
    //  readRegion(data);
    }
    for (int i = 0;i < obj2.size();++i) {
      obj2[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj3.size();++i) {
      obj3[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj4.size();++i) {
      obj4[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj5.size();++i) {
      obj5[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj6.size();++i) {
      obj6[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj7.size();++i) {
      obj7[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj8.size();++i) {
      obj8[i].receiveData(dataSourceId, dayID, data);
    }
    for (int i = 0;i < obj10.size();++i) {
      obj10[i].receiveData(dataSourceId, dayID, data);
    }
    obj9.receiveData(dataSourceId, dayID, data);
    obj13.receiveData(dataSourceId, dayID, data);
    return 0;
  }
  vector <double> predictAtrocities(int dayID)
  {/*
    tune the parameter in this way, because different regions have different situation, 
    divide regions in to 3 groups, (1)actrocity happen number<= 1, (2)actrocity happen number<= 5 and >=1 (3)actrocity happen number > 5
    and tune parameters for each groups. 
    and the reason using CPrecitor9 in this way:
    if (ans_all[i] <= ans9[i]) ans_all[i] = ans9[i];
    is that CPredictor9 is like event driven, if some rare event happen(like some special organizations do something), what will happen in the future. 
    After trying find this is the relatively good way to use it. 
   */
    vector<double> ans_all(kRegionNum, 0.0);
    vector<double> ans20 = obj2[0].predictAtrocities(dayID);
    vector<double> ans21 = obj2[1].predictAtrocities(dayID);
    vector<double> ans22 = obj2[2].predictAtrocities(dayID);
    //vector<double> ans23 = obj2[3].predictAtrocities(dayID);
    //vector<double> ans24 = obj2[4].predictAtrocities(dayID);

    vector<double> ans30 = obj3[0].predictAtrocities(dayID);
    vector<double> ans31 = obj3[1].predictAtrocities(dayID);
    vector<double> ans32 = obj3[2].predictAtrocities(dayID);
    //vector<double> ans33 = obj3[3].predictAtrocities(dayID);
    //vector<double> ans34 = obj3[4].predictAtrocities(dayID);

    vector<double> ans40 = obj4[0].predictAtrocities(dayID);
    vector<double> ans41 = obj4[1].predictAtrocities(dayID);
    vector<double> ans42 = obj4[2].predictAtrocities(dayID);
    vector<double> ans50 = obj5[0].predictAtrocities(dayID);
    vector<double> ans51 = obj5[1].predictAtrocities(dayID);
    vector<double> ans52 = obj5[2].predictAtrocities(dayID);
    vector<double> ans53 = obj5[3].predictAtrocities(dayID);
    vector<double> ans54 = obj5[4].predictAtrocities(dayID);
    vector<double> ans55 = obj5[5].predictAtrocities(dayID);
    vector<double> ans60 = obj6[0].predictAtrocities(dayID);
    vector<double> ans61 = obj6[1].predictAtrocities(dayID);
    vector<double> ans62 = obj6[2].predictAtrocities(dayID);
    vector<double> ans63 = obj6[3].predictAtrocities(dayID);
    vector<double> ans64 = obj6[4].predictAtrocities(dayID);
    vector<double> ans65 = obj6[5].predictAtrocities(dayID);
    vector<double> ans70 = obj7[0].predictAtrocities(dayID);
    vector<double> ans71 = obj7[1].predictAtrocities(dayID);
    //vector<double> ans72 = obj7[2].predictAtrocities(dayID);
    vector<double> ans80 = obj8[0].predictAtrocities(dayID);
    vector<double> ans81 = obj8[1].predictAtrocities(dayID);
    //vector<double> ans82 = obj8[2].predictAtrocities(dayID);
    
    vector<double> ans9 = obj9.predictAtrocities(dayID);
    vector<double> ans100 = obj10[0].predictAtrocities(dayID);
    vector<double> ans101 = obj10[1].predictAtrocities(dayID);
    vector<double> ans102 = obj10[2].predictAtrocities(dayID);
    vector<double> ans13 = obj13.predictAtrocities(dayID);

    for (int i = 0;i < kRegionNum;++i) {
      ans9[i] = 0.6 * ans9[i] + 0.4 * ans13[i];
    }
    for (int i = 0;i < kRegionNum;++i) {
      int hc = queue_atrocity_reg[i].getCount();
      if (hc <= 1) {//if less frequent, weight more by country or neighbor
        ans_all[i] = 0.18 * (0.8 * ans20[i] + 0.2 * ans30[i])
          + 0.4 * ans41[i]
          + 0.09 * (0.4 * ans50[0] + 0.3 * ans51[i] + 0.3 * ans52[i])
          + 0.09 * (0.4 * ans60[0] + 0.3 * ans61[i] + 0.3 * ans62[i])
          + 0.07 * (ans70[i])
          + 0.07 * (ans80[i])
          + 0.10 * (ans100[i]);
          //+ 0.05 * (ans13[i]);
        if (ans_all[i] <= ans9[i]) ans_all[i] = ans9[i];
      } else if(hc <= 5) {
        ans_all[i] = 0.4 * (0.8 * ans21[i] + 0.2 * ans31[i])
          + 0.19 * (0.4 * ans50[0] + 0.3 * ans51[i] + 0.3 * ans52[i])
          + 0.19 * (0.4 * ans60[0] + 0.3 * ans61[i] + 0.3 * ans62[i])
          + 0.19 * (ans71[i])
          + 0.01 * (ans81[i])
          + 0.02 * (ans101[i]);
        if (ans_all[i] <= ans9[i]) ans_all[i] = ans9[i];
      } else {
        ans_all[i] = 0.46 * (0.8 * ans22[i] + 0.2 * ans32[i])
          + 0.02 * ans42[i]
          + 0.09 * (0.4 * ans53[0] + 0.3 * ans54[i] + 0.3 * ans55[i])
          + 0.09 * (0.4 * ans63[0] + 0.3 * ans64[i] + 0.3 * ans65[i])
          + 0.09 * (ans71[i])
          + 0.09 * (ans81[i])
          + 0.16 * (ans102[i]);
        if (ans_all[i] <= ans9[i]) ans_all[i] = ans9[i];
      }
      //if (isnan((double)ans_all[i]) || isinf((double)ans_all[i])) ans_all[i] = 0;//if wrong
      if (ans_all[i] < 0) ans_all[i] = 0.0;
      if (ans_all[i] > 1) ans_all[i] = 1.0;
    }
    return ans_all; 
  }
};

#ifdef LOCAL
const int kMaxLine = 1000000;
const int kMaxLine1 = 1000;
map<int, vector<string> > all_atrocities_mp;
double getMean(vector<double> &s)
{
  double sum = 0;
  double count = 0.0;
  for (vector<double>::iterator it = (s).begin();it != (s).end();++it) {
    //if(!isnan((double)*it) && !isinf((double)*it)) {
      sum += *it;
      count += 1.0;
    //}
  }
  if (count > 0.0) sum /= count;
  return sum;
}

//getSquareMean:return the square mean for elements in s
double getSquareMean(vector<double> &s)
{
  double sum = 0;
  double count = 0.0;
  for (vector<double>::iterator it = (s).begin();it != (s).end();++it) {
    //if(!isnan(*it) && !isinf(*it)) {
      sum += (*it) * (*it);
      count += 1.0;
    //}
  }
  if (count > 0.0) sum /= count;
  return sum;
}

//getSigma:return the sigma for elements in s
double getSigma(vector<double> &s)
{
  if(s.size() <= 1)return 0;
  double sm = getSquareMean(s);
  double m = getMean(s);
  return sqrt(sm - m * m);
}
void prepareAtrocities()
{
  ifstream in("data/training_events.txt");
  char line[kMaxLine1];
  //int linenum = 0;
  while (true) {
    in.getline(line, kMaxLine1);
    int num = strlen(line);
    if (num <= 0) break;
    line[num - 1] = '\0';
    string tmp_str(line, num);
    istringstream sin(tmp_str);
    int tid;
    string pos1, pos2;
    string cid, rid;
    sin>>tid>>pos1>>pos2>>cid>>rid;

    string tmp = pos1 + string(" ") + pos2 + string(" ") + cid + string(" ") + rid;
    if (all_atrocities_mp.count(tid) == 0) {
      all_atrocities_mp[tid] = vector<string>(1, tmp);
    } else {
      all_atrocities_mp[tid].push_back(tmp);
    }
    //cout<<"tid:"<<tid<<' '<<tmp_str;cin.get();
  }
}

void readSociopolitical(int di, vector<string> &sociopolitical_data)
{
  ifstream in((string("data/training_data/" + toString(di) + string(".txt")).c_str()));
  if (in.fail()) {
    cout<<string("data/training_data/" + toString(di) + string(".txt")).c_str()<<endl;fflush(stdout);
    cout<<"something wrong"<<endl;
    exit(0);
  }
  //char line[kMaxLine1];
  char *data = new char[kMaxLine1];                    
  int linenum = 0;
  while (true) {
    in.getline(data, kMaxLine1, '\n');

    int num = strlen(data);
    if (num <= 0) break;
    assert(num < kMaxLine1);
    linenum++;
    data[num] = '\0';
    //cout<<data<<' ';cin.get();
    string tmp_str(data, num);
    sociopolitical_data.push_back(tmp_str);
    //cout<<tmp_str<<endl;cin.get();
  }
  delete[]data;
  //cout<<"total line:"<<linenum<<endl;
}
void readGeographical(vector<string> &vs)
{
  vs.clear();
  ifstream in;                                                                  
  in.open("data/regions.txt",ios_base::in|ios_base::binary);                    
  char *data = new char[kMaxLine];                    
  int linenum = 0;
  while (true) {
    in.getline(data, kMaxLine, '\n');
    int num = strlen(data);
    if (num <= 0) break;
    assert(num < kMaxLine - 1);
    linenum++;
    data[num - 1] = '\0';
    //cout<<data<<' ';cin.get();;
    string tmp_str(data, num);
    vs.push_back(tmp_str);
  }
  delete[]data;
}

void test()
{
  //int s1 = 11284;
  //int s2 = 14936;
  //int s3 = 15636;//learn[s1-s2] predict[s2-s3]
  int s1 = 11284;
#if 0
  int s2 = 12000;
  //int s3 = 13000;
  int s3 = 15636;//learn[s1-s2] predict[s2-s3]
#else
  int s2 = 15302;
  int s3 = 15636;//learn[s1-s2] predict[s2-s3]
#endif
  vector<string> data;
  prepareAtrocities();
  readGeographical(data);

  class MassAtrocityPredictor myobj;
  double score = 0.0;
  int region_num = 3671;
  vector<int> last_di(region_num, -1);
  vector<double> pnl;
  vector<double> cache;
  for (int i = s1;i <= s3;++i) {
    cout<<i<<":"<<s3<<endl;fflush(stdout);
    if (i == s1) {
      myobj.receiveData(2, i, data);
    }
    
    vector<string> sociopolitical_data;
    vector<string> atrocities_data;
    readSociopolitical(i, sociopolitical_data);
    if (all_atrocities_mp.count(i) > 0) {
      atrocities_data = all_atrocities_mp[i];
    }


    myobj.receiveData(1, i, sociopolitical_data);
    //for (int s = 0;s < atrocities_data.size();++s) {
    //  cout<<atrocities_data[s]<<' ';
    //}
    //cin.get();
    myobj.receiveData(0, i, atrocities_data);

    for (int k = 0;k <=0;++k) {
      if (all_atrocities_mp.count(k + i) > 0) {
        vector<string> &ts = all_atrocities_mp[k + i];
        for (int s = 0;s < ts.size();++s) {
          istringstream sin(ts[s]);
          int did, cid, rid;
          string pos1, pos2;
          sin>>did>>pos1>>pos2>>cid>>rid;
          last_di[rid] = i;
        }
      }
    }

    if (i >= s2) {
      vector<double> ans = myobj.predictAtrocities(i);
      vector<bool> ishappen(ans.size(), false);
      for (int k = 1;k <=30;++k) {
        if (all_atrocities_mp.count(k + i) > 0) {
          vector<string> &ts = all_atrocities_mp[k + i];
          for (int s = 0;s < ts.size();++s) {
            istringstream sin(ts[s]);
            int did;
            int cid;
            int rid;
            string pos1, pos2;
            sin>>did>>pos1>>pos2>>cid>>rid;
            ishappen[rid] = true;
          }
        }
      }
      double inc = 0;
      for (int j = 0;j < ans.size();++j) {
        int D = last_di[j];
        double WGH =  (D < 0)?1:(tanh(((double)(i - D + 10) / 180.0)));
        double CONF = ans[j];
        //cout<<"WGH:"<<WGH<<' ';
        //cout<<"CONF:"<<CONF<<' ';
        if (ishappen[j]) {
          score += WGH * (CONF - CONF * CONF / 2);
          inc += WGH * (CONF - CONF * CONF / 2);
        } else {
          score -= WGH * CONF * CONF / 2;
          inc -= WGH * CONF * CONF / 2;
        }
        //cout<<score<<' ';cin.get();
      }
      pnl.push_back(inc);
      cout<<score<<endl;
      cache.push_back(score);

      //for (int i = 0;i < ans.size();++i) {
      //  cout<<ans[i]<<' ';
      //}
      //cin.get();
    }

  }
  ofstream out0("score");
  for (size_t i = 0;i < cache.size();++i) {
    out0<<cache[i]<<endl;
  }
  ofstream out1("pnl");
  for (size_t i = 0;i < cache.size();++i) {
    out1<<pnl[i]<<endl;
  }
  cout<<"ir is:"<<getMean(pnl) / getSigma(pnl)<<endl;
}
int main()
{
  test();
  return 0;
}
#endif
}