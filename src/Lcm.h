#ifndef _LCM_H_
#define _LCM_H_

#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <queue>
#include <algorithm>
#include <cassert>
#include <cmath>
#include "Database.h"
#include "OccurenceDeriver.h"

using namespace std;

class Feature {
  public:
    vector<int> transactionList;
    vector<int> itemsets;
    vector<int> result;
    double param;
    double gain;

    Feature(){
      gain = 0;
    }
    void Print(ostream &out) {
      //    out << "PARAM:" << param << endl;
      out << "ITEMSETS:";
      for (unsigned int i = 0; i < (unsigned int)itemsets.size(); i++) {
        out << itemsets[i] << " ";
      }
      out << endl;
      out << "TRANSACTIONID:";
      for (unsigned int i = 0; i < (unsigned int)transactionList.size(); i++) {
        out << transactionList[i] << " ";
      }
      out << endl;
    }
    void clear() {
      transactionList.clear();
      itemsets.clear();
      result.clear();
      gain  = 0;
    }
};
/*
   bool operator<(const Feature& a, const Feature& b){
   return a.gain < b.gain;
   };
   */
class Prt {
  public:
    bool operator()(const Feature &x, const Feature &y) const {
      return x.gain > y.gain;
    }
};
class Cmp {
  public:
    bool operator()(const pair<double, int> &x, const pair<double, int> &y) const {
      return x.first > y.first;
    }
};

class Lcm {
  ostream &out;
  public:
  int min_sup;     // minimum support
  int maxItem;     // maximum item
  int minItem;     // minimum item
  int max_pat;     // maximum size of itemsets
  bool bitmapFlag; // bit bitmapFlag

  vector<int> totalItem;
  Lcm(ostream &out, int _min_sup, int _max_pat);
  ~Lcm();
  void RunLcm(Database &database);
  void LcmIter(Database &database, vector<int> &itemsets, vector<int> &transactionList, OccurenceDeriver &occ, vector<int> &freqList, int boundType);
  void PrintItemsets(const vector<int> &itemsets, const OccurenceDeriver &occ);
  int CalcurateCoreI(const vector<int> &itemsets, const vector<int> &freqList);
  bool PpcTest(const Database &database, vector<int> &itemsets, vector<int> &transactionList, int item, vector<int> &newTransactionList);
  void MakeClosure(const Database &database, vector<int> &transactionList, vector<int> &q_sets, vector<int> &itemsets, int item);
  bool CheckItemInclusion(const Database &database, vector<int> &transactionList, int item);
  vector<int> CalcTransactionList(const Database &database, const vector<int> &transactionList, int item);
  void CalcTransactionList(const Database &database, const vector<int> &transactionList, int item, vector<int> &newTransactionList);
  void UpdateTransactionList(const Database &database, const vector<int> &transactionList, const vector<int> &q_sets, int item, vector<int> &newTransactionList);
  void UpdateFreqList(const Database &database, const vector<int> &transactionList, const vector<int> &gsub, vector<int> &freqList, int freq, vector<int> &newFreq);
  void UpdateOccurenceDeriver(const Database &database, const vector<int> &transactionList, OccurenceDeriver &occurence);

  // parameters for LARS
  unsigned int l;
  double wbias;
  double maxGain;
  double nextGain;
  double d;
  vector<double> q;
  vector<double> y;
  double sigma;
  unsigned int topL;
  unsigned int maxiter;
  int countTop;
  Feature feature;
  vector<Feature> features;
  priority_queue<Feature, vector<Feature>, Prt > currentFeatures;
  int featureNum;
  string output_filename;
  string validation_filename;

  double eta0;
  double rho0;
  vector<double> uv;
  vector<double> uw;

  time_t mine_start_time;
  time_t mine_end_time;
  double total_mine_time;
  time_t numerical_start_time;
  time_t numerical_end_time;
  double numerical_time;
  int    pruning_number;

  // methods for LARS
  Database database, validation_database;
  void Init(Database &_database, Database &_validation_database, string &_output_filename, unsigned int _topL, unsigned int _iter);
  void Run();
  void LcmWrapper(int Type);
  bool GainTest(const vector<int> &transactionList, int Type);
  //  bool GainTestLars(const vector<int> &transactionList);
  void AddItem(const vector<int> itemsets, const vector<int> &transactionList, int Type);
};
#endif // _LCM_H_
