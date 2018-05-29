#include "Lcm.h"

/************************************************************************
 * Constructor
 ************************************************************************/
Lcm::Lcm(ostream &_out, int _min_sup, int _max_pat)
  : out(_out), min_sup(_min_sup), max_pat(_max_pat) {}

/************************************************************************
 * Destructor
 ************************************************************************/
Lcm::~Lcm() {}

/*************************************************************************
 * run lcm
 *************************************************************************/
void Lcm::RunLcm(Database &database) {
  //database.RemoveItemsbyFreq(min_sup);
  //database.PrintBitmap(cout);
  OccurenceDeriver occ(database);
  maxItem = database.GetMaxItem();
  minItem = database.GetMinItem();

  vector<int> itemsets;
  vector<int> transactionList;
  for (int i = 0; i < database.GetNumberOfTransaction(); i++) {
    transactionList.push_back(i);
  }

  totalItem = database.GetItemset();
  vector<int> freqList;
  int boundType = 0;
  LcmIter(database, itemsets, transactionList, occ, freqList, boundType);
}

/*************************************************************************
 * main function
 *************************************************************************/
void Lcm::LcmIter(Database &database, vector<int> &itemsets, vector<int> &transactionList, OccurenceDeriver &occ, vector<int> &freqList, int boundType) {
  //  if ((int)currentFeatures.size() == topL)
  //      return;
  bool gainTestFlag;
  //  cout << endl;
#ifdef DEBUG
  PrintItemsets(itemsets, occ); // print frequent itemset
#endif
  gainTestFlag = GainTest(transactionList, boundType);

  if (gainTestFlag == false) {
    pruning_number += 1;
    return;
  }
  nodes_number += 1;
#ifdef DEBUG
  cout << "GainTest done" << endl;
#endif
  AddItem(itemsets, transactionList, boundType);

  //  PrintItemsets(itemsets, occ); // print frequent itemset
  int core_i = CalcurateCoreI(itemsets, freqList);
  // database reduction
  // database.Reduction(transactionList, core_i);
  // database.Print(out);

  /*
     Compute the frequency of each pattern P \cup {i}, i > core_i(P)
     by Occurence deliver with P and Occ;
   */
  vector<int>::iterator iter = lower_bound(totalItem.begin(), totalItem.end(), core_i);
  vector<int> freq_i;
  for (int i = *iter; iter != totalItem.end(); iter++, i = *iter) {
    if ((int)occ.table[i].size() >= min_sup && binary_search(itemsets.begin(), itemsets.end(), i) == false)
      freq_i.push_back(i);

  }

  vector<int> newTransactionList;
  vector<int> q_sets;
  vector<int> newFreqList;
  for (vector<int>::iterator freq = freq_i.begin(); freq != freq_i.end(); freq++) {
    newTransactionList.clear();
    if (PpcTest(database, itemsets, transactionList, *freq, newTransactionList)) {
      q_sets.clear();
      MakeClosure(database, newTransactionList, q_sets, itemsets, *freq);
      if (max_pat == 0 || (int)q_sets.size() <= max_pat) {
        newTransactionList.clear();
        UpdateTransactionList(database, transactionList, q_sets, *freq, newTransactionList);
        UpdateOccurenceDeriver(database, newTransactionList, occ);
        newFreqList.clear();
        UpdateFreqList(database, transactionList, q_sets, freqList, *freq, newFreqList);
        LcmIter(database, q_sets, newTransactionList, occ, newFreqList, boundType);
      }
    }
  }
}

/************************************************************************
 * GainTest
 * tree pruning function
 ************************************************************************/
inline bool Lcm::GainTest(const vector<int> &transactionList, int boundType)
{
  if(boundType == 0){ // pls

    double gain_plus  = 0;
    double gain_minus = 0;
    wbias = 0.0;
    //  cout << "size" << (int)transactionList.size() << endl;
    for (int i = 0; i < (int)transactionList.size(); i++) {
      wbias += q[transactionList[i]]; //instead :y
      if (q[transactionList[i]] >= 0)
        gain_plus  += (double)abs(q[transactionList[i]]);
      else
        gain_minus += (double)abs(q[transactionList[i]]);
    }

#ifdef DEBUG
    cout << "gain:" << gain_plus << " " <<  gain_minus << " " << maxGain << endl;
    cout << "wbias:" << wbias << endl;
#endif // DEBUG

    return (gain_minus < maxGain) ? false: true; // pruning

  }else{ // lars
    double g = 0; double g_plus  = 0;  double g_minus = 0;
    double h = 0; double h_plus  = 0;  double h_minus = 0;

    for (int i = 0; i < (int)transactionList.size(); i++) {
      g += ug[transactionList[i]];
      h += uh[transactionList[i]];

      if (uh[transactionList[i]] >= 0){
        h_plus  += (double)abs(uh[transactionList[i]]);
      }else{
        h_minus += (double)abs(uh[transactionList[i]]);
      }

      if(ug[transactionList[i]] >= 0){
        g_plus  += (double)abs(ug[transactionList[i]]);
      }else{
        g_minus += (double)abs(ug[transactionList[i]]);
      }
    }

    double d = (lambda + g) / (lambda + g - h);
    d = max(d,0.000001);

    // double t = h / (lambda + g);
    // double diff = - h_minus - h_plus - t*(lambda + min(-g_minus,g_plus));
    double expect = (lambda > g_minus) ? (lambda - g_minus)/(lambda - g_minus + h_minus) : 0.1;

#ifdef DEBUG
    for (int i = 0; i < (int)transactionList.size(); i++)
      cout << transactionList[i] << " ";
    cout << endl;
    cout << "g " << g << endl;
    cout << "h " << h << endl;
    cout << "nextGain " << nextGain << endl;
    cout << "d " << d << " t " << t << " lambda " << lambda << " diff " << diff << endl;
#endif

    if (g_plus < -lambda || h_minus == 0){
      return false; // pruning
    // }else if(diff > 0){
    //   return false
    }else if (expect > maxGain){
      return false;
    }else{
      return true;
    }

  }
}


/*************************************************************************
 * Update Freq List
 *************************************************************************/
inline void Lcm::UpdateFreqList(const Database &database, const vector<int> &transactionList, const vector<int> &gsub, vector<int> &freqList, int freq, vector<int> &newFreq) {
  int iter = 0;
  if (freqList.size() > 0) {
    for (; iter < (int)gsub.size(); iter++) {
      if (gsub[iter] >= freq) break;
      newFreq.push_back(freqList[iter]);
    }
  }

  vector<int> newList;
  for (int i = 0; i < (int)transactionList.size(); i++)
    newList.push_back(transactionList[i]);

  vector<int> newnewList;
  for (int i = iter; i < (int)gsub.size(); i++) {
    int item = gsub[i];
    int freqCount = 0;

    for (int list = 0; list < (int)newList.size(); list++) {
      const Transaction &transaction = database.database[newList[list]];
      if (binary_search(transaction.itemsets.begin(), transaction.itemsets.end(), item) == true) {
        freqCount += 1;
        newnewList.push_back(newList[list]);
      }
    }
    newFreq.push_back(freqCount);
    newList = newnewList;
    newnewList.clear();
  }
}

/*************************************************************************
 * Update Transaction List
 *************************************************************************/
inline void Lcm::UpdateTransactionList(const Database &database, const vector<int> &transactionList, const vector<int> &q_sets, int item, vector<int> &newTransactionList) {
  for (int i = 0; i < (int)transactionList.size(); i++) {
    const Transaction &transaction = database.database[transactionList[i]];
    int iter;
    for (iter = 0; iter < (int)q_sets.size(); iter++) {
      int q = q_sets[iter];
      if (q >= item) {
        if (binary_search(transaction.itemsets.begin(), transaction.itemsets.end(), q) == false)
          break;
      }
    }
    if (iter == (int)q_sets.size())
      newTransactionList.push_back(transactionList[i]);
  }
}

/***************************************************************************
 * Print itemsets
 ***************************************************************************/
inline void Lcm::PrintItemsets(const vector<int> &itemsets, const OccurenceDeriver &occ) {
  if ((int)itemsets.size() > 0) {
    for (int i = 0; i < (int)itemsets.size(); i++)
      out << itemsets[i] << " ";

    out << endl;

    out << "(";
    const vector<int> &table = occ.GetTable(itemsets[itemsets.size() - 1]);
    for (int i = 0; i < (int)table.size() - 1; i++)
      out << table[i] << " ";

    //     out << table[table.size() - 1] << ") " << occ.GetNumOcc(itemsets[itemsets.size()-1]) << endl;
    out << table[table.size() - 1] << ") " << endl;

  }
}

/*****************************************************************************
 * calculrate core_i
 *****************************************************************************/
inline int Lcm::CalcurateCoreI(const vector<int> &itemsets, const vector<int> &freqList) {
  if (itemsets.size() > 0) {
    int current = freqList[freqList.size() - 1];
    for (int i = (int)freqList.size() - 2; i >= 0; i--) {
      if (current != freqList[i]) return itemsets[i+1];
    }
    return itemsets[0];
  }
  else
    return 0;
}


/**************************************************************************
 * Prefix Preseaving Test
 * Test whether p(i-1) is equal to q(i-1) or not
 **************************************************************************/
inline bool Lcm::PpcTest(const Database &database, vector<int> &itemsets, vector<int> &transactionList, int item, vector<int> &newTransactionList) {
  // j_sets: set not including items which are included in itemsets.
  // make transactionList pointing to the indexes of database including P \cup item
  CalcTransactionList(database, transactionList, item, newTransactionList);

  // check j s.t j < i, j \notin P(i-1) is included in every transaction of T(P \cup {i})
  for (vector<int>::iterator j = totalItem.begin(); *j < item; j++) {
    if (binary_search(itemsets.begin(), itemsets.end(), *j) == false &&
        CheckItemInclusion(database, newTransactionList, *j) == true) {
      return false;
    }
  }
  return true;
}

/****************************************************************************
 * Make closure
 * make Q = Clo(P \cup {i}) subject to Q(i-1) = P(i-1)
 *****************************************************************************/
inline void Lcm::MakeClosure(const Database &database, vector<int> &transactionList, vector<int> &q_sets, vector<int> &itemsets, int item) {
  // make transactionList pointting to the indexes of database inclding {P \cup item}
  // vector<int> &newTransactionList = CalcTransactionList(database, transactionList, item);
  // make Clo(P \cup {i})
  for (int i = 0; i < (int)itemsets.size() && itemsets[i] < item; i++) {
    q_sets.push_back(itemsets[i]);
  }
  q_sets.push_back(item);

  vector<int>::iterator i = lower_bound(totalItem.begin(), totalItem.end(), item + 1);

  for (int iter = *i; i != totalItem.end(); i++, iter = *i) {
    if (CheckItemInclusion(database, transactionList, iter) == true) {
      q_sets.push_back(iter);
    }
  }
}

/********************************************************************************
 * CheckItemInclusion
 * Check whther item is included in the transactions pointed to transactionList
 ********************************************************************************/
inline bool Lcm::CheckItemInclusion(const Database &database, vector<int> &transactionList, int item) {
  for (vector<int>::iterator iter = transactionList.begin(); iter != transactionList.end(); iter++) {
    const Transaction &transaction = database.database[*iter];
    if (binary_search(transaction.itemsets.begin(), transaction.itemsets.end(), item) == false) return false;
  }
  return true;
}


/*********************************************************************************
 *  Calcurate new transaction list
 *********************************************************************************/
inline void Lcm::CalcTransactionList(const Database &database, const vector<int> &transactionList, int item, vector<int> &newTransactionList) {

  for (int list = 0; list < (int)transactionList.size(); list++) {
    const Transaction &transaction = database.database[transactionList[list]];
    if (binary_search(transaction.itemsets.begin(), transaction.itemsets.end(), item) == true)
      newTransactionList.push_back(transactionList[list]);
  }
}

/***********************************************************************************
 * Update Occurence Deriver
 ***********************************************************************************/
inline void Lcm::UpdateOccurenceDeriver(const Database &database, const vector<int> &transactionList, OccurenceDeriver &occurence) {
  occurence.Clear();
  for (int i = 0; i < (int)transactionList.size(); i++) {
    const Transaction &transaction = database.database[transactionList[i]];
    const vector<int> &itemsets = transaction.itemsets;
    for (int j = 0; j < (int)itemsets.size(); j++) {
      occurence.table[itemsets[j]].push_back(transactionList[i]);
    }
  }
}
