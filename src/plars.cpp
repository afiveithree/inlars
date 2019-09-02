#include "Lcm.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Database.h"
#include <algorithm>

using namespace Eigen;
template <typename T> inline T sgn(T a) { return (a < T(0)) ? T(-1) : T(1); }

/************************************************************************
 * Init
 ************************************************************************/
void Lcm::Init(Database &_database, Database &_validation_database, string &_output_filename, unsigned int _topL, unsigned int _maxiter) {
  database        = _database;
  validation_database        = _validation_database;
  output_filename = _output_filename;
  topL            = _topL;
  maxiter            = _maxiter;
}

/*************************************************************************
 * Run
 *************************************************************************/
void Lcm::Run() {
  clock_t allstart = clock();

  l = database.GetNumberOfTransaction(); // # of transaction
  OccurenceDeriver occ(database);
  maxItem = database.GetMaxItem();
  minItem = database.GetMinItem();

  pruning_number  = 0;
  nodes_number    = 0;
  total_mine_time = 0;
  numerical_time         = 0;

  y.resize(l);
  for (unsigned int i = 0; i < l; i++)
    y[i] = database.database[i].label;

  // try {

    // Normalize target response values
    double sumY = 0;
    for (unsigned int i = 0; i < l; i++)
      sumY += y[i];
    double meanY = sumY/l;
    //    double stdy = 0.0;
    //    for (unsigned int i = 0; i < l; i++)
    //      stdy += (y[i]-meanY)*(y[i]-meanY);
    //    stdy = sqrt(stdy/(l-1));
    //    for (unsigned int i = 0; i < l; i++)
    //      y[i] = (y[i]-meanY)/stdy;

    //    for (unsigned int i = 0; i < l; i++)
    //      y[i] = (y[i]-meanY);

    q.resize(l);    uw.resize(l); uv.resize(l);
    for (unsigned int i = 0; i < l; i++){
      q[i] = -y[i];
    }

    mine_start_time = clock();
    LcmWrapper(0); // pls bound
    mine_end_time   = clock();

    Feature &f = features[0];
    cout << "---variable ";
    for (size_t j = 0; j < f.itemsets.size(); j++){
      cout << f.itemsets[j] << " ";
    }
    cout << "added---" << endl;

    rho0 = 0; eta0 = 0;
    MatrixXd Xa(l,1);
    Map<VectorXi > xrow(&(features[0].result[0]),l);
    Xa.col(0) = xrow.cast<double>();

    lambda = xrow.transpose().cast<double>()*y;
    VectorXd x(1); x.setZero();
    VectorXd z(1);
    uh.resize(l);
    ug.resize(l);

    unsigned int removed = 0;

#ifdef DEBUG
    for(unsigned int i = 0; i < features.size(); i++){
      Feature &f = features[i];
      cout << "No:" << i << " gain: " << features[i].gain << " Itemsets: ";
      for (size_t j = 0; j < f.itemsets.size(); j++)
        cout << f.itemsets[j] << " ";
      cout << endl;
      cout << "occurrance:";
      for (int j = 0; j < (int)features[i].result.size(); j++)
        cout << features[i].result[j];
      cout << endl << endl;;
    }
#endif

    VectorXd pred(l);
    pred = Xa * x;
    double var = 0.0;
    double rss = 0.0;
    for(unsigned int i = 0; i < l; i++){
      var += (y[i] - meanY)*(y[i] - meanY);
      rss += (y[i] - pred[i])*(y[i] - pred[i]);
    }
    double Q2 = 1.0 - rss/var;
    cout << "ITERATION: 0 Train RSS: " << rss << " Q2: " << Q2 << " var " << var << endl << endl;
    double Q2_prev = Q2;
    double Q2_val_prev = 0;

    //MAIN ITERATION
    for(unsigned int itr = 0; itr < maxiter && lambda > 0; itr++){
      if(itr >= l*3){
        cerr << "Iteration limit reached" << endl;
        End_clock(allstart);
        exit(0);
      }
      unsigned int n = x.size();

      z.resize(n);
      z = (Xa.transpose() * Xa).ldlt().solve(Xa.transpose() * y); // use llt().solve for large matrix
      ug = Xa * x - y;
      uh = Xa * z - y;

#ifdef DEBUG
      cout <<"Xat x Xa" << endl << Xa.transpose() * Xa << endl;
      cout <<"(Xat x Xa) -1" << endl << (Xa.transpose() * Xa).inverse() << endl;
      cout << "y:" << y.transpose() << endl;
      cout << "ug:" << ug.transpose() << endl;
      cout << "z: " << z.transpose() << endl;
      cout << "uh: " << uh.transpose() << endl;

      vector<double> adv; adv.resize(n);
      double madv = 0;
      for(unsigned int i = 0; i < n; i++){
        adv[i] = 0;
        for(unsigned int j = 0; j < l; j++){
          adv[i] += Xa(i,j)*ug[j];
        }
        adv[i] = fabs(adv[i]);
        madv += adv[i];
      }
      madv /= n;
      double ss = 0;
      for(unsigned int i = 0; i < n; i++){
        ss += (adv[i]-madv)*(adv[i]-madv);
      }
      ss /= n;
      if(sqrt(ss) > 0.00000000001){
        cout << "SS " << ss << endl;
        //  exit(0);
      }
#endif

      mine_start_time = clock();
      LcmWrapper(1); // lars bound
      mine_end_time   = clock();
      total_mine_time += (mine_end_time - mine_start_time);
      numerical_start_time = clock();

      // Search in active set
      double mind2 = 1000000000;
      int minid2 = -1;
      for(unsigned int j = 0; j < n; j++){
        double d2 = 0;
        if(z[j] < 0){
          d2 = x[j]/(x[j]-z[j]);
        }
#ifdef DEBUG
        cout << "x[j] " << x[j] << " z[j] " << z[j] << " d2=x[j]/(x[j]-z[j]) " << d2 << endl;
#endif
        if( (d2 > 0.000001) && (d2 < mind2) ){
          mind2 = d2;
          minid2 = j;
        }
      }
      double mind1 = maxGain;
#ifdef DEBUG
      cout << "mind1:" << mind1 << " mind2:" << mind2 << endl;
#endif

      if(mind1 > 1e4 && mind2 > 1e4){
        cout << "NO MORE FEATURES" << endl;
        End_clock(allstart);
        exit(0);
      }

      if(mind1 < mind2){
        x += mind1*(z - x); // int n = x.size();
        x.conservativeResize(n+1); x(n) = 0;

        new (&xrow) Map<VectorXi > (&(features[itr+1-removed].result[0]), l);
        Xa.conservativeResize(NoChange, n+1); Xa.col(Xa.cols()-1).setZero();
        Xa.col(Xa.cols()-1) = xrow.cast<double>();
#ifdef DEBUG
        cout << "added" << endl;
        cout << "Xa" << endl << Xa << endl;
#endif
        lambda *= (1-mind1);

        Feature &f = features[itr+1-removed];
        cout << endl << "---variable ";
        for (size_t j = 0; j < f.itemsets.size(); j++){
          cout << f.itemsets[j] << " ";
        }
        cout << "added---" << endl;


      }else{ // mind1 > mind2
        Xa.rightCols(n-minid2).leftCols(n-1-minid2) = Xa.rightCols(n-1-minid2);
        Xa.conservativeResize(NoChange, n-1); // better with permutation?

        x += mind2*(z - x); // int n = x.size();
        x.segment(minid2, n-1-minid2) = x.tail(n-1-minid2);
        x.conservativeResize(n-1);

        z.segment(minid2, n-1-minid2) = z.tail(n-1-minid2);
        z.conservativeResize(n-1);

        // h.segment(i, n-1-minid2) = h.tail(n-1-minid2);
        // h.conservativeResize(n-1);
#ifdef DEBUG
        cout << "removed" << endl;
        cout << "Xa" << endl << Xa << endl;
        cout << "x" << endl << x << endl;
        cout << "z" << endl << z << endl;
#endif
        lambda *= (1-mind2);

        vector<Feature>::iterator _itr5 = features.begin() + minid2;
        Feature &f = features[minid2];
        cout << endl << "---variable ";
        for (size_t k = 0; k < f.itemsets.size(); k++)
          cout << f.itemsets[k] << " ";
        cout << "removed---" << endl;

        features.erase(_itr5);
        ++removed;
        //  cout << "No." << minid2 << ": *** REMOVED ***" << endl;
      }

#ifdef DEBUG
      for(unsigned int i = 0; i < features.size(); i++){
        Feature &f = features[i];
        cout << "No:" << i << " gain: " << features[i].gain << " Itemsets: ";
        for (size_t j = 0; j < f.itemsets.size(); j++)
          cout << f.itemsets[j] << " ";
        cout << endl;
        for (int j = 0; j < (int)features[i].result.size(); j++)
          cout << features[i].result[j];
        cout << endl << endl;;
      }
#endif

      n = x.size();

#ifdef DEBUG
      cout << "ug" << ug.transpose() << endl;
      cout << "xrow " << xrow.transpose() << endl;
      cout << "ITERATION: " << itr+1 << " mind1: " << mind1 << " mind2: " << mind2 << endl;
      cout << "lambda:" << lambda << endl;
#endif

      numerical_end_time = clock();
      numerical_time += numerical_end_time - numerical_start_time;

      VectorXd pred(l);
      pred = Xa * x;
      double var = 0.0;
      double rss = 0.0;
      for(unsigned int i = 0; i < l; i++){
        var += (y[i] - meanY)*(y[i] - meanY);
        rss += (y[i] - pred[i])*(y[i] - pred[i]);
#ifdef DEBUG
        cout << "y " << y[i] << " p " << pred[i] << endl;
#endif
      }
      double Q2 = 1.0 - rss/var;

      cout << "ITERATION: " << itr+1;
      cout << " Train RSS: " << rss << " Q2: " << Q2 << endl;

      stringstream s_itr;
      s_itr << itr+1;
      string outfilename = output_filename + s_itr.str();
      std::ofstream ofs((outfilename).c_str());
      for (size_t i = 0; i < features.size(); i++) {
        Feature &f = features[i];
        ofs << "NO:" << i << std::endl;
        ofs << "WEIGHT:";
        //  for(unsigned int i = 0; i < beta.size(); i++){
        ofs << x[i] << " ";
        //  }
        ofs << std::endl;
        ofs << "ITEMSETS:";
        for (size_t j = 0; j < f.itemsets.size(); j++)
          ofs << f.itemsets[j] << " ";
        ofs << std::endl;
        ofs << "TRANSACTIONLIST:";
        for (size_t j = 0; j < f.transactionList.size(); j++)
          ofs << f.transactionList[j] << " ";
        ofs << std::endl;
      }
      ofs.close();
      // NOTE: train-side diff check
      if(abs((Q2-Q2_prev)/Q2) < 0.001){
        End_clock(allstart);
        exit(0);
      }
      Q2_prev=Q2;

      unsigned int m = validation_database.GetNumberOfTransaction();
      if(m > 0){
        vector<double> val_y(m);
        vector<double> val_p(m);
        for(unsigned int i = 0; i < m; i++){
          val_y[i] = validation_database.database[i].label;

          Transaction &transaction = validation_database.database[i];
          vector<int> &itemset = transaction.itemsets;
          sort(itemset.begin(), itemset.end());

          double p = 0.0;
          for(size_t j = 0; j < features.size(); j++){
            Feature &f = features[j];
            size_t k;
            for (k = 0; k < f.itemsets.size(); k++){
              if( !binary_search(itemset.begin(), itemset.end(), f.itemsets[k])) break;
            }
            if( k == f.itemsets.size() ){
              p += x[j];
            }
          }// end for j
          val_p[i] = p; // l1nls fix : no normalize
//          val_p[i] = p * stdy; // l1nls fix : var normalize only
#ifdef DEBUG
          cout << "y " << val_y[i];
          cout << " p " << val_p[i] << endl;
#endif
        } // end for i

        double mean = 0;
        for(unsigned int i = 0; i < m; i++) mean += val_y[i];
        mean /= m;
        for(unsigned int i = 0; i < m; i++) val_y[i] -= mean;
        double var = 0.0;
        double rss = 0.0;
        for(unsigned int i = 0; i < m; i++){
          var += (val_y[i] - mean)*(val_y[i] - mean);
          rss += (val_y[i] - val_p[i])*(val_y[i] - val_p[i]);
        }
        double Q2 = 1.0 - rss/var;
        cout << "ITERATION: " << itr+1;
        cout << " Validation RSS: " << rss << " Q2: " << Q2 << endl;

      } //end if m > 0
// // NOTE: validation-side diff check
//      if(abs((Q2-Q2_val_prev)/Q2) < 0.001){
//        End_clock(allstart);
//        exit(0);
//      }
//      Q2_val_prev=Q2;
    } // end MAIN ITERATION

    End_clock(allstart);


  // } // end try

  // GMM_STANDARD_CATCH_ERROR;
}

/***************************************************************************
 * LcmWrapper
 ***************************************************************************/
void Lcm::LcmWrapper(int boundType) {
  OccurenceDeriver occ(database);
  maxItem = database.GetMaxItem();
  minItem = database.GetMinItem();

  int prev_nodes_number   = nodes_number;
  int prev_pruning_number = pruning_number;

  vector<int> itemsets;
  vector<int> transactionList;
  for (unsigned int i = 0; i < (unsigned int)database.GetNumberOfTransaction(); i++) {
    transactionList.push_back(i);
  }

  totalItem = database.GetItemset();
  vector<int> freqList;
  countTop = 0;
  if(boundType == 0){ // first time pls
    //MAXIMIZE maxGain
    maxGain  = -1000000000;
    nextGain = -1000000000;

    LcmIter(database, itemsets, transactionList, occ, freqList, boundType);
    //  sort(currentFeatures.begin(), currentFeatures.end(), Prt());

    while (currentFeatures.size() > 1) currentFeatures.pop();
    features.push_back(currentFeatures.top()); // choose the last (largest) one

  }else{ // lars
    //*MINIMIZE* maxGain
    maxGain = 1000000000;
    nextGain = 1000000000;
    wbias = 0.0;

    LcmIter(database, itemsets, transactionList, occ, freqList, boundType);
    //  sort(currentFeatures.begin(), currentFeatures.end(), Prt());

    features.push_back(currentFeatures.top()); // choose the first (smallest) one
    while (!currentFeatures.empty()) currentFeatures.pop();

  }
  nodes_number   = max(nodes_number,prev_nodes_number);
  pruning_number = max(pruning_number,prev_pruning_number);
}

/********************************************************************************
 * AddItem
 ********************************************************************************/
void Lcm::AddItem(const vector<int> itemsets, const vector<int> &transactionList, int Type) {
  vector<int> result(l, 0);
  for (unsigned int i = 0; i < (unsigned int)transactionList.size(); i++){
    result[transactionList[i]] = +1;
  }

  double gain;
  if(Type == 0){ // pls
    gain = 0.0;
    for (unsigned int i = 0; i < l; i++)
      gain += y[i] * result[i];
    // gain = fabs(gain); // l1nls search only positive area
#ifdef DEBUG
    cout << "gain " << gain << endl;
#endif
    if (itemsets.size() == 0)
      return;
    if (gain > maxGain){
      maxGain = gain;

      feature.clear();
      feature.transactionList = transactionList;
      feature.result          = result;
      feature.param           = 1;
      feature.itemsets        = itemsets;
      feature.gain            = gain;

      currentFeatures.push(feature);
    }

    if ((unsigned int)currentFeatures.size() == (unsigned int)topL && gain < nextGain)
      return;


  }else{ // lars
    double g = 0; double h = 0;
    for (unsigned int i = 0; i < l; i++){
      g += ug[i] * result[i];
      h += uh[i] * result[i];
    }
    if (h >= 0)
      return;
    double d = (lambda + g) / (lambda + g - h);

    if (itemsets.size() == 0)
      return;

    //check redundancy
    for (size_t i = 0; i < features.size(); i++) {
      Feature &f = features[i];
      if(transactionList.size() == f.transactionList.size()){
        double sum = 0;
        for(unsigned int j = 0; j < l; j++)
          sum += f.result[j]*result[j];

        if(abs(sum)==transactionList.size()){
#ifdef DEBUG
          cout << "Itemset " << i << ": ";
          for (size_t j = 0; j < f.itemsets.size(); j++)
            cout << f.itemsets[j] << " ";
          cout << "is SAME:" << sum << endl;
#endif
          return;
        }
      }
    }
#ifdef DEBUG
    cout << "additem d " << d << " maxGain " << maxGain << endl;
#endif
    if(d > 0.001){
      gain = d;
      if(d < maxGain){
        maxGain = d;
        feature.clear();
        feature.transactionList = transactionList;
        feature.result          = result;
        feature.param           = 1;
        feature.itemsets        = itemsets;
        feature.gain            = gain;
        currentFeatures.push(feature);
      }
    }

    if ((unsigned int)currentFeatures.size() == (unsigned int)topL && d > nextGain)
      return;

  } // end lars

  nextGain = maxGain;
#ifdef DEBUG
  cout << "new maxGain " << maxGain << endl;
#endif
  //  cout << gain << " nextGain:" << nextGain << " maxGain:" << maxGain << endl;
  //  if ((unsigned int)currentFeatures.size() > (unsigned int)topL) {
  //      cout << "min:" << currentFeatures.top().gain << endl;
  //    currentFeatures.pop();
  //  }
  //  features.push_back(feature);
  //  nextGain = currentFeatures.top().gain;
}

inline void Lcm::End_clock(clock_t allstart){
  clock_t allend = clock();
  cout << endl << "***END***" << endl;
  cout << "ALL TIME: " << (double)(allend -   allstart)/(double)CLOCKS_PER_SEC << endl;
  cout << "NUMERICAL TIME: " << (double)numerical_time/(double)CLOCKS_PER_SEC << endl;
  cout << "MINING TIME: " << (double)total_mine_time/(double)CLOCKS_PER_SEC << endl;
  cout << "Inf::" << "maxpat: " << max_pat << " minsup: " << min_sup << " pruning number:" << pruning_number << " nodes number:" << nodes_number<< endl;
}
