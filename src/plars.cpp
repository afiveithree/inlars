#include "Lcm.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "gmm.h"
#include "gmm_dense_qr.h"
#include "Database.h"
#include <algorithm>

using namespace gmm;

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
  total_mine_time = 0;
  numerical_time         = 0;

  y.resize(l);
  for (unsigned int i = 0; i < l; i++)
    y[i] = database.database[i].label;

  try {

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

    for (unsigned int i = 0; i < l; i++)
      y[i] = (y[i]-meanY);

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
    vector<vector<unsigned int> > Xa;
    vector<unsigned int> xrow;
    xrow.resize(l); fill(xrow.begin(),xrow.end(),0);

    double sumGAinv = 0;
    const vector<int> &result = features[0].result;
    for (unsigned int i = 0; i < l; i++) {
      xrow[i]=result[i]; // binary occurrance vector
      sumGAinv += xrow[i];
    }
    Xa.push_back(xrow);

    double corr = 0.0;
    for(unsigned int i = 0; i < l; i++){
      corr -= xrow[i]*y[i];
    }
    vector<double> gamma; gamma.push_back(-sgn(corr));
    vector<double> s; s.push_back(-sgn(corr));
    vector<double> beta; beta.push_back(0.0);
    unsigned int removed = 0;

#ifdef DEBUG
    for(unsigned int i = 0; i < features.size(); i++){
      Feature &f = features[i];
      cout << "No:" << i << " gain: " << features[i].gain << " Itemsets: ";
      for (size_t j = 0; j < f.itemsets.size(); j++)
        cout << f.itemsets[j] << " ";
      cout << endl;
      cout << "occurrance:" << features[i].result << endl << endl;;
    }

    vector<double> pred; pred.resize(l);
    fill(pred.begin(), pred.end(), 0);
    for(unsigned int i = 0; i < l; i++)
      for(unsigned int j = 0; j < 1; j++)
        pred[i] += Xa[j][i]*beta[j];
    double var = 0.0;
    double rss = 0.0;
    for(unsigned int i = 0; i < l; i++){
      var += (y[i] - meanY)*(y[i] - meanY);
      rss += (y[i] - pred[i])*(y[i] - pred[i]);
    }
    double Q2 = 1.0 - rss/var;
    cout << "Train Q2 " << Q2 << " rss " << rss << " var " << var << endl;
#endif

    gmm::dense_matrix<double> C;

    //MAIN ITERATION
    for(unsigned int itr = 0; itr < maxiter; itr++){
      if(itr >= l*3){
        cerr << "Iteration limit reached" << endl;
        exit(0);
      }
      unsigned int n = beta.size();
      double maxabsuw = 0;
      int maxi = -1;

      for(unsigned int j = 0; j < n; j++)
        gamma[j] /= sqrt(sumGAinv);

      for(unsigned int i = 0; i < l; i++){
        uw[i] = -y[i];
        for(unsigned int j = 0; j < n; j++){
          uw[i] += beta[j]*Xa[j][i];
        }
        if(abs(uw[i]) > maxabsuw){
          maxi = i;
          maxabsuw = abs(uw[i]);
        }
        uv[i] = 0;
        for(unsigned int j = 0; j < n; j++){
          uv[i] -= gamma[j]*Xa[j][i];
        }
      }

      int min_size = 10000000;
      int min_id = -1;
      for(unsigned int j = 0; j < n; j++){
        int cur_size = 0;
        for(unsigned int i = 0; i < l; i++){
          cur_size += Xa[j][i];
        }
        if(cur_size < min_size) min_size = cur_size;
        min_id = j;
      }

      rho0 = 0; eta0 = 0;
      for(unsigned int i = 0; i < l; i++){
        rho0 += uw[i]*Xa[min_id][i];
        eta0 += uv[i]*Xa[min_id][i];
      }


#ifdef DEBUG
      cout << "y:" << y << endl;
      cout << "uw:" << uw << endl;
      cout << "MAXCORR : " << maxabsuw << " at " << maxi << endl;
      cout << "uv: " << uv << endl;
      cout << "rho0:" << rho0 << "eta0:" << eta0 << endl;

      vector<double> adv; adv.resize(n);
      double madv = 0;
      for(unsigned int i = 0; i < n; i++){
        adv[i] = 0;
        for(unsigned int j = 0; j < l; j++){
          adv[i] += Xa[i][j]*uw[j];
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
        if(gamma[j] != 0){
          d2 = -beta[j]/gamma[j];
        }
#ifdef DEBUG
        cout << "beta[j] " << beta[j] << " gamma[j] " << gamma[j] << " d2=-beta[j]/gamma[j] " << d2 << endl;
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
        exit(0);
      }

      if(mind1 < mind2){
        for(unsigned int j = 0; j < n; j++){
          beta[j] += mind1*gamma[j];
        }
        beta.push_back(0.0);

        fill(xrow.begin(), xrow.end(), 0);
        const vector <int> &result = features[itr+1-removed].result;
        for (unsigned int i = 0; i < l; i++) {
          xrow[i] = result[i]; // binary occurance vector
        }
        Xa.push_back(xrow);

        corr = 0;
        for (unsigned int i = 0; i < l; i++) {
          corr += uw[i]*xrow[i];
        }
        s.push_back(-sgn(corr));


        Feature &f = features[itr+1-removed];
        cout << "---variable ";
        for (size_t j = 0; j < f.itemsets.size(); j++){
          cout << f.itemsets[j] << " ";
        }
        cout << "added---" << endl;


      }else{ // mind1 > mind2
        for(unsigned int j = 0; j < n; j++){
          beta[j] += mind2*gamma[j];
        }

        vector<vector<unsigned int> >::iterator _itr1 = Xa.begin() + minid2;
        Xa.erase(_itr1);
        vector<double>::iterator _itr2 = beta.begin() + minid2;
        beta.erase(_itr2);
        vector<double>::iterator _itr3 = s.begin() + minid2;
        s.erase(_itr3);

        vector<Feature>::iterator _itr5 = features.begin() + minid2;
        Feature &f = features[minid2];
        cout << "---variable ";
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
        cout << "occurrance:" << features[i].result << endl;
      }
#endif

      n = beta.size();
      C.resize(n,n);
      for(unsigned int i = 0; i < n; i++){
        for(unsigned int j = 0; j < n; j++){
          C(i,j) = 0;
          for(unsigned int k = 0; k < l; k++){
            C(i,j) += Xa[i][k]*Xa[j][k];
          }
        }
      }

      gamma.resize(n);
      //      gmm::lu_solve(C,gamma,s); // Use this if contition number is larg
      gmm::lu_inverse(C);
      gmm::mult(C,s,gamma);

#ifdef DEBUG
      cout << "C " << C << endl;
      cout << "condition number " << gmm::condition_number(C) << endl;
      cout << "uw" << uw << endl;
      cout << "x " << xrow << endl;
      cout << "s " << s << endl;
      cout << "ITERATION: " << itr+1 << " mind1: " << mind1 << " mind2: " << mind2 << endl;
      cout << "gamma:" << gamma << endl;
#endif

      for(unsigned int i = 0; i < gamma.size(); i++){
        if(isnan(gamma[i]) || isinf(gamma[i])){
          cout << "STOP" << endl;
          exit(0);
        }
      }

      sumGAinv = 0;
      if(n > 1){
        //  gmm::lu_inverse(C);
        for(unsigned int a = 0; a < n; a++)
          for(unsigned int b = 0; b < n; b++)
            sumGAinv += s[a]*s[b]*C(a,b); // Note that C is inversed
      }

      numerical_end_time = clock();
      numerical_time += numerical_end_time - numerical_start_time;

      vector<double> pred; pred.resize(l);
      fill(pred.begin(), pred.end(), 0);
      for(unsigned int i = 0; i < l; i++)
        for(unsigned int j = 0; j < n; j++)
          pred[i] += Xa[j][i]*beta[j];
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

      cout << endl << "ITERATION: " << itr+1;
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
        ofs << beta[i] << " ";
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
              p += beta[j];
            }
          }// end for j
          val_p[i] = p;
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
    } // end MAIN ITERATION

    clock_t allend = clock();
    cout << endl << "***END***" << endl;
    cout << "ALL TIME: " << (double)(allend - allstart)/(double)CLOCKS_PER_SEC << endl;
    cout << "NUMERICAL TIME: " << (double)numerical_time/(double)CLOCKS_PER_SEC << endl;
    cout << "MINING TIME: " << (double)total_mine_time/(double)CLOCKS_PER_SEC << endl;
    //      cout << "Inf::" << "maxpat: " << max_pat << " minsup: " << min_sup << " pruning number:" << pruning_number << endl;


  } // end try

  GMM_STANDARD_CATCH_ERROR;
}

/***************************************************************************
 * LcmWrapper
 ***************************************************************************/
void Lcm::LcmWrapper(int boundType) {
  OccurenceDeriver occ(database);
  maxItem = database.GetMaxItem();
  minItem = database.GetMinItem();

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
      gain += q[i] * result[i];
    gain = fabs(gain);
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
    double rho = 0; double eta = 0;
    for (unsigned int i = 0; i < l; i++){
      rho += uw[i] * result[i];
      eta += uv[i] * result[i];
    }
    double res1 = 0; double res2 = 0;// double d = 0;
    if(eta0 - eta != 0.0){
      res1 = max( (rho0 - rho) / (eta0 - eta), 0.0);
    }
    if(eta0 + eta != 0.0){
      res2 = max( (rho0 + rho) / (eta0 + eta), 0.0);
    }
    double d = min (res1, res2);

    if (itemsets.size() == 0)
      return;

    //check redundancy
    vector<double> result(l);
    for (unsigned int i = 0; i < (unsigned int)transactionList.size(); i++){
      result[transactionList[i]] = +1;
    }
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
          return false; // pruning
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
