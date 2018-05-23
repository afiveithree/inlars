#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <strstream>
#include <cmath>

using namespace::std;

struct Feature {
  double weight;
  vector<int> itemset;
  void clear() {
    weight = 0.0;
    itemset.clear();
  }
};

void Usage(){
  cerr << endl
    << "Usage: predict MODEL_FILE TEST_FILE" << endl
    << endl;
  exit(0);
}

int main(int argc, char **argv) {
  ifstream ifs(argv[1]);

  if (!ifs) {
    Usage();
  }

  Feature feature;
  vector<Feature> f_vec;
  string line;
  int item;

  while (getline(ifs, line)) {
    if      (line.substr(0, 7) == "WEIGHT:") {
      feature.clear();
      feature.weight = atof(line.substr(7).c_str());
    }
    else if (line.substr(0, 9) == "ITEMSETS:") {
      istrstream is ((char *)line.substr(9).c_str());
      while (is >> item)
        feature.itemset.push_back(item);
      sort(feature.itemset.begin(), feature.itemset.end());
      f_vec.push_back(feature);
    }
  }

  ifstream ifs2(argv[2]);
  if (!ifs2) {
    Usage();
  }

  std::vector<double> ps;
  std::vector<double> ys;
  int n = 0;
  while (getline(ifs2, line)) {
#ifdef DEBUG
    cout << "sample no." << n << " begins" << endl;
#endif
    n++;
    double y = 0.f;
    vector<int> itemset;
    istrstream is (line.c_str());
    is >> y;
    int item;
    while (is >> item) itemset.push_back(item);
    sort(itemset.begin(), itemset.end());

    double p = 0.f;
    for (size_t i = 0; i < f_vec.size(); i++) {
      double weight = f_vec[i].weight;
      std::vector<int> &features = f_vec[i].itemset;
      size_t j;
      for (j = 0; j < features.size(); j++) {
        if (!binary_search(itemset.begin(), itemset.end(), features[j]))
          break;
      }
      if (j == features.size()){
#ifdef DEBUG
        cout << "feature";
        for (size_t i = 0; i < features.size(); i++) {
          cout << features[i] << " ";
        }
        cout << "found weight " << weight << " added" << endl;
#endif
        p += weight;
      }
    }

    //  std::cout << p_val << std::endl;

    ps.push_back(p);
    ys.push_back(y);
  }
  size_t l = ys.size();
  double mean_y = 0.f;
  //  double mean_py = 0.f;
  for (size_t i = 0; i < l; i++) {
    mean_y += ys[i];
    //  mean_py += ps[i];
  }
  mean_y = mean_y/l;
  //  mean_py = mean_py/l;
  double stdy = 0.0;
  //  double stdpy = 0.0;
  for (size_t i = 0; i < l; i++) {
    stdy += (ys[i]-mean_y)*(ys[i]-mean_y);
    //  stdpy += (ps[i]-mean_py)*(ys[i]-mean_py);
  }
  stdy = sqrt(stdy/(l-1));
  //  stdpy = sqrt(stdpy/(l-1));
  // //  for (size_t i = 0; i < l; i++) { // l1nls no-normalize
  // //     ys[i] = (ys[i]-mean_y);//stdy;
  // //     //  ps[i] = (ps[i]-mean_py)/stdpy;
  // //  }

  double rss = 0.f;
  cout << "label\t prediction" << endl;
  for (size_t i = 0; i < l; i++) {
    cout << ys[i] << "\t" << ps[i] << endl;
    rss += (ps[i] - ys[i]) * (ps[i] - ys[i]);
  }
  cout << endl;

  double var = 0.f;
  for (size_t i = 0; i < l; i++) {
    var += (ys[i] - mean_y) * (ys[i] - mean_y);
  }

  double Q2 = 1.f - rss/var;
  std::cout << "RSS:" << rss << " Q2:" << Q2 << std::endl;

}
