#include <iostream>
#include <ostream>
#include <string>
#include <cstring>
#include <vector>
#include "Lcm.h"
#include "Database.h"


/* Globals */
void Usage();
void Version();
void ParseParameters (int argc, char **argv);

vector<string> filenames;
string validation_filename;
bool           toggle_verbose       = false;
int            min_sup              = 1;    // minimum support for LCM
int            max_pat              = 0;    // maximum_itemset_size : default infinity;
unsigned int   topL                 = 1;
unsigned int   maxiter              = 1000;

using namespace std;

/*****************************************************************************************
 * The main function
 *****************************************************************************************/
int main(int argc, char **argv)
{
  Version();

  if(argc < 3) Usage();

  ParseParameters(argc, argv);

  // set dataset
  if (toggle_verbose) {
    cerr << "set database" << endl;
    cerr << "filename:" << filenames[0] << endl;
  }

  Database database, validation_database;
  string input_filename = filenames[0];
  ifstream ifs(input_filename.c_str());
  if(!ifs){
    Usage();
  }else{
    database.ReadFile(input_filename);
  }

  if( !validation_filename.empty() ){
    validation_database.ReadFile(validation_filename);
  }

  string output_filename = filenames[1];
  ofstream ofs(output_filename.c_str());
  if(!ofs){
    Usage();
  }

  // set LCM parameter
  Lcm lcm(cout, min_sup, max_pat);
  lcm.Init(database, validation_database, output_filename, topL, maxiter);
  lcm.Run();

  return 0;
}

/*****************************************************************************
 * Version
 ****************************************************************************/
void Version(){
  cerr << "pLARS - polynomial LARS - version 1.00" << endl << endl;
}

/***************************************************************************
 * Usage
 ***************************************************************************/
void Usage(){
  cerr << endl
       << "Usage: plars [OPTION] TRAIN_FILE MODEL_FILE" << endl << endl
       //       << "       where [OPTION]...  is a list of zero or more optional arguments" << endl
       //       << "             INFILE(s)    is the name of the input transaction database" << endl << endl
       << "OPTION:" << endl
       << "       -min_sup [minimum support]" << endl
       << "       -max_pat [maximum pattern]" << endl
       //       << "       -topL    [l]"               << endl
       << "       -maxiter [maximum iteration]"     << endl
       << "       -validation [VALIDATION_FILE]"     << endl
       << endl;
  exit(0);
}

/*****************************************************************************
 * ParseParameters
 *****************************************************************************/
void ParseParameters (int argc, char **argv){
  if (argc == 1) Usage();
  filenames.clear();

  for (int argno = 1; argno < argc; argno++){
    if (argv[argno][0] == '-'){
      if      (!strcmp (argv[argno], "-version")){
        Version();
      }
      else if (!strcmp (argv[argno], "-verbose")) {
        toggle_verbose = true;
      }
      else if (!strcmp (argv[argno], "-min_sup")) {
        if (argno == argc - 1) cerr << "Must specify minimum support after -min_sup" << endl;
        min_sup = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-max_pat")) {
        if (argno == argc - 1) cerr << "Must specify miximum itemset size after -max_pat" << endl;
        max_pat = atoi(argv[++argno]);
      }
      //  else if (!strcmp (argv[argno], "-topL")) {
      //    if (argno == argc - 1) cerr << "Must specify miximum itemset size after -topL" << endl;
      //    topL = atoi(argv[++argno]);
      //  }
      else if (!strcmp (argv[argno], "-maxiter")) {
        if (argno == argc - 1) cerr << "Must specify miximum itemset size after -maxiter" << endl;
        maxiter = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-validation")) {
        if (argno == argc - 1) cerr << "Must specify validation file name after -validation" << endl;
        validation_filename = argv[++argno];
      }
    }
    else {
      filenames.push_back(argv[argno]);
    }
  }
}
