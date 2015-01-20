#include "sim.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, const char* argv[])
{
  //sim(atoi(argv[1]), atoi(argv[2]), argv[3]);
  vector<int> first_pos;
  vector<int> last_pos;
  ifstream f;
  f.open("first_last_pos_all_chr_sgd.txt");
  string line;
  while(getline(f, line))
    {
      int start_index = line.find(',') + 1;
      int end_index = line.find(',', start_index);
      first_pos.push_back(atoi(line.substr(start_index, end_index - start_index).c_str()));
      start_index = end_index + 1;
      end_index = line.length();
      last_pos.push_back(atoi(line.substr(start_index, end_index - start_index).c_str()));      
    }

  int num_sims = atoi(argv[1]);
  int num_segregants = atoi(argv[2]);
  int num_per_cross = atoi(argv[3]);

  int num_qtl = atoi(argv[4]);
  vector<double> h2;
  for (int q = 0; q < num_qtl; q++)
    h2.push_back(atof(argv[5 + q]));

  string parent_file = "/net/gs/vol1/home/aclark4/qtl_mapping/grant/code/diversity/8_strains_binary_test.txt";

  int run_id = atoi(argv[5 + num_qtl]);

  stringstream ss;
  ss << '_' << num_sims;
  ss << '_' << num_segregants;
  ss << '_' << num_per_cross;
  for (int q = 0; q < num_qtl; q++)
    ss << '_' << h2[q];
  string outfile_ext = ss.str();

  cout << "Num sims: " << num_sims << '\n';
  cout << "Num segregants: " << num_segregants << '\n';
  cout << "Num per cross: " << num_per_cross << '\n';
  for (int q = 0; q < num_qtl; q++)
    cout << "QTL heritability: " << h2[q] << '\n';
  cout << "Parent file: " << parent_file << '\n';
  cout << "Run ID: " << run_id << '\n' << flush;

  sim(num_sims, num_segregants, num_per_cross, outfile_ext, h2, 
      parent_file, first_pos, last_pos);
  
  cout << "Done!\n" << flush;

  return 0;
}
