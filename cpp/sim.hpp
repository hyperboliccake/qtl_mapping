#ifndef SIM_INCLUDED
#define SIM_INCLUDED 1

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <gsl/gsl_multifit.h>


std::vector<int> mate_one_chr(std::vector<int>& p1, std::vector<int>& p2, 
			      std::default_random_engine& generator,
			      int first_pos, int last_pos, 
			      std::vector<int>& snp_list);

std::vector<int> mate(std::vector<int>& p1, std::vector<int>& p2, 
		      std::default_random_engine& generator,
		      int first_pos, int last_pos, 
		      std::map<int, std::vector<int> >& snp_list);

std::vector<int> gen_segregant(std::vector<std::vector<int> >& initial_generation, 
			      std::default_random_engine& generator,
			       std::map<int, std::vector<int> >& snp_list, 
			       int first_pos, int last_pos);

std::vector<std::vector<int> > gen_segregants(std::vector<std::vector<int> >& initial_generation, 
					      std::default_random_engine& generator,
					      std::map<int, std::vector<int>>& snp_list, 
					      int num_segregants, 
					      int first_pos, int last_pos);

std::vector<double> sim_phenotypes(std::vector<std::vector<int> >& segregants, 
				   int qtl1, int qtl2,
				   double fixed_effect, 
				   std::default_random_engine& generator,
				   std::normal_distribution<double>& n_dist_0,
				   std::normal_distribution<double>& n_dist_effect);

double predict(std::vector<double>& phenotypes, 
	       std::vector<std::vector<int> >& segregants,
	       int& max_LOD_ind_1, int& max_LOD_ind_2, double& h2_pred,
	       std::vector<std::vector<double> >& LOD_scores);

void write_stats_headers(std::map<int, std::vector<int> >& snp_list, 
			 std::ofstream& f_summary, std::ofstream& f_probs);

void write_stats(std::ofstream& f_summary, std::ofstream& f_probs, 
		 int qtl_pred_1, int qtl_pred_2, 
		 int qtl1, int qtl2, double max_LOD, double h2_pred,
		 double allele_freq_1, double allele_freq_2, 
		 std::vector<std::vector<double> >& LOD_scores);

double heritability_to_fixed_effect(double h2);

std::vector<std::vector<int> > get_parent_gen(std::string fn, 
					      std::map<int, std::vector<int> >& snp_list, 
					      std::vector<double>& allele_freqs);


void sim(int num_sims = 1, int num_segregants = 100, 
	 std::string outfile_ext = "", double h2 = .1,
	 std::string parent_file = "parent_generation_for_sim.txt", 
	 vector<int> first_pos = {1, 1}, vector<int> last_pos = {252867, ...});


#endif
