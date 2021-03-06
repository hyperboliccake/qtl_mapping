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
			      std::vector<int>& ps);

std::vector<std::vector<int> > mate(std::vector<std::vector<int> >& p1, 
				     std::vector<std::vector<int> >& p2, 
				     std::default_random_engine& generator,
				     std::vector<int>& first_pos, 
				     std::vector<int>& last_pos, 
				     std::vector<std::vector<int> >& ps);

std::vector<std::vector<std::vector<int> > > gen_segregants(std::vector<std::vector<std::vector<int> > >& initial_generation, 
			       std::default_random_engine& generator,
			       std::vector<int>& chr,
			       std::vector<std::vector<int> >& ps,
			       std::vector<int>& first_pos, std::vector<int>& last_pos, 
			       int num_segregants, int num_per_cross);

std::vector<double> sim_phenotypes(std::vector<std::vector<std::vector<int> > >& segregants, 
				   std::vector<std::pair<int, int> > qtls,
				   double fixed_effect, 
				   std::default_random_engine& generator);

double predict_one(std::vector<double>& phenotypes, 
		   std::vector<std::vector<std::vector<int> > >& segregants,
		   std::vector<std::vector<int> >& ps,
		   std::pair<int, int>& max_LOD_ind,
		   double& h2_pred,
		   std::vector<std::vector<double> >& LOD_scores);

double predict_two(std::vector<double>& phenotypes, 
		   std::vector<std::vector<std::vector<int> > >& segregants,
		   std::vector<std::vector<int> >& ps,
		   std::vector<std::pair<int, int> >& max_LOD_inds,
		   double& h2_pred,
		   std::vector<double>& LOD_scores,
		   std::vector<std::vector<int> >& LOD_scores_loc,
		   std::vector<std::pair<int, int> > ranges,
		   std::vector<int> correct_chrs);

void write_output_header(std::ofstream& f_summary, int num_qtls);

void write_LOD_header(std::ofstream& f_LOD, int num_sims, 
		      std::vector<int>& chr, std::vector<std::vector<int> >& ps);

void write_stats(std::ofstream& f_summary,
		 std::ofstream& f_LOD,
		 std::vector<int>& chr,
		 std::vector<std::vector<int> >& ps,
		 std::vector<std::pair<int, int> > qtl_pred,
		 std::vector<std::pair<int, int> > qtl,
		 double max_LOD, double h2_pred,
		 std::vector<double> af,
		 std::vector<std::vector<double> > & LOD_scores, int sim);

void write_stats_2(std::ofstream& f_summary,
		 std::ofstream& f_LOD,
		 std::vector<int>& chr,
		 std::vector<std::vector<int> >& ps,
		 std::vector<std::pair<int, int> > qtl_pred,
		 std::vector<std::pair<int, int> > qtl,
		 double max_LOD, double h2_pred,
		 std::vector<double> af,
		 std::vector<double>& LOD_scores, 
		   std::vector<std::vector<int> > LOD_scores_loc,
		   int sim);

double heritability_to_fixed_effect(double h2,
				    std::vector<double> af,
				    int n);

std::vector<std::vector<std::vector<int> > > get_parent_gen(std::string fn, 
							    std::vector<int>& chr,
							    std::vector<std::vector<int> >& ps,
							    std::vector<std::vector<double> >& af);


void sim(int num_sims, int num_segregants, int num_per_cross,
	 std::string outfile_ext, int num_qtls, double h2,
	 std::string parent_file,
	 std::vector<int> first_pos, std::vector<int> last_pos);


#endif
