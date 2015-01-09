#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <algorithm>
#include <vector>
#include <gsl/gsl_multifit.h>
#include <cassert>

using namespace std;

vector<int> mate_one_chr(vector<int>& p1, vector<int>& p2, 
			 default_random_engine& generator,
			 int first_pos, int last_pos, vector<int>& snp_list)
{
  // from mancera et al 2008 nature
  double extra_crossover_prob = 6.1 * (double)(last_pos - first_pos + 1) / pow(10, 6);
OB
  poisson_distribution<int> p_dist(extra_crossover_prob);
  int num_crossovers = p_dist(generator);

  vector<int> crossovers;

  for(int i = 0; i < num_crossovers + 1; i++)
    {
      int site = (int)(round(rand() % (last_pos - first_pos + 1))) + first_pos;
      // Returns an iterator pointing to the first element in the
      // range which does not compare less than val.
      int site_index = lower_bound(snp_list.begin(), snp_list.end(), site) - snp_list.begin();
      assert(site_index <= snp_list.size());
      assert(site_index >= 0);
      if(site_index == snp_list.size())
	site_index -= 1;
      // don't know why this was here:
      //if(site - snp_list[site_index] > snp_list[site_index + 1] - site_index)
      //site_index = site_index + 1;
      if(((double)rand()) / RAND_MAX < 0.5)
	crossovers.push_back(site_index);
    }
  sort(crossovers.begin(), crossovers.end());
  vector<int> gamete;
  vector<int>::iterator next_crossover = crossovers.begin();
  // need to tile p1 and p2, alternating based on indices in crossovers
  vector<int>::iterator it1 = p1.begin();
  vector<int>::iterator it2 = p2.begin();
  bool toggle = true;
  if(((double)rand()) / RAND_MAX < 0.5)
    toggle = false;
  assert(p1.size() == p2.size());
  for(int i = 0; i < p1.size(); i++)
    {
      if(toggle)
	gamete.push_back(*it1);
      else
	gamete.push_back(*it2);
      if(next_crossover != crossovers.end() && i == *next_crossover)
	{
	  toggle = !toggle;
	  next_crossover++;
	}
      it1++;
      it2++;
    }
  assert(gamete.size() == p1.size());
  return gamete;
}

vector<int> mate(vector<int>& p1, vector<int>& p2, 
		 default_random_engine& generator, 
		 int first_pos, int last_pos, 
		 map<int, vector<int> >& snp_list)
{
  map<int, vector<int> >::iterator it = snp_list.begin()

  int start_ind = 0;
  int end_ind = snp_list.begin()->size();
  vector<int> gamete;
  while (it != snp_list.end())
  {
    vector<int> p1_chr = vector<int>(p1.begin() + start_ind, p1.begin() + end_ind);
    vector<int> p2_chr = vector<int>(p2.begin() + start_ind, p2.begin() + end_ind);
      
    gamete_chr = mate(p1_chr, p2_chr, generator, first_pos, last_pos, it->second);
    gamete.insert(gamete_chr.begin(), gamete_chr.size());

    start_ind = end_ind;
    end_ind += it->size();
    it++;
  }
  return gamete;
}

vector<int> gen_segregant(vector<vector<int> >& initial_generation, 
			  default_random_engine& generator,			  
			  map<int, vector<int> >& snp_list, 
			  vector<int> first_pos, vector<int> last_pos)
{
  int num_gens = (int)(log(initial_generation.size())/log(2));
  // for each generation, have collection of organisms, each of which is a vector
  vector<vector<vector<int> > > generations;
  vector<vector<int> > previous_generation = initial_generation;
  generations.push_back(initial_generation);
  for(int n = 0; n < num_gens; n++)
    {
      vector<vector<int> > new_generation;
      previous_generation = generations.back();
      random_shuffle(previous_generation.begin(), previous_generation.end());
      vector<vector<int> >::iterator it = previous_generation.begin();
      for(int i = 0; i < previous_generation.size(); i += 2)
	{
	  vector<int> p1 = *it;
	  it++;
	  vector<int> p2 = *it;
	  vector<int> baby = mate(p1, p2, generator, first_pos, last_pos, snp_list);
	  it++;
	  new_generation.push_back(baby);
	}
      generations.push_back(new_generation);
    }
  return generations.back().back();
}

vector<vector<int> > gen_segregants(vector<vector<int> >& initial_generation, 
				    default_random_engine& generator,
				    vector<int>& chr, vector<int>& ps,
				    vector<double>& af,
				    vector<int>& first_pos, vector<int>& last_pos,
				    int num_segregants, int num_per_cross = 2 * num_segregants)
{
  
}

vector<vector<int> > gen_segregants(vector<vector<int> >& initial_generation, 
				    default_random_engine& generator,
				    vector<int>& chr, vector<int>& ps,
				    vector<double>& af,
				    vector<int>& first_pos, vector<int>& last_pos,
				    int num_segregants, int num_per_cross = 2 * num_segregants)
{
  vector<vector<int> > segregants(num_segregants, vector<int>());
  for(int i = 0; i < num_segregants; i++)
    {
      if(i%100 == 0)
	cout << "seg " << i << '\n' << flush;
      segregants[i] = gen_segregant(initial_generation, generator, snp_list, 
				    first_pos, last_pos);

    }
  return segregants;
}

vector<double> sim_phenotypes(vector<vector<int> >& segregants, int qtl1, int qtl2,
			      double fixed_effect, default_random_engine& generator,
			      normal_distribution<double>& n_dist_0,
			      normal_distribution<double>& n_dist_effect)
{
  
  vector<double> phenotypes;

  int epi_allele_1;
  int epi_allele_2;
  if(((double)rand()) / RAND_MAX > 0.5)
    epi_allele_1 = 0;
  else
    epi_allele_1 = 1;
  if(((double)rand()) / RAND_MAX > 0.5)
    epi_allele_2 = 0;
  else
    epi_allele_2 = 1;

  for(int i = 0; i < segregants.size(); i++)
    {
      int allele1 = segregants[i][qtl1];
      int allele2 = segregants[i][qtl2];
      if(allele1 == epi_allele_1 && allele2 == epi_allele_2)
	phenotypes.push_back(n_dist_effect(generator));
      else
	phenotypes.push_back(n_dist_0(generator));
    }
  return phenotypes;
}

double predict(vector<double>& phenotypes, vector<vector<int> >& segregants, int& max_LOD_ind_1,
	       int& max_LOD_ind_2, double& h2_pred, vector<vector<double> >& LOD_scores)
{
  int n = phenotypes.size();
  int num_sites = segregants[0].size();
  for(int i = 0; i < num_sites; i++)
    {
      vector<double> row(num_sites);
      for(int j = 0; j < num_sites; j++)
	row[j] = 0.0;
      LOD_scores.push_back(row);
    }
  double max_LOD = 0;
  int skip = 1000;

  gsl_vector  * Y = gsl_vector_alloc(n);
  for (int i = 0; i < n; i++)
    {
      gsl_vector_set(Y, i, phenotypes[i]);
    }

  // for no interaction model
  gsl_matrix * X = gsl_matrix_alloc(n, 3);
  gsl_vector * c = gsl_vector_alloc (3);
  gsl_matrix * cov = gsl_matrix_alloc (3, 3);
  double chisq;

  // for interaction model
  gsl_matrix * Xi = gsl_matrix_alloc(n, 4);
  gsl_vector * ci = gsl_vector_alloc (4);
  gsl_matrix * covi = gsl_matrix_alloc (4, 4);
  double chisqi;

  assert(num_sites > 1);
  int success_count = 0;
  for(int i = 0; i < num_sites; i++)
    {
      if(i % 100 == 0)
	cout << '*' << i << '\n' << flush;
      for(int j = i + skip; j < num_sites; j++)
	{
	  // add a row for each segregant, at the same time checking
	  // whether we see every possible combination of the two - if
	  // not, we won't bother trying to fit this model
	  bool seen00 = false;
	  bool seen01 = false;
	  bool seen10 = false;
	  bool seen11 = false;

	  for(int s = 0; s < n; s++)
	    {
	      int allele1 = segregants[s][i];
	      int allele2 = segregants[s][j];

	      // no interaction
	      gsl_matrix_set(X, s, 0, 1.0);
	      gsl_matrix_set(X, s, 1, allele1);
	      gsl_matrix_set(X, s, 2, allele2);

	      // interaction
	      gsl_matrix_set(Xi, s, 0, 1.0);
	      gsl_matrix_set(Xi, s, 1, allele1);
	      gsl_matrix_set(Xi, s, 2, allele2);
	      if(allele1 == 0)
		{
		  if(allele2 == 0) // 00
		    {
		      gsl_matrix_set(Xi, s, 3, 1);
		      seen00 = true;
		    }
		  else // 01
		    {
		      gsl_matrix_set(Xi, s, 3, 0);
		      seen01 = true;
		    }
		}
	      else
		{
		  if(allele2 == 0) // 10
		    {
		      gsl_matrix_set(Xi, s, 3, 0);
		      seen10 = true;
		    }
		  else // 11
		    {
		      gsl_matrix_set(Xi, s, 3, 0);
		      seen11 = true;
		    }
		}
	    }

	  if(seen00 && seen01 && seen10 && seen11)
	    {
	      success_count++;

	      // fit no interaction model
	      gsl_multifit_linear_workspace * work  = gsl_multifit_linear_alloc (n, 3);
	      gsl_multifit_linear(X, Y, c, cov, &chisq, work);
	      gsl_multifit_linear_free(work);

	      // fit interaction model
	      gsl_multifit_linear_workspace * worki  = gsl_multifit_linear_alloc (n, 4);
	      gsl_multifit_linear(Xi, Y, ci, covi, &chisqi, worki);
	      gsl_multifit_linear_free(worki);

#define C(i) (gsl_vector_get(c,(i)))
#define Ci(i) (gsl_vector_get(ci,(i)))
	      {
		for(int k = 0; k < 4; k++)
		  if(Ci(k) > 100)
		    {
		      cout << "LARGE\n";
		      cout << "coefficients: " <<  C(0) << ' ' << C(1) << ' ' <<C(2) << ' ' <<C(3) << '\n';
		      cout << flush;
		    }

		// calculate RSS for both models
		double RSS = 0.0;
		double fitted_Y;

		double RSSi = 0.0;
		double fitted_Yi;

		for(int k = 0; k < n; k++)
		  {
		    fitted_Y = C(0) * gsl_matrix_get(X, k, 0) +
		      C(1) * gsl_matrix_get(X, k, 1) +
		      C(2) * gsl_matrix_get(X, k, 2);
		    RSS += pow(gsl_vector_get(Y, k) - fitted_Y, 2);

		    fitted_Yi = Ci(0) * gsl_matrix_get(Xi, k, 0) +
		      Ci(1) * gsl_matrix_get(Xi, k, 1) +
		      Ci(2) * gsl_matrix_get(Xi, k, 2) +
		      Ci(3) * gsl_matrix_get(Xi, k, 3);
		    RSSi += pow(gsl_vector_get(Y, k) - fitted_Yi, 2);
		  }

		// and finally F and LOD score
		double p = 3;
		double pi = 4;
		double F = ((RSS - RSSi) / (pi - p)) / (RSSi / (n - pi));
		int df = pi - p; // is this correct?
		double LOD = n / 2.0 * log((F * df) / (n - df - 1) + 1) / log(10);
		if(LOD > max_LOD)
		  {
		    max_LOD = LOD;
		    max_LOD_ind_1 = i;
		    max_LOD_ind_2 = j;
		  }
	      }
	    }
	  //else
	  // {
	      
	  //}
	}
    }
  gsl_vector_free(Y);

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(c);

  gsl_matrix_free(Xi);
  gsl_matrix_free(covi);
  gsl_vector_free(ci);

  h2_pred = 1 - pow(10, -2.0 / n * max_LOD);

  cout << "number sites with all allele combos: " << success_count << '\n' << flush;

  return max_LOD;
}

void write_output_header(ofstream& f_summary)
{
  f_summary << "pred_qtl1 pred_qtl2 qtl1 qtl2 LOD_score h2_prediction allele_frequency1 allele_frequency2\n";
  f_summary.flush();
}

void write_stats(ofstream& f_summary, ofstream& f_probs, int qtl_pred_1, int qtl_pred_2, 
		 int qtl1, int qtl2, double max_LOD, double h2_pred,
		 double allele_freq_1, double allele_freq_2, vector<vector<double> >& LOD_scores)
{
  f_summary << qtl_pred_1 << ' ' << qtl_pred_2 << ' ' << qtl1 << ' ' << qtl2 << ' ' 
	    << max_LOD << ' ' << h2_pred << ' ' << allele_freq_1 << ' ' << allele_freq_2 <<  '\n';
  //int n = LOD_scores.size();
  //for(int i = 0; i < n; i++)
  //  for(int j = 0; j < n; j++)
  //    f_probs << LOD_scores[i][j] << ',';
  f_summary << flush;
  //f_probs << '\n' << flush;
}

vector<double> heritability_to_fixed_effect(vector<double> h2, 
					    vector<double> af, int n)
{
  vector<double> k;
  for (int i = 0; i < af.size(); i++)
    k.push_back(pow(h2[i] * (n-1) / ((h2[i]-1) * n * (af[i]-1) * af[i]), .5));
  return k;
}

vector<vector<int> > get_parent_gen(string fn, vector<int> chr, 
				    vector<int> ps, vector<double>& af)
{
  // file format is 
  // ,chr,chr
  // ,ps,ps
  // ,af,af
  // strain,allele,allele
  // strain,allele,allele

  ifstream f;
  f.open(fn);
  string line;

  // read in chromosomes
  getline(f, line);
  int index = line.find(',');
  int new_index;

  while(index < line.length())
    {
      new_index = line.find(',', index + 1);
      chr.push_back(atoi(line.substr(index + 1, new_index - index - 1).c_str()));
      index = new_index;
    }

  // read in positions
  getline(f, line);
  index = line.find(',');

  while(index < line.length())
    {
      new_index = line.find(',', index + 1);
      ps.push_back(atoi(line.substr(index + 1, new_index - index - 1).c_str()));
      index = new_index;
    }
    
  // read in allele frequencies
  getline(f, line);
  index = line.find(',');
  while(index < line.length())
    {
      new_index = line.find(',', index + 1);
      af.push_back(atof(line.substr(index + 1, new_index - index - 1).c_str()));
      index = new_index;
    }
  
  // read in genotypes
  vector<vector<int> > parent_gen;
  while(getline(f, line))
    {
      index = line.find(',');
      vector<int> parent;
      while(index < line.length())
	{
	  new_index = line.find(',', index + 1);
	  parent.push_back(atoi(line.substr(index + 1, new_index - index - 1).c_str()));
	  index = new_index;
	}
      parent_gen.push_back(parent);
    }

  return parent_gen;
}

void sim(int num_sims, int num_segregants, string outfile_ext, 
	 vector<double> h2, 
	 string parent_file, vector<int> first_pos, vector<int> last_pos,
	 int num_qtls)
{
  // print parameters for reference
  cout << num_sims << ", " << num_segregants << ", " << outfile_ext
       << ", " << h2 << ", " << parent_file << ", " << first_pos
       << ", " <<  last_pos << '\n' << flush;

  // read in parent generation genotypes and information about the
  // variant sites
  vector<int> chr;
  vector<int> ps;
  vector<double> af;
  vector<vector<int> > parents = get_parent_gen(parent_file, chr, ps, af);
  int num_snps = ps.size();
  vector<int> unique_chr; // all the unique chromosomes in the data set
  int prev_chr = -1;
  for (int i = 0; i < num_snps; i++) 
    {
      if (chr[i] != prev_chr) 
	{
	  unique_chr.push_back(chr[i]);
	  prev_chr = chr[i];
	}
    }
  assert(num_qtls <= unique_chr.size())
  
  // open output file and write column headers
  string summary_fn = "out/sim_output_summary" + outfile_ext + ".txt";
  ofstream f_summary;
  f_summary.open(summary_fn);
  write_output_header(f_summary);

  // initialize random seed (note that time resolution (1 second?) can
  // be an issue here if we're starting a bunch of simulations at
  // close to the same time)
  srand(time(NULL));

  // random number generator and distributions for the two alleles
  // (one is centered at 0 and one at a fixed effect to be determined
  // later based on the allele frequencies of the snvs we choose)
  default_random_engine generator;
  normal_distribution<double> n_dist_0(0, 1);
  normal_distribution<double> n_dist_effect;

  // vectors for storing segregants, their phenotypes, and the qtls
  // and their fixed effects for each simulation
  vector<vector<int> > segregants;
  vector<double> phenotypes;
  vector<int> qtls;
  vector<double> fixed_effects;

  // run simulations
  for(int s = 0; s < num_sims; s++)
    {
      cout << s << '\n';
      cout.flush();

      // generate segregants by simulating mating of parents in a
      // funnel cross
      segregants = gen_segregants(parents, generator, chr, ps, 
				  num_segregants, first_pos, last_pos);

      // choose qtls (all on different chromosomes)
      int num_sites_left = num_sites;
      for(int q = 0; q < num_qtls; q++)
	{
	  int index = (int)(round(rand() % num_sites_left));

	  // remove all sites on same chromosome before picking again
	}
      

      fixed_effects = heritability_to_fixed_effect(h2, af, num_segregants);

      phenotypes = sim_phenotypes(segregants, qtl1, qtl2, fixed_effect,
				  generator, n_dist_0, n_dist_effect);
      
      int qtl_pred_1 = 0;
      int qtl_pred_2 = 0;
      double h2_pred;
      vector<vector<double> > LOD_scores;
      double max_LOD = predict(phenotypes, segregants, qtl_pred_1, qtl_pred_2, h2_pred, LOD_scores);

      // store all the results
      write_stats(f_summary, f_probs, qtl_pred_1, qtl_pred_2, qtl1, qtl2, max_LOD, h2_pred,
		  allele_freqs[qtl1], allele_freqs[qtl2], LOD_scores);
    }
  f_summary.close();
  f_probs.close();
}

