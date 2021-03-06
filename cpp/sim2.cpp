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
			 int first_pos, int last_pos, 
			 vector<int>& ps)
{
  // from mancera et al 2008 nature
  double extra_crossover_prob = 6.1 * (double)(last_pos - first_pos + 1) / pow(10, 6);

  poisson_distribution<int> p_dist(extra_crossover_prob);
  int num_crossovers = p_dist(generator);

  vector<int> crossovers;
  // add 1 because we're guaranteed one crossover per chromatid pair
  for(int i = 0; i < num_crossovers + 1; i++)
    {
      int site = (int)(round(rand() % (last_pos - first_pos + 1))) + first_pos;
      // Returns an iterator pointing to the first element in the
      // range which does not compare less than val
      int site_index = lower_bound(ps.begin(), ps.end(), site) - ps.begin();
      assert(site_index <= ps.size());
      assert(site_index >= 0);
      if(site_index == ps.size())
	site_index -= 1;
      // 50% chance of ending up on chromatid pair we choose
      if(((double)rand()) / RAND_MAX < 0.5)
	crossovers.push_back(site_index);
    }
  sort(crossovers.begin(), crossovers.end());
  vector<int> gamete;
  vector<int>::iterator next_crossover = crossovers.begin();
  // need to tile p1 and p2, alternating based on indices in crossovers
  vector<int>::iterator it1 = p1.begin();
  vector<int>::iterator it2 = p2.begin();
  // toggle for which parent we're currently on
  bool toggle = true; 
  // choose which parent we start with
  if(((double)rand()) / RAND_MAX < 0.5)
    toggle = false;
  assert(p1.size() == p2.size());
  // then alternate between parents
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

vector<vector<int> > mate(vector<vector<int> >& p1, vector<vector<int> >& p2, 
			  default_random_engine& generator, 
			  vector<int>& first_pos, vector<int>& last_pos, 
			  vector<int>& chr, vector<vector<int> >& ps)
{
  assert(p1.size() == p2.size());
  assert(chr.size() == p1.size());
  assert(chr.size() == ps.size());

  vector<vector<int> > gamete;
  vector<int> gamete_chr;
  int num_chr = chr.size();

  for (int i = 0; i < num_chr; i++)
    {
      gamete_chr = mate_one_chr(p1[i], p2[i], generator, first_pos[i], 
				last_pos[i], ps[i]);
      gamete.push_back(gamete_chr);
    }
  
  return gamete;
}

vector<vector<vector<int> > > gen_segregants(vector<vector<vector<int> > >& initial_generation, 
					     default_random_engine& generator,
					     vector<int>& chr, 
					     vector<vector<int> >& ps,
					     vector<int>& first_pos, 
					     vector<int>& last_pos,
					     int num_segregants, int num_per_cross)
{
  // check that size of parent generation is a power of two
  assert((initial_generation.size() & (initial_generation.size() - 1)) == 0);

  // this is all hardcoded for cross of 8 strains right now
  assert(initial_generation.size() == 8);

  // mate pairs in order, generating a pool of num_per_cross offspring
  // for each pair
  vector<vector<vector<vector<int> > > > pools;
  for(int i = 0; i < initial_generation.size(); i += 2)
    {
      vector<vector<vector<int> > > pool;
      for(int s = 0; s < num_per_cross; s++) 
	{
	  vector<vector<int> > baby = mate(initial_generation[i], 
					   initial_generation[i + 1], 
					   generator, first_pos, last_pos, chr, ps);
	  pool.push_back(baby);
	}
      pools.push_back(pool);
    }

  // then mate those pools by randomly selecting one mate from one and
  // the other from the other
  vector<vector<vector<vector<int> > > > new_pools;
  for (int i = 0; i < pools.size(); i += 2)
    {
      // choose from each of the two pools with replacement; in
      // reality we would do this twice switching which is a and
      // alpha, but that doesn't make a difference here
      vector<vector<vector<int> > > pool1 = pools[i];
      vector<vector<vector<int> > >pool2 = pools[i + 1];
      assert(pool1.size() == num_per_cross);
      assert(pool2.size() == num_per_cross);

      vector<vector<vector<int> > > new_pool;
      for (int j = 0; j < num_per_cross; j++)
	{
	  vector<vector<int> > p1 = pool1[rand() % num_per_cross];
	  vector<vector<int> > p2 = pool2[rand() % num_per_cross];
	  new_pool.push_back(mate(p1, p2, generator, first_pos, last_pos, chr, ps));
	}
      new_pools.push_back(new_pool);
    }

  // then take the two pools from above and mate them to produce
  // num_segregants offspring
  vector<vector<vector<int> > > final_pool;
  vector<vector<vector<int> > > pool1 = new_pools[0];
  vector<vector<vector<int> > > pool2 = new_pools[1];
  assert(pool1.size() == num_per_cross);
  assert(pool2.size() == num_per_cross);
  for (int j = 0; j < num_segregants; j++)
    {
      vector<vector<int> > p1 = pool1[rand() % num_per_cross];
      vector<vector<int> > p2 = pool2[rand() % num_per_cross];
      final_pool.push_back(mate(p1, p2, generator, first_pos, last_pos, chr, ps));
    }

  return final_pool;
}

vector<double> sim_phenotypes(vector<vector<vector<int> > >& segregants, 
			      vector<pair<int, int> > qtls,
			      double fixed_effect, 
			      default_random_engine& generator)
{
  // model is that you have an effect only if you have all of the
  // qtls, no effect otherwise

  normal_distribution<double> n_dist_0(0, 1);
  vector<double> phenotypes;

  cout << "fixed effect " << fixed_effect << '\n' << flush;
  normal_distribution<double> n_dist_effect(fixed_effect, 1);

  for(int i = 0; i < segregants.size(); i++)
    {
      bool all_minor = true;
      for (int q = 0; q < qtls.size(); q++)
	{
	  int chr_ind = qtls[q].first;
	  int ps_ind = qtls[q].second;
	  // minor allele is always 0, major 1
	  if (segregants[i][chr_ind][ps_ind] == 1)
	    {
	      all_minor = false;
	      break;
	    }
	}

      if (all_minor)
	{
	  phenotypes.push_back(n_dist_effect(generator));
	}
      else
	{
	  phenotypes.push_back(n_dist_0(generator));
	}
    }
  return phenotypes;
}

double predict_one(vector<double>& phenotypes, 
		   vector<vector<vector<int> > >& segregants, 
		   vector<vector<int> >& ps,
		   pair<int, int>& max_LOD_ind,
		   double& h2_pred, 
		   vector<vector<double> >& LOD_scores)
{
  int n = phenotypes.size();
  int num_chr = segregants[0].size();
  int num_sites = 0;
  LOD_scores = vector<vector<double> >();
  for(int i = 0; i < num_chr; i++)
    {
      num_sites += ps[i].size();
      LOD_scores.push_back(vector<double>(ps[i].size(), 0));
    }
  double max_LOD = -numeric_limits<double>::infinity();
  max_LOD_ind = pair<int, int>(-1, -1);
  
  gsl_vector  * Y = gsl_vector_alloc(n);
  for (int i = 0; i < n; i++)
    {
      gsl_vector_set(Y, i, phenotypes[i]);
    }

  // association model
  // intercept + allele
  gsl_matrix * X = gsl_matrix_alloc(n, 2); 
  // coefficients for above
  gsl_vector * c = gsl_vector_alloc (2);
  // covariance matrix for coefficients
  gsl_matrix * cov = gsl_matrix_alloc (2, 2);
  double chisq;

  double fitted_Y_null = accumulate(phenotypes.begin(), phenotypes.end(), 0.0) / n;
  double RSS_null = 0;
  for (int k = 0; k < n; k++)
    {
      RSS_null += pow(gsl_vector_get(Y, k) - fitted_Y_null, 2);
      gsl_matrix_set(X, k, 0, 1.0);
    }

  for (int chr_ind = 0; chr_ind < num_chr; chr_ind++)
    {
      int num_ps = ps[chr_ind].size();
      for (int ps_ind = 0; ps_ind < num_ps; ps_ind++)
	{
	  // set alleles for all strains
	  for (int s = 0; s < n; s++)
	    {
	      int allele = segregants[s][chr_ind][ps_ind];
	      gsl_matrix_set(X, s, 1, allele);
	    } 

	  // test for association
	  gsl_multifit_linear_workspace * work  = gsl_multifit_linear_alloc (n, 2);
	  gsl_multifit_linear(X, Y, c, cov, &chisq, work);
	  gsl_multifit_linear_free(work);

#define C(i) (gsl_vector_get(c,(i)))
	  {
	    // calculate RSS for association mode
	    double RSS = 0.0;
	    double fitted_Y;
			
	    for(int k = 0; k < n; k++)
	      {
		fitted_Y = C(0) * gsl_matrix_get(X, k, 0) +
		  C(1) * gsl_matrix_get(X, k, 1);
		RSS += pow(gsl_vector_get(Y, k) - fitted_Y, 2);
	      }

	    // and finally F and LOD score
	    double p_null = 1;
	    double p = 2;
	    double F = ((RSS_null - RSS) / (p - p_null)) / (RSS / (n - p));
	    int df = p - p_null; // is this correct?
	    double LOD = n / 2.0 * log((F * df) / (double)(n - df - 1) + 1) / log(10);

	    // check whether this is the strongest
	    // associaiton we've found so far
	    if(LOD > max_LOD)
	      {
		max_LOD = LOD;
		max_LOD_ind = pair<int, int>(chr_ind, ps_ind);
	      }

	    // and save it no matter how good it is
	    LOD_scores[chr_ind][ps_ind] = LOD;
	  }
	} 
    }

  gsl_vector_free(Y);

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(c);

  h2_pred = 1 - pow(10, -2.0 / n * max_LOD);

  return max_LOD;
}

double predict_two(vector<double>& phenotypes, 
		   vector<vector<vector<int> > >& segregants, 
		   vector<vector<int> >& ps,
		   vector<pair<int, int> >& max_LOD_inds,
		   double& h2_pred, 
		   vector<double>& LOD_scores,
		   vector<vector<int> >& LOD_scores_loc,
		   vector<pair<int, int> > ranges, vector<int> correct_chrs)
{
  // simplest strategy is to test every combination of qtls on
  // different chromosomes (but this will explode very quickly and
  // might not even be reasonable for two loci)

  int n = phenotypes.size();
  int num_chr = segregants[0].size();

  /*
  for (int i = 0; i < 2; i++) 
    {
      vector<vector<double> > qtl_scores;
      for (int c = 0; c < num_chr; c++)
	{
	  qtl_scores.push_back(vector<double>(ps[c].size(), 0));
	}
      LOD_scores.push_back(qtl_scores);
    }
  */

  double max_LOD = 0;
  pair<int, int> max_LOD_ind_1(-1, -1);
  pair<int, int> max_LOD_ind_2(-1, -1);
  

  gsl_vector  * Y = gsl_vector_alloc(n);
  for (int i = 0; i < n; i++)
    {
      gsl_vector_set(Y, i, phenotypes[i]);
    }

  // for no interaction model
  // intercept + each of 2 individual alleles
  gsl_matrix * X = gsl_matrix_alloc(n, 3); 
  // coefficients for above
  gsl_vector * c = gsl_vector_alloc (3);
  // covariance matrix for coefficients
  gsl_matrix * cov = gsl_matrix_alloc (3, 3);
  double chisq;

  // for interaction model
  // same as above model but with an extra interaction term
  gsl_matrix * Xi = gsl_matrix_alloc(n, 4);
  gsl_vector * ci = gsl_vector_alloc (4);
  gsl_matrix * covi = gsl_matrix_alloc (4, 4);
  double chisqi;

  int success_count = 0;
  //for (int chr1 = 0; chr1 < num_chr; chr1++)
  //{
  //  for (int chr2 = chr1 + 1; chr2 < num_chr; chr2++)
  //	{
	  //int num_ps1 = ps[chr1].size();
	  //int num_ps2 = ps[chr2].size();
  int chr1 = correct_chrs[0];
  int chr2 = correct_chrs[1];
	  for (int ps1 = ranges[0].first; ps1 < ranges[0].second; ps1++)
	    {
	      if (ps1 % 100 == 0)
		{
		  cout << ps1 << '\n' << flush;
		}
	      
	      for (int ps2 = ranges[1].first; ps2 < ranges[1].second; ps2++)
		{
		  // add a row for each segregant, at the same time checking
		  // whether we see every possible combination of the two - if
		  // not, we won't bother trying to fit this model
		  bool seen00 = false;
		  bool seen01 = false;
		  bool seen10 = false;
		  bool seen11 = false;

		  // iteratate over all strains
		  for (int s = 0; s < n; s++)
		    {
		      // alleles at the two loci we're testing
		      int allele1 = segregants[s][chr1][ps1];
		      int allele2 = segregants[s][chr2][ps2];

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
		    } // end iterating over strains

		  // if we've seen all combinations of alleles, test
		  // for association
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
			  {
			    if(Ci(k) > 100)
			      {
				cout << "LARGE\n";
				cout << "coefficients: " <<  C(0) << ' ' << C(1) << ' ' << C(2) << ' ' << C(3) << '\n';
				cout << flush;
			      }
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

			// check whether this is the strongest
			// associaiton we've found so far
			if(LOD > max_LOD)
			  {
			    max_LOD = LOD;
			    max_LOD_ind_1 = pair<int, int>(chr1, ps1);
			    max_LOD_ind_2 = pair<int, int>(chr2, ps2);
			  }

			// store LOD score regardless
			vector<int> loc = vector<int>();
			loc.push_back(chr1);
			loc.push_back(ps1);
			loc.push_back(chr2);
			loc.push_back(ps2);
			LOD_scores_loc.push_back(loc);
			LOD_scores.push_back(LOD);
		      }
		      
		    } 
		} // end iterating over ps2
	    }
	  //	}
//    }

  gsl_vector_free(Y);

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(c);

  gsl_matrix_free(Xi);
  gsl_matrix_free(covi);
  gsl_vector_free(ci);

  max_LOD_inds.push_back(max_LOD_ind_1);
  max_LOD_inds.push_back(max_LOD_ind_2);

  h2_pred = 1 - pow(10, -2.0 / n * max_LOD);

  cout << "number sites with all allele combos: " << success_count << '\n' << flush;

  return max_LOD;
}

void write_output_header(ofstream& f_summary, int num_qtls)
{
  for (int i = 1; i <= num_qtls; i++)
    {
      f_summary << "pred_qtl_chr" << i << ' ' << "pred_qtl_ps" << i << ' ';
    }
  for (int i = 1; i <= num_qtls; i++)
    {
      f_summary << "qtl_chr" << i << ' ' << "qtl_ps" << i << ' ';
    }

  f_summary << "LOD_score h2_prediction ";

  for (int i = 1; i <= num_qtls; i++)
    {
      f_summary << "allele_frequency" << i << ' ';
    }

  f_summary << '\n';

  f_summary.flush();
}

void write_LOD_header(ofstream& f_LOD, int num_sims, 
		      vector<int>& chr, vector<vector<int> >& ps)
{
  for (int c = 0; c < chr.size(); c++)
    {
      for (int p = 0; p < ps[c].size(); p++)
	{
	  f_LOD << ',' << chr[c];
	}
    }
  f_LOD << '\n';

  for (int c = 0; c < chr.size(); c++)
    {
      for (int p = 0; p < ps[c].size(); p++)
	{
	  f_LOD << ',' << ps[c][p];
	}
    }
  f_LOD << '\n';
  f_LOD << flush;
}


void write_stats(ofstream& f_summary,
		 ofstream& f_LOD,
		 vector<int>& chr,
		 vector<vector<int> >& ps,
		 vector<pair<int, int> > qtl_pred, 
		 vector<pair<int, int> > qtl,
		 double max_LOD, double h2_pred,
		 vector<double> af, vector<vector<double> >& LOD_scores,
		 int sim)
{
  for (int i = 0; i < qtl_pred.size(); i++)
    {
      int chr_ind = qtl_pred[i].first;
      int ps_ind = qtl_pred[i].second;
      f_summary << chr[chr_ind] << ' ' << ps[chr_ind][ps_ind] << ' ';
    }

  for (int i = 0; i < qtl_pred.size(); i++)
    {
      int chr_ind = qtl[i].first;
      int ps_ind = qtl[i].second;
      f_summary << chr[chr_ind] << ' ' << ps[chr_ind][ps_ind] << ' ';
    }

  f_summary << max_LOD << ' ' << h2_pred << ' ';

  for (int i = 0; i < qtl_pred.size(); i++)
    {
      f_summary << af[i] << ' ';
    }

  f_summary << '\n';
  f_summary << flush;

  f_LOD << sim;
  for (int c = 0; c < chr.size(); c++)
    {
      for (int p = 0; p < ps[c].size(); p++)
	{
	  f_LOD << ',' << LOD_scores[c][p];
	}
    }
  f_LOD << '\n';
  f_LOD << flush;
}

void write_stats_2(ofstream& f_summary,
		 ofstream& f_LOD,
		 vector<int>& chr,
		 vector<vector<int> >& ps,
		 vector<pair<int, int> > qtl_pred, 
		 vector<pair<int, int> > qtl,
		 double max_LOD, double h2_pred,
		 vector<double> af, vector<double>& LOD_scores,
		   vector<vector<int> > LOD_scores_loc,
		 int sim)
{
  for (int i = 0; i < qtl_pred.size(); i++)
    {
      int chr_ind = qtl_pred[i].first;
      int ps_ind = qtl_pred[i].second;
      f_summary << chr[chr_ind] << ' ' << ps[chr_ind][ps_ind] << ' ';
    }

  for (int i = 0; i < qtl_pred.size(); i++)
    {
      int chr_ind = qtl[i].first;
      int ps_ind = qtl[i].second;
      f_summary << chr[chr_ind] << ' ' << ps[chr_ind][ps_ind] << ' ';
    }

  f_summary << max_LOD << ' ' << h2_pred << ' ';

  for (int i = 0; i < qtl_pred.size(); i++)
    {
      f_summary << af[i] << ' ';
    }

  f_summary << '\n';
  f_summary << flush;

  // chr1
  for (int i = 0; i < LOD_scores.size(); i++)
    {
      f_LOD << ',' << chr[LOD_scores_loc[i][0]];
    }
  f_LOD << '\n';
  // ps1
  for (int i = 0; i < LOD_scores.size(); i++)
    {
      f_LOD << ',' << ps[LOD_scores_loc[i][0]][LOD_scores_loc[i][1]];
    }
  f_LOD << '\n';
  // chr2
  for (int i = 0; i < LOD_scores.size(); i++)
    {
      f_LOD << ',' << chr[LOD_scores_loc[i][2]];
    }
  f_LOD << '\n';
  // ps2
  for (int i = 0; i < LOD_scores.size(); i++)
    {
      f_LOD << ',' << ps[LOD_scores_loc[i][2]][LOD_scores_loc[i][3]];
    }
  f_LOD << '\n';
  // LOD
  for (int i = 0; i < LOD_scores.size(); i++)
    {
      f_LOD << ',' << LOD_scores[i];
    }
  f_LOD << '\n';
  f_LOD << flush;
}

double heritability_to_fixed_effect(double h2, 
				    vector<double> af, int n)
{
  double avg_af = ((double)accumulate(af.begin(), af.end(), 0.0))/af.size();
  return pow(h2 * (n-1) / ((h2-1) * n * (avg_af-1) * avg_af), .5);
}

vector<vector<vector<int> > > get_parent_gen(string fn, vector<int>& chr, 
					     vector<vector<int> >& ps, 
					     vector<vector<double> >& af)
{
  // chr,ps,af,strain,strain,strain
  // chr,ps,af,allele,allele,allele
  // chr,ps,af,allele,allele,allele
  ifstream f;
  f.open(fn);
  string line;
  getline(f, line);
  int index = 0;
  int new_index = -1;
  int col_ind = 0;
  while (index < line.length())
    {
      new_index = line.find(',', index + 1);
      // could store strain names here if we cared about them...
      index = new_index;      
      col_ind++;
    }
  int num_parents = col_ind - 3;

  int prev_chr = -1;
  int current_chr = -1;
  int chr_ind = -1;
  vector<vector<vector<int> > > parents;
  for (int i = 0; i < num_parents; i++)
    {
      parents.push_back(vector<vector<int> >());
    }
  while(getline(f, line))
    {
      index = 0;
      new_index = line.find(',', index + 1);
      current_chr = atoi(line.substr(index, new_index - index).c_str());
      index = new_index;

      if (current_chr != prev_chr)
	{
	  for (int i = 0; i < num_parents; i++)
	    {
	      parents[i].push_back(vector<int>());
	    }
	  ps.push_back(vector<int>());
	  af.push_back(vector<double>());
	  chr.push_back(current_chr);
	  prev_chr = current_chr;
	  chr_ind++;
	}

      new_index = line.find(',', index + 1);
      ps[chr_ind].push_back(atoi(line.substr(index + 1, new_index - index - 1).c_str()));
      index = new_index;

      new_index = line.find(',', index + 1);
      af[chr_ind].push_back(atof(line.substr(index + 1, new_index - index - 1).c_str()));
      index = new_index;

      for (int i = 0; i < num_parents; i++)
	{
	  new_index = line.find(',', index + 1);
	  parents[i][chr_ind].push_back(atoi(line.substr(index + 1, new_index - index - 1).c_str()));
	  index = new_index;
	}
    }

  return parents;
}

void sim(int num_sims, int num_segregants, int num_per_cross, string outfile_ext, 
	 int num_qtls, double h2, 
	 string parent_file, vector<int> first_pos, vector<int> last_pos)
{
  // read in parent generation genotypes and information about the
  // variant sites
  vector<int> chr;
  vector<vector<int> > ps;
  vector<vector<double> > af;
  // minor allele is always 0, major allele 1
  vector<vector<vector<int> > > parents = get_parent_gen(parent_file, chr, ps, af);
  int num_snps = 0;
  for (int i = 0; i < chr.size(); i++)
    {
      num_snps += ps[i].size();
    }

  // open output files and write column headers
  string summary_fn = "out/sim_output_summary" + outfile_ext + ".txt";
  ofstream f_summary;
  f_summary.open(summary_fn);
  write_output_header(f_summary, num_qtls);

  string LOD_fn = "out/sim_output_LOD" + outfile_ext + ".txt";
  ofstream f_LOD;
  f_LOD.open(LOD_fn);
  if (num_qtls == 1)
    write_LOD_header(f_LOD, num_sims, chr, ps);

  // initialize random seed (note that time resolution (1 second?) can
  // be an issue here if we're starting a bunch of simulations at
  // close to the same time)
  srand(time(NULL));

  // random number generator and distributions for the two alleles
  // (one is centered at 0 and one at a fixed effect to be determined
  // later based on the allele frequencies of the snvs we choose)
  default_random_engine generator;
  normal_distribution<double> n_dist_0(0, 1);

  // vectors for storing segregants, their phenotypes, and the qtls
  // and their fixed effects for each simulation
  vector<vector<vector<int> > > segregants;
  vector<double> phenotypes;
  vector<pair<int, int> > qtls;
  double fixed_effect;

  // run simulations
  for(int s = 0; s < num_sims; s++)
    {
      cout << "simulation " << s << '\n';
      cout.flush();

      // generate segregants by simulating mating of parents in a
      // funnel cross
      segregants = gen_segregants(parents, generator, chr, ps, 
				  first_pos, last_pos, num_segregants, 
				  num_per_cross);
      cout << "generated segregants\n" << flush;

      // choose qtls (all on different chromosomes)
      int num_indices_to_choose = num_snps;
      vector<int> chosen_chrs;
      //vector<int>::iterator start_it = chosen_chrs.begin();
      //vector<int>::iterator end_it = chosen_chrs.end();
      
      qtls.clear();
      for(int q = 0; q < num_qtls; q++)
	{
	  int index = (int)(round(rand() % num_indices_to_choose));
	  int chosen_chr = 0;
	  int chosen_ps = -1;
	  while (true)
	    {
	      if (find(chosen_chrs.begin(), chosen_chrs.end(), chosen_chr) !=
		  chosen_chrs.end())
		{
		  chosen_chr++;
		}
	      else if (index >= ps[chosen_chr].size())
		{
		  index -= ps[chosen_chr].size();
		  chosen_chr++;
		}
	      else
		{
		  chosen_ps = index;
		  break;
		}
	    }
	  num_indices_to_choose -= ps[chosen_chr].size();
	  chosen_chrs.push_back(chosen_chr);
	  qtls.push_back(pair<int, int>(chosen_chr, chosen_ps));
	}

      vector<double> qtl_af;
      for(int i = 0; i < num_qtls; i++)
	{
	  qtl_af.push_back(af[qtls[i].first][qtls[i].second]);
	}
      fixed_effect = heritability_to_fixed_effect(h2, qtl_af, num_segregants);

      phenotypes = sim_phenotypes(segregants, qtls, fixed_effect, generator);

      vector<pair<int, int> > max_LOD_inds;
      double h2_pred;
      vector<vector<double> > LOD_scores;
      double max_LOD;
      if (num_qtls == 1)
	{
	  pair<int, int> max_LOD_ind;
	  max_LOD = predict_one(phenotypes, segregants, ps,
				max_LOD_ind,
				h2_pred, LOD_scores);
	  max_LOD_inds.push_back(max_LOD_ind);
	  // store all the results
	  write_stats(f_summary, f_LOD, chr, ps, max_LOD_inds, qtls,  max_LOD, h2_pred,
		      qtl_af, LOD_scores, s);
	}
      else if (num_qtls == 2)
	{
	  vector<double> LOD_scores_2;
	  vector<vector<int> > LOD_scores_loc;
	  vector<pair<int, int> > ranges;
	  int window = 300;
	  int qtl1_range_start = max(0, qtls[0].second - window);
	  int qtl1_range_end = min(ps[qtls[0].first].size(), (size_t)(qtls[0].second + window));
	  int qtl2_range_start = max(0, qtls[1].second - window);
	  int qtl2_range_end = min(ps[qtls[1].first].size(), (size_t)(qtls[1].second + window));
	  ranges.push_back(pair<int, int>(qtl1_range_start, qtl1_range_end));
	  ranges.push_back(pair<int, int>(qtl2_range_start, qtl2_range_end));
	  vector<int> correct_chrs;
	  correct_chrs.push_back(qtls[0].first);
	  correct_chrs.push_back(qtls[1].first);
	  max_LOD = predict_two(phenotypes, segregants, ps,
				max_LOD_inds,
				h2_pred, LOD_scores_2, LOD_scores_loc,
				ranges, correct_chrs);
	  write_stats_2(f_summary, f_LOD, chr, ps, max_LOD_inds, qtls,  
			max_LOD, h2_pred,
			qtl_af, LOD_scores_2, LOD_scores_loc, s);
	}
      else
	{
	  cout << "don't know how to deal with " << num_qtls << " qtls\n" << flush;
	}

    }
  f_summary.close();
  f_LOD.close();
}

