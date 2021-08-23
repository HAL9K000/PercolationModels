#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
using namespace std::chrono;
static std::mt19937_64 rng(time(NULL));
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <omp.h>

#ifndef CLUSTER_DYNAMICS_MC_H
#define CLUSTER_DYNAMICS_MC_H

//---------------------Some Primitives ------------------------------------------------------ //

struct coordinates {
   int x;
   int y;
};

struct cluster {
	int label;
	std::vector<int> coords;
};

struct transformation {
	int before;
	int after;
};

struct f_coordinates {
   float x;
   float y;
};

struct zd_coordinates {
   double x;
   double y;
   double z;
};

double random_real(double initial, double last);
int random_int(int initial, int last);
float mean_of_array(float array[],int size);
float standard_deviation_of_array(float array[],int size);
float mean_of_vector(vector<float> array,int size);
void random_frame(int frame[], int grid_size);
void random_frame_of_density(float density, int frame[], int grid_size);
void solitary_droplet(int frame[], int grid_size);
void zeros(int frame[], int grid_size);
float calculate_density(int frame[], int grid_size);
int number_of_elem_of_array(int frame[],int grid_size);
void increase_stack_limit(int stack_size);
int disjoint_vectors(vector<int> vec_1, vector<int> vec_2);
vector<int> common_elements(vector<int> vec_1, vector<int> vec_2);
template<typename T>std::vector<double> linspace(T start_in, T end_in, int num_in);
void print_vector(std::vector<double> vec);

void output_dump(int frame[], int grid_size, int i, double p, int labels[], vector <cluster>& clusters, vector<int>& spanning_cluster_labels);

// -------------------------------Percolation Models-------------------------------------------//


coordinates select_neighbor_of_site(coordinates site, int grid_size);
coordinates select_neighbor_of_pair(coordinates site, coordinates neighbour, int grid_size);
void np_update(int frame[], int grid_size, float birth_probability);
void dp_update(int frame[], int grid_size, float birth_probability);
void tp_update(int frame[], int grid_size, float birth_probability, float feedback_strength);
void simulate_np(int frame[], int grid_size, float birth_probability, int updates_per_site);
void simulate_dp(int frame[], int grid_size, float birth_probability, int updates_per_site);
void simulate_tp(int frame[], int grid_size, float birth_probability, float feedback_strength, int updates_per_site);
float equilibrium_density_dp(int grid_size, float birth_probability, int r_init=1, int number_of_census=5, int lag=2000, int updates_per_site=100000, int collect_frames=0);
float equilibrium_density_tp(int grid_size, float birth_probability, float feedback_strength, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);

// -----------------------------------Cluster Statics -----------------------------------------//

vector<int> spanning_cluster_label(int frame[], int grid_size);
void find_clusters_without_spanning(int frame[], int labels[],vector<cluster>& clusters, int grid_size);
void find_patch_sizes_dp(vector<int>& cluster_sizes,int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);
void find_patch_sizes_np(vector<int>& cluster_sizes,int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);

// ---------------------------------Critical Exponent Tau ------------------------------------//
bool sort_s_cluster(coordinates a, coordinates b);
void tau_patch_size_find_np(int grid_size,vector<zd_coordinates>& tau_data, double p, int r_init=5, int number_of_census=25, int lag=12);
void iter_patch_sizes_np(int grid_size, float p_start, float p_end, int divisions, int r_init=5, int number_of_census=25, int lag=12);

//-----------------------------------Critical Exponent Tau (TCP)--------------------------------//

void tau_patch_size_find_tcp(int grid_size, vector<zd_coordinates>& tau_data, double p, double q, int r_init, int number_of_census, int lag);

//------------------------------------Cluster Dynamics----------------------------------------//

void find_neighbours_of_site(coordinates neighbours[4], coordinates site, int grid_size);
void depth_first_search(coordinates site, vector<int>& sites, int frame[], int labels[], int id, int grid_size);
void find_clusters(int frame[], int labels[],vector<cluster>& clusters, int grid_size);
int check_presence(std::vector<int> vec,int element);
void find_transformations_single_shot(vector<transformation>& transformations, int previous_frame[], int previous_frame_labels[],vector<cluster>& previous_frame_clusters,  int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size);
void find_equilibrium_single_shot_transformations_np(int grid_size, float birth_probability, int time_to_equilibrium, int how_many, vector<transformation>& transformations);
void find_equilibrium_single_shot_transformations_dp(int grid_size, float birth_probability, int time_to_equilibrium, int how_many, vector<transformation>& transformations);
void find_equilibrium_single_shot_transformations_tp(int grid_size, int divisions, int time_to_equilibrium, int how_many, vector<transformation>& transformations, vector<f_coordinates>& pq);
void TCP_scattershot(vector <transformation>& trans, int previous_frame[], int previous_frame_labels[], vector<cluster>& previous_frame_clusters, int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size, int how_many, double p, double q);

// -----------------------------------This and That ------------------------------------------------//

void parallelizer(int processes, int grid_size, float birth_probabilities[], float feedback_strength, int number_of_census, int lag, int updates_per_site, int collect_frames);
void find_neighbours_of_site_free_boundary(coordinates neighbours[4], coordinates site, int grid_size);
void depth_first_search_free_boundary(coordinates site, vector<int>& sites, int frame[], int labels[], int id, int grid_size);
void find_clusters_free_boundary(int frame[], int labels[],vector<cluster>& clusters, int grid_size);
int is_spanning_vertical(int frame[], int grid_size);
int is_spanning_horizontal(int frame[], int grid_size);
int is_spanning(int frame[], int grid_size);
vector<int> spanning_cluster_coordinates(int frame[], int grid_size);
vector<int> spanning_cluster_label_id(int frame[], int grid_size, int labels[], vector<cluster>& clusters);
float how_many_red_sites(int frame[], int grid_size);
float find_average_cluster_size(int frame[], int grid_size);
float percolation_probability_np(int grid_size, float birth_probability, int number_of_census, int lag);
float percolation_probability_dp(int grid_size, float birth_probability, int number_of_census, int lag);
void percolation_probabilities_np(int grid_size, float p_start, float p_end, int divisions, int number_of_census, int lag);
void percolation_probabilities_dp(int grid_size, float p_start, float p_end, int divisions, int number_of_census, int lag);
void average_cluster_size_np(int grid_size,vector<zd_coordinates> &avg_clust_size, float birth_probability,int r_init, int number_of_census, int lag);
float average_cluster_size_dp(int grid_size, float birth_probability, int number_of_census, int lag);

// -----------------------------------Theoretical Percolation (NP) ------------------------------------------------//

int size_spanning_vertical(int frame[], int grid_size);
int size_spanning_horizontal(int frame[], int grid_size);
int size_spanning_2D(int frame[], int grid_size);
float theoretical_percolation_probability_np(int grid_size, float birth_probability, int r_init, int number_of_census, int lag);


// --------------------------------Theoretical Percolation Graphs & P_C estimation (DP) ---------------------------------??

void theoretical_percolation_probabilities_np(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag);
float theoretical_percolation_probability_dp(int grid_size, float birth_probability, int r_init, int number_of_census, int lag);
void theoretical_percolation_probabilities_dp(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag);
double theoret_percol_prob_denovo_dp(int grid_size, float birth_probability, int r_init, int number_of_census, int lag);
void calculate_pc_dp(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag);
void pavg_map_pc_dp(int grid_size, int r_init, int number_of_census, int lag);
//void pavg_map_pc_dp(int grid_size,float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag);

//----------------------------- P_C Estimation & Theoretical Percol (TCP)--------------------------------------------------//
void pavg_map_pc_tcp(double q, int grid_size, int r_init, int number_of_census, int lag);

// ---------------------------------- ACF NP ---------------------------------------------------------------------//

float xdynamic(int frame[], int grid_size);
void acf_np(int grid_size, vector<f_coordinates>& acf_data, float p, int length, int lag);
void acf_np_custom(int grid_size, float p, int divisions, int length, int lag);

// --------------------------------- Crtical Exponents (Beta, Gamma etc) [Finite Scaling] [TCP]------------------------//

zd_coordinates binsearch_p_c_TCP(double p, double q, int frame[], int grid_size, int num, int seed);
void crtexp_gamma_TCP(int grid_size,vector<zd_coordinates> &comp_data, double p, double q, int r_init, int number_of_census, int lag=5000);
void finite_scaling_crtexp_TCP(int grid_sizes[], double p, double q, string type, int divisions, int r_init, int number_of_census, int lag=5000);

// --------------------------------- Crtical Exponents (Beta, Gamma etc) [Finite Scaling] [NP/DP]------------------------//

zd_coordinates binsearch_p_c(double p, int frame[], int grid_size, int num, int seed);
void crtexp_nu(int grid_size,vector<zd_coordinates> &comp_data, int r_init, int number_of_census);
void crtexp_gamma(int grid_size,vector<zd_coordinates> &comp_data, double p, int r_init, int number_of_census);
void crtexp_beta_gamma(int grid_size,vector<zd_coordinates> &comp_data, double p, int r_init, int number_of_census);
void finite_scaling_crtexp(int grid_sizes[], double p, string type, int divisions, int r_init, int number_of_census);


//-------------------------------- Critical Exponent Calculation [DP Model]---------------------------------------//

void crtexp_DP_Basic(int grid_size,vector<zd_coordinates> &comp_data, double p, int r_init, int length);
void crtexp_dynamo_dp(int grid_size, double p_start, double p_end, int divisions, int r_init, int length);

#endif
