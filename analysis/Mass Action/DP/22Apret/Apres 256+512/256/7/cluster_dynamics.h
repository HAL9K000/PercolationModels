#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <cstdlib>
#include <random>
#include <vector>
#include <algorithm>
static std::mt19937_64 rng(time(NULL));
#include <chrono>
using namespace std::chrono;
#include<bits/stdc++.h>
using namespace std;
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
#include <fstream>
#include <sstream>

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
void zeros(int frame[], int grid_size);
float calculate_density(int frame[], int grid_size);
f_coordinates calculate_SD(int frame[],int grid_size);
void increase_stack_limit(int stack_size);
int disjoint_vectors(vector<int> vec_1, vector<int> vec_2);
vector<int> common_elements(vector<int> vec_1, vector<int> vec_2);
template<typename T>std::vector<double> linspace(T start_in, T end_in, int num_in);
void print_vector(std::vector<double> vec);

// -------------------------------Percolation Models-------------------------------------------//


coordinates select_neighbor_of_site(coordinates site, int grid_size);
coordinates select_neighbor_of_pair(coordinates site, coordinates neighbour, int grid_size);
void np_update(int frame[], int grid_size, float birth_probability);
void dp_update(int frame[], int grid_size, float birth_probability);
void tp_update(int frame[], int grid_size, float birth_probability, float feedback_strength);
void simulate_np(int frame[], int grid_size, float birth_probability, int updates_per_site);
void simulate_dp(int frame[], int grid_size, float birth_probability, int updates_per_site);
void simulate_tp(int frame[], int grid_size, float birth_probability, float feedback_strength, int updates_per_site);
float equilibrium_density_dp(int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);
float equilibrium_density_tp(int grid_size, float birth_probability, float feedback_strength, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);

// -----------------------------------Cluster Statics -----------------------------------------//

vector<int> spanning_cluster_label(int frame[], int grid_size);
void find_clusters_without_spanning(int frame[], int labels[],vector<cluster>& clusters, int grid_size);
void find_patch_sizes_dp(vector<int>& cluster_sizes,int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);
void find_patch_sizes_np(vector<int>& cluster_sizes,int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site=100000, int collect_frames=0);

//------------------------------------Cluster Dynamics----------------------------------------//

void find_neighbours_of_site(coordinates neighbours[4], coordinates site, int grid_size);
void depth_first_search(coordinates site, vector<int>& sites, int frame[], int labels[], int id, int grid_size);
void find_clusters(int frame[], int labels[],vector<cluster>& clusters, int grid_size);
int check_presence(std::vector<int> vec,int element);
void find_transformations_single_shot(vector<transformation>& transformations, int previous_frame[], int previous_frame_labels[],vector<cluster>& previous_frame_clusters,  int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size);
void find_transformations_multi_shot(vector<transformation>& transformations, int previous_frame[], int previous_frame_labels[],vector<cluster>& previous_frame_clusters, int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size);
void find_equilibrium_single_shot_transformations_np(int grid_size, float birth_probability, int time_to_equilibrium, int how_many, vector<transformation>& transformations);
void find_equilibrium_single_shot_transformations_dp(int grid_size, float birth_probability, int time_to_equilibrium, int how_many, vector<transformation>& transformations);

void find_equilibrium_multi_shot_transformations_DP(int grid_size, int divisions, int time_to_equilibrium, int how_many, vector<transformation>& transformations, vector<f_coordinates>& pL);
void DP_scattershot(vector <transformation>& trans, int previous_frame[], int previous_frame_labels[], vector<cluster>& previous_frame_clusters, int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size, int how_many, double p, double lag);
void find_equilibrium_single_shot_transformations_tp(int grid_size, float birth_probability, float feedback_strength, int time_to_equilibrium, int how_many, vector<transformation>& transformations);
void TCP_scattershot(vector <transformation>& trans, int previous_frame[], int previous_frame_labels[], vector<cluster>& previous_frame_clusters,
int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size, int how_many, double p, double q);

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
float average_cluster_size_np(int grid_size, float birth_probability, int number_of_census, int lag);
float average_cluster_size_dp(int grid_size, float birth_probability, int number_of_census, int lag);

//---------------------------Correlation Function-----------------------------------------------//

zd_coordinates crosscorrelation_2D(int frame1[], int frame2[], int grid_size);
void CrossCol_DP(int grid_size, float p, int number_of_census, int r_init, int lag);

#endif
