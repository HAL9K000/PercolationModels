#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(1024);

	int grid_size;
	float birth_probability;
	int r_init;
	int updates_per_site;
	int lag;
	int number_of_census;

	int t_to_eq;

	/** cout << "Enter grid size: ";
  cin >> grid_size;
  cout << "Enter p value of contention: ";
  cin >> birth_probability;
  cout << "Enter number of census: ";
  cin >> number_of_census;
  cout << "Enter lag: ";
  cin >> lag;
	cout << "Enter number of random trials: ";
  cin >> r_init; */

	number_of_census = 200; grid_size =256; birth_probability =0.65; r_init =8; lag=1;

  auto start = high_resolution_clock::now();
  // Auto variable deduces type of variable by itself.


  CrossCol_DP(grid_size, birth_probability, number_of_census, r_init, lag);

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);

  cout << endl << "Total CPU Time: " << duration.count() << " seconds" << endl;
  }
