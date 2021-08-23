#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(1024);

	int grid_size;
	float birth_probability;
	float density;
	int updates_per_site;
	int how_many;
	int replicate;

	int t_to_eq;

	/* cout << "Enter grid size: ";
  cin >> grid_size;

  cout << "Enter time taken to arrive at equilibrium:\t ";
  cin >> t_to_eq;

  cout << "How many transformations do you want? Enter a number: ";
  cin >> how_many; */

	grid_size=128; t_to_eq =100000; how_many=50000;
	int divisions=23;


	f_coordinates PQ[divisions] = {{0.5, 0.04}, {0.52, 0.04}, {0.54, 0.04}, {0.56, 0.04}, {0.58, 0.04}, {0.6, 0.04}, {0.62, 0.04}, {0.64, 0.04},
	{0.66, 0.04}, {0.68, 0.04}, {0.70, 0.04}, {0.72, 0.04}, {0.722279867, 0.04}, {0.74, 0.04}, {0.76, 0.04}, {0.78, 0.04},
	{0.8, 0.04}, {0.82, 0.04}, {0.84, 0.04}, {0.86, 0.04}, {0.88, 0.04}, {0.9, 0.04}, {0.92, 0.04}};

	std::vector<f_coordinates> pq; pq.assign(PQ, PQ + divisions);
	std::vector<transformation> transformations;

	cout << "p\t|  q\n";
	for (int i=0; i < divisions; i++)
	{
			cout << pq[i].x << "  " << pq[i].y <<endl;
	}

	auto start = high_resolution_clock::now();
  // Auto variable deduces type of variable by itself.


	find_equilibrium_single_shot_transformations_tp(grid_size, divisions, t_to_eq, how_many, transformations, pq);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Total CPU Time: " << duration.count() << " seconds" << endl;
}
