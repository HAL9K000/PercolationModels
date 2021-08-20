#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(64);

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

	grid_size=128; t_to_eq =100000; how_many=10000;
	int divisions=8;


	f_coordinates PQ[divisions] = {{0.50, 0.25}, {0.55, 0.25}, {0.60, 0.25}, {0.65, 0.25}, {0.7, 0.25}, {0.75, 0.25}, {0.45, 0.25}, {0.40, 0.25},
	{0.70, 0.02}, {0.70, 0.03}, {0.70, 0.04}, {0.70, 0.05}, {0.7, 0.06}, {0.7, 0.07}, {0.7, 0.08}, {0.7, 0.09},
	{0.74, 0.02}, {0.74, 0.03}, {0.74, 0.04}, {0.74, 0.05}, {0.74, 0.06}, {0.74, 0.07}, {0.74, 0.08}, {0.74, 0.09}};

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
