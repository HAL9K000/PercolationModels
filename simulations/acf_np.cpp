#include "cluster_dynamics.h"

int main()
{

	increase_stack_limit(64);

  int grid_size;
	float p;
	int divisions;     // Number of random trials to be initiated and averaged over to calculate ACF.
	int length;      //Extent of time period over which ACF will be calculated.
	int lag=1;       // Time intervals at which ACF will be sampled.


  cout << "Enter grid size: ";
  cin >> grid_size;

  cout << "Enter p: ";
  cin >> p;

  cout << "Enter number of trials (frame initialisations) over which ACF will be averaged: ";
  cin >> divisions;

  cout << "Enter time duration over which ACF will be calculated (+ve integral value only)";
  cin >> length;

  cout << "Enter discrete time intervals at which ACF will be sampled (default is 1) ";
  cin >> lag;

  cout << endl;

  auto start = high_resolution_clock::now();
  // Auto variable deduces type of variable by itself.

  acf_np_custom(grid_size, p, divisions, length, lag);
  // Outputs ACF data.

  auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "CPU Time: " << duration.count() << " seconds" << endl;

  return 0;

}
