#include "cluster_dynamics.h"

int main()
{

	increase_stack_limit(256);

  float g1; float g2;
	double p;
	int divisions;     // Number of random trials to be initiated and averaged over to calculate ACF.
	int r_init=1;       // Number of random trials to be initiated per value of grid size.
  int number_of_census=1; string s = "";

  /* cout << "Enter p (must be very, very close to p_c): ";
  cin >> p;

  cout << "Enter starting grid size (as a power of 2, such as 5 ~ 32): ";
  cin >> g1;

  cout << "Enter ending grid size (as a power of 2, such as 9 ~ 512): ";
  cin >> g2;

  cout << "Enter number of points at which to sample Grid Length: ";
  cin >> divisions;

  cout << "Enter number of random trials to be inititated per value of the Grid Length: ";
  cin >> r_init;

  cout << "Enter number of censuses per random trial (given this is SP, default is 1): ";
  cin >> number_of_census;

	cout << "Enter type of scaling experiment to be performed (choose b/w 'Beta', 'Nu' or 'Gam' (default is Beta)): ";
  cin >> s; */

  cout << endl;

	p= 0.728; // Percolation Threshold For DP Universality Class
	r_init = 18; number_of_census = 12;
	divisions= 12; g1= 5; g2= 6.375;

  vector<double> grid_pow = linspace(g1, g2, divisions);

  int grid_sizes[divisions]; //Stores Grid Sizes at which P[p] and S[p] will be evaluated

  for( int i=0; i< grid_pow.size(); i++)
  {
    grid_sizes[i] = (int) pow(2, grid_pow[i]);
  }

  auto start = high_resolution_clock::now();
  // Auto variable deduces type of variable by itself.

  finite_scaling_crtexp(grid_sizes, p, s, divisions, r_init, number_of_census);

  auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "CPU Time: " << duration.count() << " seconds" << endl;

  return 0;

}
