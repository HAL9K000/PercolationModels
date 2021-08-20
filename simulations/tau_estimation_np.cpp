#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(512);

	int grid_size;
	float p_start;
  float p_end;
	float density;
	//int update_per_site; Going to default to 10,000 for NP.
	int number_of_census;
	int divisions;
  int r_init;
	int lag;
	int collect_frames = 1;
	vector <int> cluster_sizes;

  cout << "The point of this module is to generate ns(p) distributions to estimate tau for 1 << s << s_e for p ---> pc-\n";

	cout << "Enter grid size: ";
  cin >> grid_size;

  cout << "Enter starting p (should be close to p_c = 0.592746): ";
  cin >> p_start;

  cout << "Enter ending p: ( p ----> p_c-): ";
  cin >> p_end;

  cout << "Enter divisions: ";
  cin >> divisions;

  cout << "Enter number of random censuses to be initiated (suggested choice b/w 3-5): ";
  cin >> r_init;

  cout << "Enter number of census: ";
  cin >> number_of_census;

  cout << "Enter lag in terms of frames: ";
  cin >> lag;

  cout << endl;

  auto start = high_resolution_clock::now();

  iter_patch_sizes_np(grid_size, p_start, p_end, divisions, r_init, number_of_census, lag);
  /* Will store ns(p) distributions for all specified p in a CSV file. Booyaa */

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);

  cout << endl << "CPU Time: " << duration.count() << " seconds" << endl;
}
