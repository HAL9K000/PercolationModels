#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(64);

	int grid_size;
	float p_start;
	float p_end;
	int divisions;
	int r_init;
	int length_of_census;

	cout << "Enter grid size: ";
    cin >> grid_size;

    cout << "Enter starting p: ";
    cin >> p_start;

  	cout << "Enter ending p: ";
    cin >> p_end;

    cout << "Enter divisions: ";
    cin >> divisions;

		cout << "Enter number of random trials for each value of p: ";
    cin >> r_init;

    cout << "Enter trial length for a given p: ";
    cin >> length_of_census;

    cout << endl;

	auto start = high_resolution_clock::now();

	crtexp_dynamo_dp(grid_size,p_start,p_end,divisions,r_init, length_of_census);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "CPU Time: " << duration.count() << " seconds" << endl;
}
