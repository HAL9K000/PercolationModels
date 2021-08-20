#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(64);

	int grid_size;
	float p_start;
	float p_end;
	int divisions;
	int number_of_census;
	int lag;

	cout << "Enter grid size: ";
    cin >> grid_size;

    cout << "Enter starting p: ";
    cin >> p_start;

  	cout << "Enter ending p: ";
    cin >> p_end;

    cout << "Enter divisions: ";
    cin >> divisions;

    cout << "Enter number of census: ";
    cin >> number_of_census;

    cout << "Enter lag in terms of frames: ";
    cin >> lag;

    cout << endl;

	auto start = high_resolution_clock::now(); 

	percolation_probabilities_dp(grid_size,p_start,p_end,divisions,number_of_census,lag);

	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<seconds>(stop - start); 

	cout << endl << "CPU Time: " << duration.count() << " seconds" << endl;
}