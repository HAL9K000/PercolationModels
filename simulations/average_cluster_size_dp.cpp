#include "cluster_dynamics.h"

int main(){

	//OpenMP array reduce idiom.
	//Source: https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
	//Credit: https://stackoverflow.com/users/2542702/z-boson
	//Used with modification 
	
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

    vector<double> birth_probabilities_dp = linspace(p_start, p_end, divisions);

    float avg_cluster_size_dp[divisions] = {0};

    auto start = high_resolution_clock::now(); 

	#pragma omp parallel
	{
	    float avg_cluster_sizes_dp_private[divisions] = {0};

	    #pragma omp for
		for (int i=0; i < divisions; i++){

			int seed = std::random_device{}();
			rng.seed(seed); 
			
			float birth_probability_dp = birth_probabilities_dp[i];

			avg_cluster_sizes_dp_private[i] = average_cluster_size_dp(grid_size, birth_probability_dp, number_of_census,lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            avg_cluster_size_dp[n] += avg_cluster_sizes_dp_private[n];
	        }
	    }
	}

	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<seconds>(stop - start); 

	for (int i=0; i< divisions; i++){
		cout << "p: " << setprecision(3) << birth_probabilities_dp[i] << " S: " << setprecision(3) << avg_cluster_size_dp[i] << endl;
	}

	cout << endl << "CPU Time: " << duration.count() << " Seconds" << endl;

}