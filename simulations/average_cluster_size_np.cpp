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

    vector<double> birth_probabilities_np = linspace(p_start, p_end, divisions);

    float avg_cluster_size_np[divisions] = {0};

    auto start = high_resolution_clock::now(); 

	#pragma omp parallel
	{
	    float avg_cluster_sizes_np_private[divisions] = {0};

	    #pragma omp for
		for (int i=0; i < divisions; i++){

			int seed = std::random_device{}();
			rng.seed(seed); 
			
			float birth_probability_np = birth_probabilities_np[i];

			avg_cluster_sizes_np_private[i] = average_cluster_size_np(grid_size, birth_probability_np, number_of_census,lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            avg_cluster_size_np[n] += avg_cluster_sizes_np_private[n];
	        }
	    }
	}

	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<seconds>(stop - start); 

	for (int i=0; i< divisions; i++){
		cout << "p: " << setprecision(3) << birth_probabilities_np[i] << " S: " << setprecision(3) << avg_cluster_size_np[i] << endl;
	}

	cout << endl << "CPU Time: " << duration.count() << " Seconds" << endl;

}