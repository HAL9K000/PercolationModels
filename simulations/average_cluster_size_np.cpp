#include "cluster_dynamics.h"

int main(){

	increase_stack_limit(64);

	int grid_size;
	float p_start;
	float p_end;
	int divisions;
	int r_init;
	int number_of_census;
	int lag;

	cout << "Enter grid size: ";
    cin >> grid_size;

    cout << "Enter starting p (should be close to p_c = 0.592746): ";
    cin >> p_start;

  	cout << "Enter ending p: ";
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

		stringstream p_st, p_en;

	  p_st << setprecision(3) << p_start;
	  p_en << setprecision(3) << p_end;
	  // setprecision() is a stream manipulator that sets the decimal precision of a variable.

		ofstream outputgamma;
	  // Creating a file instance called outputsigma to store output data as CSV.

		outputgamma.open("CrtExp/Gamma_SP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_R_" + std::to_string(r_init) + "_Cen_"+ std::to_string(number_of_census) + ".csv");
	  // Creating CSV file in "CrtExp" sub-directory to store output data

    vector<double> birth_probabilities_np = linspace(p_start, p_end, divisions);

    vector<zd_coordinates> avg_cluster_size_np;

    auto start = high_resolution_clock::now();

	#pragma omp parallel
	{
	    vector<zd_coordinates> avg_cluster_sizes_np_private;

	    #pragma omp for nowait schedule(static)
			for (int i=0; i < divisions; i++)
			{

			int seed = std::random_device{}();
			rng.seed(seed);

			float birth_probability_np = birth_probabilities_np[i];
			std::vector<zd_coordinates> gamma_data;

		  average_cluster_size_np(grid_size, gamma_data, birth_probability_np, r_init, number_of_census,lag);
			avg_cluster_sizes_np_private.insert(avg_cluster_sizes_np_private.end(), gamma_data.begin(), gamma_data.end());

		 }

		#pragma omp for schedule(static) ordered
		for(int i=0; i< omp_get_num_threads(); i++)
		{
			#pragma omp ordered
				avg_cluster_size_np.insert(avg_cluster_size_np.end(), avg_cluster_sizes_np_private.begin(), avg_cluster_sizes_np_private.end());
				// Inserting average cluster size data for each trial in order.
		}
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	/*
	double p_c = 0.592746;
	//Percolation point for 2D Grid Network.
	std::vector<double> lnbase; 	// Stores ln(p - pc) values for various values of p.

	// Storing ln| p - p_c | for p ----> p_c
	for(int i=0; i<divisions; i++)
	{
		double x =0.0;
		if (birth_probabilities_np[i] < p_c)
		{
			x= log(p_c - birth_probabilities_np[i]);
		}
		else
		{
			x= log(birth_probabilities_np[i] - p_c);
		}

		lnbase.push_back(x);
	} */

	cout << "The vector elements are: "<< endl;
  cout << "# Trial Number, \t p, \t S(p)\n";
	outputgamma << "# Trial Number, p, S(p)\n";

	// Printing out and writing home obtained results to CSV file.

	for (int i=0; i< avg_cluster_size_np.size(); i++)
	{
		cout << setprecision(6) << avg_cluster_size_np[i].x << "  " << setprecision(7) << avg_cluster_size_np[i].y << "  " << setprecision(8) << avg_cluster_size_np[i].z << endl;
		outputgamma << setprecision(7) << avg_cluster_size_np[i].x << "," << setprecision(8) << avg_cluster_size_np[i].y << "," << setprecision(10) << avg_cluster_size_np[i].z << endl;
	}
	outputgamma.close();



	cout << endl << "CPU Time: " << duration.count() << " Seconds" << endl;

}
