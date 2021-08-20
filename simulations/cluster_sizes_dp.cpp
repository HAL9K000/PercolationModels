#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(64);

	int grid_size;
	float birth_probability;
	float density;
	int update_per_site;
	int number_of_census;
	int replicate;
	int lag;
	int collect_frames = 1;
	vector <int> cluster_sizes;

	cout << "Enter grid size: ";
    cin >> grid_size;

    cout << "Enter birth probability: ";
	cin >> birth_probability;

	cout << "What is the corresponding density? : ";
	cin >> density;


    cout << "Enter updates per site: ";
    cin >> update_per_site;

    cout << "Enter number of census: ";
    cin >> number_of_census;

    cout << "Enter lag between frames: ";
    cin >> lag;

    cout << "Which replicate is this? : ";
    cin >> replicate;

    cout << "Do you want to collect frames? Answer with 1 or 0: ";
    cin >> collect_frames;


	find_patch_sizes_dp(cluster_sizes, grid_size, birth_probability, number_of_census, lag, update_per_site, collect_frames);

	ofstream outdata;

	stringstream d;
	d << setprecision(3) << density;


	outdata.open("dump/dp_cluster_sizes_"+std::to_string(grid_size)+"_"+d.str()+"_"+std::to_string(number_of_census)+"_"+std::to_string(replicate)+".txt");
	if( !outdata ) {
		cerr << "File error, try again." << endl;
		exit(1);
	}
	for (int i=0; i< cluster_sizes.size(); i++){
		outdata << cluster_sizes[i] << " ";
	}
	outdata.close();
}