#include "cluster_dynamics.h"

int main() {

	increase_stack_limit(64);

	int grid_size;
	float birth_probability;
	float density;
	int updates_per_site;
	int how_many;
	int replicate;
	vector <transformation> transformations;

	cout << "Enter grid size: ";
    cin >> grid_size;

    cout << "Enter birth probability: ";
	cin >> birth_probability;

	cout << "What is the corresponding density? : ";
	cin >> density;


    cout << "Enter updates per site: ";
    cin >> updates_per_site;

    cout << "How many transformations do you want? Enter a number: ";
    cin >> how_many;

    cout << "Which replicate is this? : ";
    cin >> replicate;



	find_equilibrium_single_shot_transformations_dp(grid_size, birth_probability, updates_per_site, how_many, transformations);

	ofstream outdata;

	stringstream d;
	d << setprecision(3) << density;


	outdata.open("dump/dp_transformations_"+std::to_string(grid_size)+"_"+d.str()+"_"+std::to_string(how_many)+"_"+std::to_string(replicate)+".txt");
	if( !outdata ) {
		cerr << "File error, try again." << endl;
		exit(1);
	}
	for (int i=0; i< transformations.size(); i++){
		outdata << transformations[i].before << " " << transformations[i].after << endl;
	}
	outdata.close();
}