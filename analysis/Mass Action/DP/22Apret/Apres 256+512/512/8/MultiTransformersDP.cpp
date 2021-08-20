#include "cluster_dynamics.h"

int main() {

	auto start = high_resolution_clock::now();
  // Auto variable deduces type of variable by itself.

	increase_stack_limit(1024);

	int grid_size;
	//float birth_probability;
	//int r_init;
	int how_many;
	//int lag;
	//int number_of_census;
	int divisions;
	int t_to_eq;

	/* cout << "Enter grid size: ";
  cin >> grid_size;
  cout << "Enter p value of contention: ";
  cin >> birth_probability;
  cout << "Enter number of census: ";
  cin >> number_of_census;
  cout << "Enter lag: ";
  cin >> lag;
	cout << "Enter number of random trials: ";
  cin >> r_init; */

	t_to_eq = 50000; grid_size = 512; //birth_probability =0.704; r_init =16; lag=1;
	how_many = 1000000;

	//Taken from: http://www.cplusplus.com/forum/unices/112048/

	float data[40][7];
  std::ifstream file("dump/15_16_KungF---U.csv");

	file.ignore(140, '\n'); //ignore the first 140 characters, or until first \n, whichever is met first

  for(int row = 0; row < 40; ++row)
  {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < 7; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> data[row][col];
        }
  }

	std::vector<f_coordinates> pL; std::vector<transformation> transformations;

	for(int i=32 ; i< 36; i++)
	{
		for(int j=1; j <7; j++)
		{
			f_coordinates temp; temp.x = data[i][0]; temp.y =data[i][j]*(grid_size*grid_size); pL.push_back(temp);
		}
	}
	divisions = 4*6;
	cout << "p\t|  lag\n";
	for (int i=0; i < divisions; i++)
	{
			cout << pL[i].x << "  " << pL[i].y <<endl;
	}



	find_equilibrium_multi_shot_transformations_DP(grid_size, divisions, t_to_eq, how_many, transformations, pL);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Total CPU Time: " << duration.count() << " seconds" << endl;
  }
