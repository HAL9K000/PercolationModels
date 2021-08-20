#include "cluster_dynamics.h"

//---------------------Some Primitives -------------------------------------------//

double random_real(double initial, double last) {

    std::uniform_real_distribution<double> distribution(initial, last);
    return distribution(rng);  // Use rng as a generator
}

int random_int(int initial, int last) {

    std::uniform_int_distribution<int> distribution(initial, last);
    return distribution(rng);  // Use rng as a generator
}

float mean_of_array(float array[],int size){

	float sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(float)size;
}

float standard_deviation_of_array(float array[],int size){

	float mean = mean_of_array(array,size);
	float sum = 0.0;

	for (int i = 0; i < size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}

	float variance = sum/(float)size;

	return sqrt(variance);
}

float mean_of_vector(vector<float> array,int size){

	float sum = 0.0;

	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(float)size;
}

void random_frame(int frame[], int grid_size) {

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {

            frame[i*grid_size + j] = random_int(0, 1);

        }
    }
}

void random_frame_of_density(float density, int frame[], int grid_size) {

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {

        	if (random_real(0,1)<=density){
        		frame[i*grid_size + j] = 1;
        	}
        	else{
        		frame[i*grid_size + j] = 0;
        	}
        }
    }
}

void zeros(int frame[], int grid_size) {

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            frame[i*grid_size + j] = 0;
        }
    }
}

float calculate_density(int frame[], int grid_size){

	float occupancy = 0;

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            occupancy += frame[i*grid_size + j];
        }
    }

    float density = occupancy/(grid_size*grid_size);
	return density;
}

f_coordinates calculate_SD(int frame[],int grid_size)
{
  //Returns tuple with mean and SD.
  float mean = calculate_density(frame, grid_size);
	float sum = 0.0;

	for (int i = 0; i < grid_size*grid_size; ++i){
		sum += (frame[i]-mean)*(frame[i]-mean);
	}

	float variance = sum/(float)(grid_size*grid_size);

  f_coordinates mean_sd;
  mean_sd.x = mean; mean_sd.y = sqrt(variance);

	return mean_sd;
}

void increase_stack_limit(int stack_size){

	//Source: https://stackoverflow.com/questions/2279052/increase-stack-size-in-linux-with-setrlimit
	//Credit: https://stackoverflow.com/users/253056/paul-r
	//Used with modification

	// This becomes necessary for recursive depth first search when finding clusters.

	const rlim_t kStackSize = 64 * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }
}

int disjoint_vectors(vector<int> vec_1, vector<int> vec_2){

	//https://www.faceprep.in/c/check-if-the-given-arrays-are-disjoint-in-c-c-java-and-python-faceprep/

	// returns 1 if the two arrays have a common element, 0 otherwise

    int size_vec_1 = vec_1.size();
    int size_vec_2 = vec_2.size();

    for(int i = 0;i<size_vec_1;i++)
    {
        for(int j=0;j<size_vec_2;j++)
        {
            if(vec_1[i] == vec_2[j])
                return 1;
        }
    }
    return 0;
}

vector<int> common_elements(vector<int> vec_1, vector<int> vec_2){

	//Credits: https://www.faceprep.in/c/check-if-the-given-arrays-are-disjoint-in-c-c-java-and-python-faceprep/

	//returns the unique common elements of two vectors

	vector<int> common_elements;

    int size_vec_1 = vec_1.size();
    int size_vec_2 = vec_2.size();

    for(int i = 0;i<size_vec_1;i++)
    {
        for(int j=0;j<size_vec_2;j++)
        {
            if(vec_1[i] == vec_2[j])
                common_elements.push_back(vec_1[i]);
        }
    }

    sort( common_elements.begin(), common_elements.end() );
	common_elements.erase( unique( common_elements.begin(), common_elements.end() ), common_elements.end() );

    return common_elements;
}

template<typename T>std::vector<double> linspace(T start_in, T end_in, int num_in){

	//Source: https://stackoverflow.com/questions/27028226/python-linspace-in-c/27030598#27030598
	//Credits: https://stackoverflow.com/users/1078084/akavall


	//Equivalent of numpy's linspace method.

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);
  return linspaced;
}

void print_vector(std::vector<double> vec){

  std::cout << "size: " << vec.size() << std::endl;
  for (double d : vec)
    std::cout << d << " ";
  std::cout << std::endl;
}

// -------------------------------Percolation Models---------------------------------------//

coordinates select_neighbor_of_site(coordinates site, int grid_size){

	// Selects at random one of the four neighbours in von-neumann radius of 1 with periodic boundary conditions

	int x = site.x;
	int y = site.y;

	coordinates neighbours[4];

	neighbours[0].x = (x-1)%grid_size, neighbours[0].y = y; //left
	neighbours[1].x = (x+1)%grid_size, neighbours[1].y = y; //right
	neighbours[2].x = x, neighbours[2].y = (y-1)%grid_size; //bottom
	neighbours[3].x = x, neighbours[3].y = (y+1)%grid_size; //top

	int who = random_int(0,3);

	coordinates neighbour = neighbours[who];

	if (neighbour.x < 0){
		neighbour.x = grid_size-1; //correction if selected site is the left of first column
	}

	if (neighbour.y < 0){
		neighbour.y = grid_size-1; //correction if selected site is the top of first row
	}

	return neighbour;
}

coordinates select_neighbor_of_pair(coordinates site, coordinates neighbour, int grid_size){

	// Selects at random one of the six neighbours of a pair of sites. Here too we have periodic boundary conditions.

	int diff_x = site.x - neighbour.x;
	int diff_y = site.y - neighbour.y;

	coordinates neighbours_of_pair[6];

	if ((diff_x ==1||diff_x==-(grid_size-1)) && diff_y ==0){ // neighbour is to the left of focal site

		neighbours_of_pair[0].x = (site.x+1)%grid_size, neighbours_of_pair[0].y = site.y;
		neighbours_of_pair[1].x = (neighbour.x-1)%grid_size, neighbours_of_pair[1].y = site.y;
		neighbours_of_pair[2].x = site.x, neighbours_of_pair[2].y = (site.y+1)%grid_size;
		neighbours_of_pair[3].x = site.x, neighbours_of_pair[3].y = (site.y-1)%grid_size;
		neighbours_of_pair[4].x = neighbour.x, neighbours_of_pair[4].y = (neighbour.y+1)%grid_size;
		neighbours_of_pair[5].x = neighbour.x, neighbours_of_pair[5].y = (neighbour.y-1)%grid_size;
	}
	if ((diff_x ==-1 ||diff_x==(grid_size-1)) && diff_y ==0){ // neighbour is to the right of focal site

		neighbours_of_pair[0].x = (site.x-1)%grid_size, neighbours_of_pair[0].y = site.y;
		neighbours_of_pair[1].x = (neighbour.x+1)%grid_size, neighbours_of_pair[1].y = site.y;
		neighbours_of_pair[2].x = site.x, neighbours_of_pair[2].y = (site.y+1)%grid_size;
		neighbours_of_pair[3].x = site.x, neighbours_of_pair[3].y = (site.y-1)%grid_size;
		neighbours_of_pair[4].x = neighbour.x, neighbours_of_pair[4].y = (neighbour.y+1)%grid_size;
		neighbours_of_pair[5].x = neighbour.x, neighbours_of_pair[5].y = (neighbour.y-1)%grid_size;
	}
	if (diff_x ==0 && (diff_y ==1 ||diff_y==-(grid_size-1))){ // neighbour is below the focal site

		neighbours_of_pair[0].x = site.x, neighbours_of_pair[0].y = (site.y+1)%grid_size;
		neighbours_of_pair[1].x = neighbour.x, neighbours_of_pair[1].y = (neighbour.y-1)%grid_size;
		neighbours_of_pair[2].x = (site.x+1)%grid_size, neighbours_of_pair[2].y = site.y;
		neighbours_of_pair[3].x = (site.x-1)%grid_size, neighbours_of_pair[3].y = site.y;
		neighbours_of_pair[4].x = (neighbour.x+1)%grid_size, neighbours_of_pair[4].y = neighbour.y;
		neighbours_of_pair[5].x = (neighbour.x-1)%grid_size, neighbours_of_pair[5].y = neighbour.y;
	}
	if (diff_x ==0 && (diff_y ==-1 ||diff_y==(grid_size-1))){ // neighbour is above the focal site

		neighbours_of_pair[0].x = site.x, neighbours_of_pair[0].y = (site.y-1)%grid_size;
		neighbours_of_pair[1].x = neighbour.x, neighbours_of_pair[1].y = (neighbour.y+1)%grid_size;
		neighbours_of_pair[2].x = (site.x+1)%grid_size, neighbours_of_pair[2].y = site.y;
		neighbours_of_pair[3].x = (site.x-1)%grid_size, neighbours_of_pair[3].y = site.y;
		neighbours_of_pair[4].x = (neighbour.x+1)%grid_size, neighbours_of_pair[4].y = neighbour.y;
		neighbours_of_pair[5].x = (neighbour.x-1)%grid_size, neighbours_of_pair[5].y = neighbour.y;
	}

	int who = random_int(0,5);

	coordinates neighbour_of_pair = neighbours_of_pair[who];

	if (neighbour_of_pair.x < 0){
		neighbour_of_pair.x = grid_size-1; //correction if selected site is the left of first column
	}

	if (neighbour_of_pair.y < 0){
		neighbour_of_pair.y = grid_size-1; //correction if selected site is the top of first row
	}

	return neighbour_of_pair;
}

void np_update(int frame[], int grid_size, float birth_probability){

	// performs single asynchronous update for P at a randomly selected site on the lattice

	int x = random_int(0, grid_size-1);
	int y = random_int(0, grid_size-1);

	double chance = random_real(0, 1);

	if (frame[x*grid_size+y]==0){ // selected site is occupied

		if (chance < birth_probability){
			frame[x*grid_size+y]=1; // birth at empty focal site with probability p
		}
	}
	else {

		if (chance < 1-birth_probability){
			frame[x*grid_size+y]=0; // death at occupied focal site with probability (1-p)
		}
	}
}

void dp_update(int frame[], int grid_size, float birth_probability){

	// performs single asynchronous update for DP at a randomly selected site on the lattice

	int x = random_int(0, grid_size-1);
	int y = random_int(0, grid_size-1);

	coordinates site;
	site.x = x;
	site.y = y;

	if (frame[x*grid_size+y]==1){ // selected site is occupied

		coordinates neighbour = select_neighbor_of_site(site,grid_size);

		double chance = random_real(0, 1);

		if (chance<birth_probability){

			frame[neighbour.x*grid_size+neighbour.y] = 1; //birth at neighbour of occupied focal site with probability p
		}
		else {

			frame[x*grid_size+y]=0; //death at focal site with probability (1-p_)
		}
	}
}

void tp_update(int frame[], int grid_size, float birth_probability, float feedback_strength){

	// performs single asynchronous update for TP at a randomly selected site on the lattice

	int x = random_int(0, grid_size-1);
	int y = random_int(0, grid_size-1);

	coordinates site;
	site.x = x;
	site.y = y;

	if (frame[x*grid_size+y]==1){ // selected site is occupied

		coordinates neighbour = select_neighbor_of_site(site,grid_size);

		if (frame[neighbour.x*grid_size+neighbour.y]==0) { // neighbour of focal site is empty. DP like updates in this case

			double chance = random_real(0, 1);

			if (chance < birth_probability){
				frame[neighbour.x*grid_size+neighbour.y] = 1; //birth at neighbour of occupied focal site with probability p
			}
			else {
				frame[x*grid_size+y]=0; //death at focal site with probability (1-p_)
			}
		}
		else { // neighbour of focal site is occupied. Facilitation taken into account

			double chance = random_real(0, 1);
			double another_chance = random_real(0, 1);

			if (chance < feedback_strength){

				coordinates neighbour_of_pair = select_neighbor_of_pair(site,neighbour,grid_size);
				frame[neighbour_of_pair.x*grid_size+neighbour_of_pair.y] = 1; // birth at neighbour of pair with enhanced probability q
			}
			else if (another_chance < 1-birth_probability){

				frame[x*grid_size+y]=0; // death at focal site with diminished death probability (1-q)(1-p)
			}
		}
	}
}

void simulate_np(int frame[], int grid_size, float birth_probability, int updates_per_site){

	// Takes a frame and simulates P for the specified number of updates per site with the specified birth probability

	while (updates_per_site > 0) {

		for (int i = 0; i < grid_size*grid_size; i++){ // loop performs on average a single update per site

			np_update(frame, grid_size, birth_probability);
		}
		updates_per_site -= 1;
	}
}

void simulate_dp(int frame[], int grid_size, float birth_probability, int updates_per_site){

	// Takes a frame and simulates DP for the specified number of updates per site with the specified birth probability

	while (updates_per_site > 0) {

		for (int i = 0; i < grid_size*grid_size; i++){ // loop performs on average a single update per site

			dp_update(frame, grid_size, birth_probability);
		}

		updates_per_site -= 1;
	}
}

void simulate_tp(int frame[], int grid_size, float birth_probability, float feedback_strength, int updates_per_site){

	// Takes a frame and simulates TP for the specified number of updates per site with the specified birth probability

	while (updates_per_site > 0) {

		for (int i = 0; i < grid_size*grid_size; i++){ // loop performs on average a single update per site

			tp_update(frame, grid_size, birth_probability, feedback_strength);
		}
		updates_per_site -= 1;
	}
}

float equilibrium_density_dp(int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site, int collect_frames){

	// Simulates DP with specified parameters and prints the mean and standard deviation of vegetation cover to the terminal. The frames from
	// which vegetation cover is calculated can also be captured. Finally it returns the mean vegetation cover.

	float densities[number_of_census];
	int frame[grid_size*grid_size];
    random_frame(frame, grid_size);

    simulate_dp(frame,grid_size,birth_probability,updates_per_site);

    for (int i=0; i<number_of_census; i++) {

    	simulate_dp(frame,grid_size,birth_probability,lag);
    	densities[i] = calculate_density(frame,grid_size);

    	if (collect_frames ==1){ // conditional for collecting frames

    		ofstream outdata;

			stringstream p, census, Lag;
			p << setprecision(3) << birth_probability;
			census << i;
			Lag << lag;

			outdata.open("dump/dp_frame_"+std::to_string(grid_size)+"_"+p.str()+"_"+Lag.str()+'_'+census.str()+".txt");
			if( !outdata ) {
				cerr << "File error, try again." << endl;
				exit(1);
			}
			for (int i = 0; i < grid_size; ++i) {
		        for (int j = 0; j < grid_size; ++j) {
		            outdata << frame[i*grid_size + j];
		        }
		        outdata << '\n';
		    }
			outdata.close();
    	}
    }

    float mean_density = mean_of_array(densities,number_of_census);
    float standard_deviation_density = standard_deviation_of_array(densities,number_of_census);

    cout << endl;
    cout << "p: " << setprecision(4) << birth_probability << " Mean: " << setprecision(2) << mean_density << " Standard Deviation: " << setprecision(4) << standard_deviation_density << endl;
    return mean_density;
}

float equilibrium_density_tp(int grid_size, float birth_probability, float feedback_strength, int number_of_census, int lag, int updates_per_site, int collect_frames){

	// Simulates TP with specified parameters and prints the mean and standard deviation of vegetation cover to the terminal. The frames from
	// which vegetation cover is calculated can also be captured. Finally it returns the mean vegetation cover.

	float densities[number_of_census];
	int frame[grid_size*grid_size];
    random_frame(frame, grid_size);

    simulate_tp(frame,grid_size,birth_probability,feedback_strength, updates_per_site);

    for (int i=0; i<number_of_census; i++) {

    	simulate_tp(frame, grid_size, birth_probability, feedback_strength, lag);
    	densities[i] = calculate_density(frame,grid_size);

    	if (collect_frames ==1){

    		ofstream outdata;

			stringstream p, q, census, Lag;
			p << setprecision(3) << birth_probability;
			q << setprecision(2) << feedback_strength;
			census << i;
			Lag << lag;

			outdata.open("dump/tp_frame_"+std::to_string(grid_size)+"_"+q.str()+'_'+p.str()+"_"+Lag.str()+"_"+census.str()+".txt");
			if( !outdata ) {
				cerr << "File error, try again." << endl;
				exit(1);
			}
			for (int i = 0; i < grid_size; ++i) {
		        for (int j = 0; j < grid_size; ++j) {
		            outdata << frame[i*grid_size + j];
		        }
		        outdata << '\n';
		    }
			outdata.close();
    	}
    }

    float mean_density = mean_of_array(densities,number_of_census);
    float standard_deviation_density = standard_deviation_of_array(densities,number_of_census);

    cout << endl;
    cout << "p: " << setprecision(4) << birth_probability << " q: " << feedback_strength <<" Mean: " << setprecision(2) << mean_density << " Standard Deviation: " << setprecision(2) << standard_deviation_density << endl;
    return mean_density;
}

// -----------------------------------Cluster Statics -----------------------------------------//

vector<int> spanning_cluster_label(int frame[], int grid_size){

	// Takes a frame, labels it, checks if a spanning cluster is present and returns its label.

	vector<int> up; // empty vector to store upper border
	vector<int> down; // empty vector to store lower border
	vector <int> left; // empty vector to store left border
	vector <int> right; // empty vector to store right border

	int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

	vector<int> spanning_cluster_label;
	vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

	find_clusters_free_boundary(frame, labels, clusters, grid_size);
	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

	for (int i = 0; i<grid_size; i++){ // loop to populate up, down, left and right.

		if (labels[i] != 0){
			up.push_back(labels[i]);
		}

		if (labels[grid_size*(grid_size-1)+i] != 0){
			down.push_back(labels[grid_size*(grid_size-1)+i]);
		}

		if (labels[grid_size*i] != 0){
			left.push_back(labels[grid_size*i]);
		}

		if (labels[grid_size*(i+1)-1] != 0){
			right.push_back(labels[grid_size*(i+1)-1]);
		}
	}

	if (is_spanning_horizontal(frame,grid_size)){ // conditional to check if the frame has a cluster that spans the lattice horizontally

		spanning_cluster_label = common_elements(left, right); // label of the horizontally spanning cluster which must be present in both left and right

	}
	else if (is_spanning_vertical(frame,grid_size)){ // conditional to check if the frame has a cluster that spans the lattice vertically

		spanning_cluster_label = common_elements(up, down); // label of the vertically spanning cluster which must be present in both up and down

	}
	else{ // the frame does not have a spanning cluster

		spanning_cluster_label.push_back(-1); // Unique identification of absence of spanning cluster since label must be non-negative integer

	}
	return spanning_cluster_label;
}

void find_clusters_without_spanning(int frame[], int labels[],vector<cluster>& clusters, int grid_size){

	// Accumulates all the finite clusters from the frame

 	int id = 1;

 	vector <int> spanning_label = spanning_cluster_label(frame, grid_size); // find the label of spanning cluster if present

	for(int i=0; i<grid_size;i++){
     	for (int j=0; j<grid_size; j++){  //traverse every site on the lattice once

         	coordinates site;
         	site.x = i;
         	site.y = j;
         	if((frame[i*grid_size+j]==1) && labels[i*grid_size+j]==0){ // if the site is occupied but unlabelled.

             	std::vector<int> sites; //to accumulate all sites that are part of the cluster to which the focal site belongs

             	depth_first_search_free_boundary(site,sites,frame,labels,id++,grid_size); //recursive depth first search to find the cluster's sites

             	cluster this_cluster; // declaring a cluster structure
             	this_cluster.label = labels[i*grid_size+j]; // storing the label of the cluster
             	this_cluster.coords = sites; // storing the coordinates of the cluster

             	if (this_cluster.label != spanning_label[0]){ // accumulate finite clusters only
             		clusters.push_back(this_cluster);
             	}
         	}
     	}
 	}
}

void find_patch_sizes_dp(vector<int>& cluster_sizes,int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site, int collect_frames){

	// Accumulates the cluster sizes from the specified number of frames with the specified parameters

	int frame[grid_size*grid_size];
    random_frame(frame, grid_size);

    simulate_dp(frame,grid_size,birth_probability,updates_per_site);

    for (int i=0; i<number_of_census; i++) {

		simulate_dp(frame,grid_size,birth_probability,lag);

		int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

		vector<cluster> clusters;

		find_clusters_without_spanning(frame, labels, clusters, grid_size); //accumulate all finite clusters from the current frame

		for (int i=0; i<clusters.size(); i++){
			cluster_sizes.push_back(clusters[i].coords.size()); // append the clusters from the current frame to the main cluster sizes vector
		}

		if (collect_frames ==1){ // conditional to collect frames

			ofstream outdata;

			stringstream p, census, Lag;
			p << setprecision(3) << birth_probability;
			census << i;
			Lag << lag;

			outdata.open("dump/dp_frame_"+std::to_string(grid_size)+"_"+p.str()+"_"+Lag.str()+'_'+census.str()+".txt");
			if( !outdata ) {
				cerr << "File error, try again." << endl;
				exit(1);
			}
			for (int i = 0; i < grid_size; ++i) {
		        for (int j = 0; j < grid_size; ++j) {
		            outdata << frame[i*grid_size + j];
		        }
		        outdata << '\n';
		    }
			outdata.close();
		}
    }
}

void find_patch_sizes_np(vector<int>& cluster_sizes,int grid_size, float birth_probability, int number_of_census, int lag, int updates_per_site, int collect_frames){

	// Accumulates the cluster sizes from the specified number of frames with the specified parameters

	int frame[grid_size*grid_size];
    random_frame(frame, grid_size);

    simulate_np(frame,grid_size,birth_probability,updates_per_site);

    for (int i=0; i<number_of_census; i++) {

		simulate_np(frame,grid_size,birth_probability,lag);

		int labels[grid_size*grid_size]; // initialize all sites of label lattice to 0

		for (int i=0; i<grid_size; i++){
			for (int j=0; j<grid_size; j++){
				labels[grid_size*i+j] = 0;
			}
		}

		vector<cluster> clusters;

		find_clusters_without_spanning(frame, labels, clusters, grid_size); //accumulate all finite clusters from the current frame

		for (int i=0; i<clusters.size(); i++){
			cluster_sizes.push_back(clusters[i].coords.size()); // append the clusters from the current frame to the main cluster sizes vector
		}

		if (collect_frames ==1){ // conditional to collect frames

			ofstream outdata;

			stringstream p, census, Lag;
			p << setprecision(3) << birth_probability;
			census << i;
			Lag << lag;

			outdata.open("dump/np_frame_"+std::to_string(grid_size)+"_"+p.str()+"_"+Lag.str()+'_'+census.str()+".txt");
			if( !outdata ) {
				cerr << "File error, try again." << endl;
				exit(1);
			}
			for (int i = 0; i < grid_size; ++i) {
		        for (int j = 0; j < grid_size; ++j) {
		            outdata << frame[i*grid_size + j];
		        }
		        outdata << '\n';
		    }
			outdata.close();
		}
    }
}

//------------------------------------Cluster Dynamics----------------------------------------//

void find_neighbours_of_site(coordinates neighbours[4], coordinates site, int grid_size){

	// returns all the four neighbours in the von Newmann neighbourhood of the focal site as an array of coordinates. See cluster_dynamics.h for the structure coordinates

	int x = site.x;
	int y = site.y;

	neighbours[0].x = (x-1)%grid_size, neighbours[0].y = y; // left neighbour
	neighbours[1].x = (x+1)%grid_size, neighbours[1].y = y; // right neighbour
	neighbours[2].x = x, neighbours[2].y = (y-1)%grid_size; // bottom neighbour
	neighbours[3].x = x, neighbours[3].y = (y+1)%grid_size; // top neighbour

	for(int i=0; i<4; i++){
		if (neighbours[i].x < 0){
		neighbours[i].x = grid_size-1; //correction if selected site is the left of first column
		}

		if (neighbours[i].y < 0){
		neighbours[i].y = grid_size-1; //correction if selected site is the top of first row
		}
	}
}

void depth_first_search(coordinates site, vector<int>& sites, int frame[], int labels[], int id, int grid_size){

	//Source: https://stackoverflow.com/questions/22051069/how-do-i-find-the-connected-components-in-a-binary-image
	//Credits: https://stackoverflow.com/users/2185825/shole
	//Used with modification

	// Recursively labels a cluster and accumulates its sites

    labels[grid_size*site.x+site.y] = id; // Label the focal site on the labels lattice

    sites.push_back(grid_size*site.x+site.y); // Add the focal site to the sites vector

    coordinates neighbours[4];

    find_neighbours_of_site(neighbours, site, grid_size); // find the neighbours of the focal site with periodic boundary conditions

    for(int i=0; i<4;i++){ // iterate through the neighbours of the focal site

        if((frame[grid_size*neighbours[i].x + neighbours[i].y]==1) && (labels[grid_size*neighbours[i].x+neighbours[i].y]==0)){

    		//conditional entered if the neighbour is occupied but unlaballed

         	depth_first_search(neighbours[i],sites,frame,labels,id,grid_size);
         	// recursively call depth first search so that the neighbour becomes the focal site
        }
    }
}

void find_clusters(int frame[], int labels[],vector<cluster>& clusters, int grid_size){

	// Find all the clusters (labels and coordinates) in the frame

 	int id = 1;

 	for(int i=0; i<grid_size;i++){
     	for (int j=0; j<grid_size; j++){ //traverse the lattice

         	coordinates site;
         	site.x = i;
         	site.y = j;
         	if((frame[i*grid_size+j]==1) && labels[i*grid_size+j]==0){ // if the site is occupied but unlabelled

             	std::vector<int> sites;

             	depth_first_search(site,sites,frame,labels,id++,grid_size); //label the cluster to which the current site belongs and accumulate its coordinates

             	cluster this_cluster;
             	this_cluster.label = labels[i*grid_size+j]; // label of the cluster to which the current site belongs
             	this_cluster.coords = sites; // coordinates of the cluster to which the current site belongs

             	clusters.push_back(this_cluster);
         	}
     	}
 	}
}

int check_presence(std::vector<int> vec,int element){

	// Returns 1 if the element is present, 0 otherwise

	std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), element);

	if (it != vec.end()){
		return 1;
	}
	else {
		return 0;
	}
}

void find_transformations_single_shot(vector<transformation>& transformations,
	int previous_frame[], int previous_frame_labels[],vector<cluster>& previous_frame_clusters,
	int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size){

	// takes two consecutive frames with the difference of a single update and finds the size of the cluster before and after the update.
	// The change in size is the difference of these two cluster sizes.

	int diff[grid_size*grid_size];
	std::vector<int> changes;

	for (int i=0; i<grid_size*grid_size; i++){ // scan the lattice to find the change
	    diff[i] = current_frame[i] - previous_frame[i];
	    if (diff[i] != 0){
	        changes.push_back(i); // collect the site at which the change has occured
	    }
	}


	if (changes.size()==0 || changes.size() >1){ // exit function if no change or more than one change has occured
     	return;
 	}
 	else{

     	std::vector<int> finished_labels_current_frame;
     	std::vector<int> finished_labels_previous_frame;

     	if (diff[changes[0]]>0){ // the change was a birth

     		// it may have: a) created a cluster of size 1 b) caused a cluster of size n to become a cluster of size n+1 c) caused two or clusters to merge

     		if (!check_presence(finished_labels_current_frame,current_frame_labels[changes[0]])){

     			// conditional is entered if the label of the site which appeared in the current frame is not present in finished labels of the current frame

     			finished_labels_current_frame.push_back(current_frame_labels[changes[0]]); // add the label to finished labels of current frame

                cluster focal_cluster = current_frame_clusters[current_frame_labels[changes[0]]-1]; // get the cluster to which the newly appeared site belongs

                std::vector<int> comparison_clusters;

                std::vector<int> finished_labels_comparison_clusters;

                for (int j=0; j<focal_cluster.coords.size(); j++){ // loop through the sites of the focal cluster
                    if ((previous_frame_labels[focal_cluster.coords[j]] != 0)&&
                    	(!check_presence(finished_labels_comparison_clusters,previous_frame_labels[focal_cluster.coords[j]]))){

                    	// conditional is entered if the label of the cluster in the previous frame at the site
                    	// corresponding to the current site of the focal cluster (from the current frame)
                    	// is absent in finished_labels_comparison_clusters

                        finished_labels_comparison_clusters.push_back(previous_frame_labels[focal_cluster.coords[j]]);
                        //add the label from the previous frame to finished_labels_comparison_clusters

                        comparison_clusters.push_back(previous_frame_clusters[previous_frame_labels[focal_cluster.coords[j]]-1].coords.size());
                        // add the cluster from the previous frame to comparison_clusters
                    }
                }

                transformation transformation_merge; // See cluster_dynamics.h for the definition of the transformation structure

                if (comparison_clusters.size() == 0){ //This means the a cluster of size 1 has appeared.

                    transformation_merge.before = 0;
                    transformation_merge.after = focal_cluster.coords.size();
                    transformations.push_back(transformation_merge);

                }
                else { // This means that either a cluster has grown by size 1 or two or more clusters have merged.

                	// If two or more clusters have merged we only consider the largest of the comparison clusters from the previous frame
                	// as the cluster which has transformed to the single cluster in the current frame

                    transformation_merge.before = *max_element(comparison_clusters.begin(), comparison_clusters.end()); //finds the largest among the comparison clusters
                    transformation_merge.after = focal_cluster.coords.size();
                    transformations.push_back(transformation_merge);

                }
            }
     	}
     	else { // the change was a death

     		// it may have: a) removed a cluster of size 1 b) caused a cluster of size n to become a cluster of size n-1 c) caused two or clusters to split


     		if (!check_presence(finished_labels_previous_frame,previous_frame_labels[changes[0]])){

     			// conditional is entered if the label of the site which disappeared in the current frame is not present in finished labels of the previous frame

     			finished_labels_previous_frame.push_back(previous_frame_labels[changes[0]]); // add the label to finished labels of previous frame

     			cluster focal_cluster = previous_frame_clusters[previous_frame_labels[changes[0]]-1]; // get the cluster in the previous frame from which a site has disappeared

     			std::vector<int> comparison_clusters;

                std::vector<int> finished_labels_comparison_clusters;

                for (int j=0; j<focal_cluster.coords.size(); j++){ // loop through the sites of the focal cluster

	                if ((current_frame_labels[focal_cluster.coords[j]] != 0)&&
	                	(!check_presence(finished_labels_comparison_clusters,current_frame_labels[focal_cluster.coords[j]]))){

	                	// conditional is entered if the label of the cluster in the current frame at the site
                    	// corresponding to the current site of the focal cluster (from the previous frame)
                    	// is absent in finished_labels_comparison_clusters

	                    finished_labels_comparison_clusters.push_back(current_frame_labels[focal_cluster.coords[j]]);
	                	//add the label from the current frame to finished_labels_comparison_clusters

	                    comparison_clusters.push_back(current_frame_clusters[current_frame_labels[focal_cluster.coords[j]]-1].coords.size());
	                    // add the cluster from the current frame to comparison_clusters
	                }
                }

                transformation transformation_primary_split; // See cluster_dynamics.h for the definition of the transformation structure

                if (comparison_clusters.size() == 0){ //This means the a cluster of size 1 has disappeared.

	                transformation_primary_split.before = focal_cluster.coords.size();
	                transformation_primary_split.after = 0;
	                transformations.push_back(transformation_primary_split);

               	}
				else { // This means that either a cluster has shrunk by size 1 or a cluster has split into two or more clusters

					// If a cluster has split into two or more clusters we only consider the largest of the comparison clusters
					//from the curernt frame as the cluster to which the single cluster from the previous frame has transformed.

					transformation_primary_split.before = focal_cluster.coords.size();
					transformation_primary_split.after = *max_element(comparison_clusters.begin(), comparison_clusters.end()); //finds the largest among the comparison clusters
					transformations.push_back(transformation_primary_split);
				}
     		}
     	}
    }
}

void find_equilibrium_single_shot_transformations_np(int grid_size, float birth_probability, int time_to_equilibrium,
	int how_many, vector<transformation>& transformations){

	// Simulates P for the specified parameters and collects the specified number of transformations from the steady state

	int seed = std::random_device{}();
	rng.seed(seed);


	int previous_frame[grid_size*grid_size];
	int previous_frame_labels[grid_size*grid_size];
	vector<cluster> previous_frame_clusters;

	int current_frame[grid_size*grid_size];
	int current_frame_labels[grid_size*grid_size];
	vector<cluster> current_frame_clusters;

	random_frame(current_frame, grid_size); // initialize the current frame
	zeros(current_frame_labels,grid_size); // initialize labels for the current frame

	simulate_np(current_frame, grid_size, birth_probability, time_to_equilibrium); // simulate P to reach steady state

	auto start = high_resolution_clock::now();

	find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size); // find clusters from the current frame

	while (transformations.size() < how_many){ // loop till we have the specified number of transformations

		for (int i=0; i<grid_size; i++){
			for (int j=0; j<grid_size; j++){ // scan the current frame and duplicate it into the previous frame
				previous_frame[grid_size*i+j] = current_frame[grid_size*i+j];
				previous_frame_labels[grid_size*i+j] = current_frame_labels[grid_size*i+j];
			}
		}

		previous_frame_clusters = current_frame_clusters; // duplicate the clusters too

		np_update(current_frame, grid_size, birth_probability); // make a single update to the current frame

		current_frame_clusters.clear(); // wash the stale clusters of the current frame (since it has been updated now)
		zeros(current_frame_labels,grid_size); // wash the labels too

		find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size);
		// find the clusters and labels of the current frame (after update) again

		find_transformations_single_shot(transformations,
			previous_frame, previous_frame_labels, previous_frame_clusters,
			current_frame, current_frame_labels, current_frame_clusters, grid_size); // find the transformation and accumulate into the transformations vector
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Transformation Finding Time: " << duration.count() << " seconds" << endl;
}

void find_equilibrium_single_shot_transformations_dp(int grid_size, float birth_probability, int time_to_equilibrium,
	int how_many, vector<transformation>& transformations){

	// Simulates DP for the specified parameters and collects the specified number of transformations from the steady state

	int seed = std::random_device{}();
	rng.seed(seed);


	int previous_frame[grid_size*grid_size];
	int previous_frame_labels[grid_size*grid_size];
	vector<cluster> previous_frame_clusters;

	int current_frame[grid_size*grid_size];
	int current_frame_labels[grid_size*grid_size];
	vector<cluster> current_frame_clusters;

	random_frame(current_frame, grid_size); // initialize the current frame
	zeros(current_frame_labels,grid_size); // initialize labels for the current frame

	auto start = high_resolution_clock::now();

	simulate_dp(current_frame, grid_size, birth_probability, time_to_equilibrium); // simulate DP to reach steady state

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Simulation Time: " << duration.count() << " seconds" << endl;

	start = high_resolution_clock::now();

	find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size); // find clusters from the current frame

	while (transformations.size() < how_many){ // loop till we have the specified number of transformations

		for (int i=0; i<grid_size; i++){ // scan the current frame and duplicate it into the previous frame
			for (int j=0; j<grid_size; j++){
				previous_frame[grid_size*i+j] = current_frame[grid_size*i+j];
				previous_frame_labels[grid_size*i+j] = current_frame_labels[grid_size*i+j];
			}
		}

		previous_frame_clusters = current_frame_clusters; // duplicate the clusters too

		dp_update(current_frame, grid_size, birth_probability); // make a single update to the current frame

		current_frame_clusters.clear(); // wash the stale clusters of the current frame (since it has been updated now)
		zeros(current_frame_labels,grid_size); // wash the labels too

		find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size);
		// find the clusters and labels of the current frame (after update) again

		find_transformations_single_shot(transformations,
			previous_frame, previous_frame_labels, previous_frame_clusters,
			current_frame, current_frame_labels, current_frame_clusters, grid_size); // find the transformation and accumulate into the transformations vector
	}

	stop = high_resolution_clock::now();
	duration = duration_cast<seconds>(stop - start);

	cout << endl << "Transformation Finding Time: " << duration.count() << " seconds" << endl;
}

void find_equilibrium_single_shot_transformations_tp(int grid_size, float birth_probability, float feedback_strength, int time_to_equilibrium,
	int how_many, vector<transformation>& transformations){

	// Simulates TP for the specified parameters and collects the specified number of transformations from the steady state

	int seed = std::random_device{}();
	rng.seed(seed);


	int previous_frame[grid_size*grid_size];
	int previous_frame_labels[grid_size*grid_size];
	vector<cluster> previous_frame_clusters;

	int current_frame[grid_size*grid_size];
	int current_frame_labels[grid_size*grid_size];
	vector<cluster> current_frame_clusters;

	random_frame(current_frame, grid_size); // initialize the current frame
	zeros(current_frame_labels,grid_size); // initialize labels for the current frame

	simulate_tp(current_frame, grid_size, birth_probability, feedback_strength, time_to_equilibrium); // simulate TP to reach steady state

	auto start = high_resolution_clock::now();

	find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size); // find clusters from the current frame

	while (transformations.size() < how_many){ // loop till we have the specified number of transformations

		for (int i=0; i<grid_size; i++){
			for (int j=0; j<grid_size; j++){ // scan the current frame and duplicate it into the previous frame
				previous_frame[grid_size*i+j] = current_frame[grid_size*i+j];
				previous_frame_labels[grid_size*i+j] = current_frame_labels[grid_size*i+j];
			}
		}

		previous_frame_clusters = current_frame_clusters; // duplicate the clusters too

		tp_update(current_frame, grid_size, birth_probability, feedback_strength); // make a single update to the current frame

		current_frame_clusters.clear(); // wash the stale clusters of the current frame (since it has been updated now)
		zeros(current_frame_labels,grid_size); // wash the labels too

		find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size);
		// find the clusters and labels of the current frame (after update) again

		find_transformations_single_shot(transformations,
			previous_frame, previous_frame_labels, previous_frame_clusters,
			current_frame, current_frame_labels, current_frame_clusters, grid_size); // find the transformation and accumulate into the transformations vector
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << endl << "Transformation Finding Time: " << duration.count() << " seconds" << endl;
}

// -----------------------------------This and That ------------------------------------------------//

void parallelizer(int processes, int grid_size, float birth_probabilities[], float feedback_strength, int number_of_census, int lag, int updates_per_site, int collect_frames) {

    // A parallel loop to find densities at specified birth probalities for TP. This structure can be used to parallelize other functions.

	#pragma omp parallel for
	for (int i=0; i<processes; i++) {

		int seed = std::random_device{}();
		rng.seed(seed);

		int frame[grid_size*grid_size];
        random_frame(frame, grid_size);

    	equilibrium_density_tp(grid_size, birth_probabilities[i], feedback_strength, number_of_census, lag, updates_per_site, collect_frames);
	}
}


void find_neighbours_of_site_free_boundary(coordinates neighbours[4], coordinates site, int grid_size){

	// returns all the four neighbours in the von Newmann neighbourhood of the focal site as an array of coordinates. See cluster_dynamics.h for the structure coordinates
	// Similar to find_neighbours_of_site but with free boundary conditions instead of periodic boundary conditions

	int x = site.x;
	int y = site.y;

	int _site = x*grid_size + y;

	if (_site%grid_size==0){ // if the focal site is on the left border
		neighbours[0].x = -1, neighbours[0].y = -1; // the left neighbour is undefined
	}
	else{
		neighbours[0].x = x, neighbours[0].y = y-1; // the left neighbour
	}

	if (_site%grid_size==grid_size-1){ // if the focal site is on the bottom border
		neighbours[1].x = -1, neighbours[1].y = -1; // the bottom neighbour is undefined
	}
	else{
		neighbours[1].x = x, neighbours[1].y = y+1; // the bottom neighbour
	}

	if (_site < grid_size){ // if the focal site is on the top border
		neighbours[2].x = -1, neighbours[2].y = -1; // the top neighbour is undefined
	}
	else{
		neighbours[2].x = x-1, neighbours[2].y = y; // the top neighbour
	}

	if (_site >= (grid_size-1)*grid_size){ // if the focal site is on the right border
		neighbours[3].x = -1, neighbours[3].y = -1; // the right neighbour is undefined
	}
	else{
		neighbours[3].x = x+1, neighbours[3].y = y; // the right neighbour
	}
}

void depth_first_search_free_boundary(coordinates site, vector<int>& sites, int frame[], int labels[], int id, int grid_size){

	//Source: https://stackoverflow.com/questions/22051069/how-do-i-find-the-connected-components-in-a-binary-image
	//Credits: https://stackoverflow.com/users/2185825/shole
	//Used with modification

	// Recursively labels a cluster and accumulates its sites
	// Similar to depth first search but with free boundary conditions

    labels[grid_size*site.x+site.y] = id; // Label the focal site on the labels lattice

    sites.push_back(grid_size*site.x+site.y); // Add the focal site to the sites vector

    coordinates neighbours[4];

    find_neighbours_of_site_free_boundary(neighbours, site, grid_size); // find the neighbours of the focal site with free boundary conditions

    for(int i=0; i<4;i++){

    	if ((neighbours[i].x >= 0)&&(neighbours[i].y >= 0)){ // iterate through the neighbours of the focal site

    		if((frame[grid_size*neighbours[i].x + neighbours[i].y]==1) && (labels[grid_size*neighbours[i].x+neighbours[i].y]==0)){

    			//conditional entered if the neighbour is occupied but unlaballed

         		depth_first_search_free_boundary(neighbours[i],sites,frame,labels,id,grid_size);
         		// recursively call depth first search so that the neighbour becomes the focal site

    		}
        }
    }
}

void find_clusters_free_boundary(int frame[], int labels[],vector<cluster>& clusters, int grid_size){

	// Find all the clusters (labels and coordinates) in the frame
	// Similar to find_clusters but with free boundary conditions

 	int id = 1;

	for(int i=0; i<grid_size;i++){
     	for (int j=0; j<grid_size; j++){ //traverse the lattice

         	coordinates site;
         	site.x = i;
         	site.y = j;
         	if((frame[i*grid_size+j]==1) && labels[i*grid_size+j]==0){ // if the site is occupied but unlabelled

             	std::vector<int> sites;

             	depth_first_search_free_boundary(site,sites,frame,labels,id++,grid_size); //label the cluster to which the current site belongs and accumulate its coordinates

             	cluster this_cluster;
             	this_cluster.label = labels[i*grid_size+j]; // label of the cluster to which the current site belongs
             	this_cluster.coords = sites; // coordinates of the cluster to which the current site belongs

             	clusters.push_back(this_cluster);
         	}
     	}
 	}
}

int is_spanning_vertical(int frame[], int grid_size){

	// Returns 1 if a cluster spans the lattice vertically, 0 otherwise

	vector<int> up;  // empty vector to store upper border
	vector<int> down;  // empty vector to store bottom border

	int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

	vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

	find_clusters_free_boundary(frame, labels, clusters, grid_size);
	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions


	for (int i = 0; i<grid_size; i++){ // loop to populate up and down

		if (labels[i] != 0){
			up.push_back(labels[i]);
		}

		if (labels[grid_size*(grid_size-1)+i] != 0){
			down.push_back(labels[grid_size*(grid_size-1)+i]);
		}
	}

	return disjoint_vectors(up, down); // check if a cluster spans the lattice vertically and return result
}

int is_spanning_horizontal(int frame[], int grid_size){

	// Returns 1 if a cluster spans the lattice horizontally, 0 otherwise

	vector <int> left; // empty vector to store left border
	vector <int> right; // empty vector to store right border

	int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

	vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

	find_clusters_free_boundary(frame, labels, clusters, grid_size);
	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions


	for (int i = 0; i<grid_size; i++){ // loop to populate left and right

		if (labels[grid_size*i] != 0){
			left.push_back(labels[grid_size*i]);
		}

		if (labels[grid_size*(i+1)-1] != 0){
			right.push_back(labels[grid_size*(i+1)-1]);
		}

	}

	return disjoint_vectors(left, right); // check if a cluster spans the lattice horizontally and return result
}

int is_spanning(int frame[], int grid_size){

	// Returns 1 if a cluster spans the lattice both horizontally and vertically, 0 otherwise

	int span_vertical = is_spanning_vertical(frame, grid_size);
	int span_horizontal = is_spanning_horizontal(frame, grid_size);

	if ((span_horizontal==1)&&(span_vertical==1)){
		return 1;
	}
	else{
		return 0;
	}
}

vector<int> spanning_cluster_coordinates(int frame[], int grid_size){

	// Returns a vector with coordinates of the spanning cluster (if present)

	vector<int> up; // empty vector to store upper border
	vector<int> down; // empty vector to store bottom border
	vector <int> left; // empty vector to store left border
	vector <int> right; // empty vector to store right border

	int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

	vector<int> spanning_cluster_coords;
	vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

	find_clusters_free_boundary(frame, labels, clusters, grid_size);
	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

	for (int i = 0; i<grid_size; i++){ // loop to populate up down left and right

		if (labels[i] != 0){
			up.push_back(labels[i]);
		}

		if (labels[grid_size*(grid_size-1)+i] != 0){
			down.push_back(labels[grid_size*(grid_size-1)+i]);
		}

		if (labels[grid_size*i] != 0){
			left.push_back(labels[grid_size*i]);
		}

		if (labels[grid_size*(i+1)-1] != 0){
			right.push_back(labels[grid_size*(i+1)-1]);
		}
	}

	if (is_spanning_horizontal(frame,grid_size)){ // conditional to check if the frame has a cluster that spans the lattice horizontally

		vector<int> spanning_cluster_label = common_elements(left, right); // label of the horizontally spanning cluster which must be present in both left and right

		for (int i=0; i < clusters[spanning_cluster_label[0]-1].coords.size(); i++){ // accumulate spanning cluster coordinates
			spanning_cluster_coords.push_back(clusters[spanning_cluster_label[0]-1].coords[i]);
		}

	}
	else if (is_spanning_vertical(frame,grid_size)){ // conditional to check if the frame has a cluster that spans the lattice vertically

		vector<int> spanning_cluster_label = common_elements(up, down); // label of the vertically spanning cluster which must be present in both up and down

		for (int i=0; i < clusters[spanning_cluster_label[0]-1].coords.size(); i++){ // accumulate spanning cluster coordinates
			spanning_cluster_coords.push_back(clusters[spanning_cluster_label[0]-1].coords[i]);
		}

	}
	else{ // the frame does not have a spanning cluster

		spanning_cluster_coords.push_back(-1);
		// Unique identification of absence of spanning cluster since coordinates must be non-negative integer

	}
	return spanning_cluster_coords;
}

vector<int> spanning_cluster_label_id(int frame[], int grid_size, int labels[], vector<cluster>& clusters)
{
  // Returns a vector with coordinates of the spanning cluster (if present), else 0 (which is not a valid label id)

  vector<int> up; // empty vector to store upper border
	vector<int> down; // empty vector to store bottom border
	vector <int> left; // empty vector to store left border
	vector <int> right; // empty vector to store right border

  //Only perform DFS search if one hasn't been performed on current frame before hand.

  if(clusters.empty())
  {
    //cluster data for current frame is NOT available. DFS needs be performed.

    find_clusters_free_boundary(frame, labels, clusters, grid_size);
  	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions
  }

  for (int i = 0; i<grid_size; i++)
  { // loop to populate up down left and right

		if (labels[i] != 0){
			up.push_back(labels[i]);
		}

		if (labels[grid_size*(grid_size-1)+i] != 0){
			down.push_back(labels[grid_size*(grid_size-1)+i]);
		}

		if (labels[grid_size*i] != 0){
			left.push_back(labels[grid_size*i]);
		}

		if (labels[grid_size*(i+1)-1] != 0){
			right.push_back(labels[grid_size*(i+1)-1]);
		}
	}

  std::vector<int> spanning_cluster_label = common_elements(left, right);
  // label of the horizontally spanning cluster which must be present in both left and right

  if(spanning_cluster_label.empty() == false)
  {
    // There is at least a horizontal spanning cluster.
    return spanning_cluster_label;
  }
  else if(common_elements(up, down).empty()==false)
  {
    // There is at least one spanning vertical cluster.
    spanning_cluster_label = common_elements(up, down);
    return spanning_cluster_label;
  }

  /*std::vector<int> h_spanning_cluster_label = common_elements(left, right);
  std::vector<int> v_spanning_cluster_label = common_elements(up, down);
  std::vector<int> spanning_cluster_label = common_elements(h_spanning_cluster_label, v_spanning_cluster_label);

  if(spanning_cluster_label.empty() == false)
  {
    // There is a 2D spanning cluster.
    return spanning_cluster_label;
  } */
  else
  {
    //No spanning cluster of any sort.
    spanning_cluster_label.push_back(-1);
    return spanning_cluster_label;
    // Returns -1 as ID if no spanning cluster exists.
  }


}

float how_many_red_sites(int frame[], int grid_size){

	// Gives a count of the red sites in the frame

	float red_sites = 0;

	vector<int> spanning_cluster_coords = spanning_cluster_coordinates(frame, grid_size); // get the coordinates of the spanning clusters

	if (spanning_cluster_coords[0]==-1){ // spanning cluster absent hence no red sites
		red_sites = 0;
	}
	else{ // spanning cluster present. Looking for red sites...
		for (int i=0; i<spanning_cluster_coords.size(); i++){ // traverse the coordinates of the spanning cluster

			frame[spanning_cluster_coords[i]] = 0; // flip the current coordinate to to check if spanning cluster is disrupted

			if (!is_spanning(frame,grid_size)){ // conditional is entered is spanning cluster is disrupted
				red_sites += 1; // the site is a red site because it disrupts the spanning cluster
			}

			frame[spanning_cluster_coords[i]] = 1; // restore the integrity of the frame

		}
	}
	return red_sites;
}

float find_average_cluster_size(int frame[], int grid_size){

	// Returns the average cluster size of the frame

	int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0
	vector <cluster> clusters;

	find_clusters_free_boundary(frame, labels,clusters, grid_size);
	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

	vector<float> cluster_sizes; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

	for (int i=0; i<clusters.size(); i++){
		cluster_sizes.push_back((float)clusters[i].coords.size()); // accumulate the cluster sizes
	}

	if (is_spanning(frame,grid_size)){ // conditional to check and remove the spanning cluster from cluster sizes

		vector<int> spanning_cluster_coords =  spanning_cluster_coordinates(frame, grid_size); // find coordinates of the spanning cluster
		int largest_cluster = spanning_cluster_coords.size(); // size of the spanning cluster
		cluster_sizes.erase(std::remove(cluster_sizes.begin(), cluster_sizes.end(), largest_cluster), cluster_sizes.end());
		// remove size of the spanning cluster from cluster sizes

	}

	float average_cluster_size = mean_of_vector(cluster_sizes,cluster_sizes.size()); // find the average of the cluster sizes
	return average_cluster_size;
}

float percolation_probability_np(int grid_size, float birth_probability, int number_of_census, int lag){

	// Finds an estimate of the percolation probability for NP

	int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 10000; // NP reaches steady state by 10000 for all values of p
	float percolation_probabilities[number_of_census];


	simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state

	for (int i = 0; i < number_of_census; ++i){

		simulate_np(frame,grid_size,birth_probability,lag); // simulate NP for an average of (lag) number of update per site

        if (is_spanning(frame,grid_size)){ // the frame has a spanning cluster

        	percolation_probabilities[i] = 1.0; // the lattice is percolating
        }
        else{ // the frame does not have a spanning cluster
        	percolation_probabilities[i] = 0.0; // the lattice is not percolating
        }
	}

	float percolation_probability = mean_of_array(percolation_probabilities,number_of_census);
	// fraction of frames from sample that had a spanning cluster

	return percolation_probability;
}

float percolation_probability_dp(int grid_size, float birth_probability, int number_of_census, int lag){

	// Finds an estimate of the percolation probability for DP

	int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site;
	float percolation_probabilities[number_of_census];


	if (birth_probability > 0.65){
		updates_per_site = 25000; // DP reaches steady state by 25000 for all values of p > 0.65
	}
	else{
		updates_per_site = 100000; // close to the critical point, dynamics are slow.
	}

	simulate_dp(frame,grid_size,birth_probability,updates_per_site); // simulate DP till it reaches steady state


    for (int i=0; i<number_of_census; i++) {

     	simulate_dp(frame,grid_size,birth_probability,lag); // simulate DP for an average of (lag) number of update per site

     	if (is_spanning(frame,grid_size)){ // the frame has a spanning cluster

        	percolation_probabilities[i] = 1.0; // the lattice is percolating
        }
        else{ // the frame does not have a spanning cluster
        	percolation_probabilities[i] = 0.0; // the lattice is not percolating
        }
     }

    float percolation_probability = mean_of_array(percolation_probabilities,number_of_census);
    // fraction of frames from sample that had a spanning cluster

	return percolation_probability;
}

void percolation_probabilities_np(int grid_size, float p_start, float p_end, int divisions, int number_of_census, int lag){

	//OpenMP array reduce idiom.
	//Source: https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
	//Credit: https://stackoverflow.com/users/2542702/z-boson
	//Used with modification

	// Find percolation probabilities for NP over a range of value of p by running parallel simulations

	vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

	float percolation_probabilities[divisions] = {0}; // initialize percolation probabilities

	#pragma omp parallel // The implementation below is to obtain an order of percolation probabilities that shadows the order of birth_probabilities
	{
	    float percolation_probabilities_private[divisions] = {0};

	    #pragma omp for
		for (int i=0; i < divisions; i++){

			int seed = std::random_device{}();
			rng.seed(seed);

			percolation_probabilities_private[i] = percolation_probability_np(grid_size, birth_probabilities[i], number_of_census, lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            percolation_probabilities[n] += percolation_probabilities_private[n];
	        }
	    }
	}

	for (int i=0; i< divisions; i++){
		cout << "p: " << setprecision(3) << birth_probabilities[i] << " P: " << setprecision(3) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}
}

void percolation_probabilities_dp(int grid_size, float p_start, float p_end, int divisions, int number_of_census, int lag){

	//OpenMP array reduce idiom.
	//Source: https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
	//Credit: https://stackoverflow.com/users/2542702/z-boson
	//Used with modification

	// Find percolation probabilities for DP over a range of value of p by running parallel simulations

	vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

	float percolation_probabilities[divisions] = {0}; // initialize percolation probabilities

	#pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
	{
	    float percolation_probabilities_private[divisions] = {0};

	    #pragma omp for
		for (int i=0; i < divisions; i++){

			int seed = std::random_device{}();
			rng.seed(seed);

			percolation_probabilities_private[i] = percolation_probability_dp(grid_size, birth_probabilities[i], number_of_census, lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            percolation_probabilities[n] += percolation_probabilities_private[n];
	        }
	    }
	}

	for (int i=0; i< divisions; i++){
		cout << "p: " << setprecision(3) << birth_probabilities[i] << " P: " << setprecision(3) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}
}

float average_cluster_size_np(int grid_size, float birth_probability, int number_of_census, int lag){

	int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 10000; // NP reaches steady state by 10000 for all values of p
	float average_cluster_size[number_of_census];


	simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state


    for (int i=0; i<number_of_census; i++) {

     	simulate_np(frame,grid_size,birth_probability,lag); // simulate NP for an average of (lag) number of update per site
     	average_cluster_size[i] = find_average_cluster_size(frame, grid_size); // find the average cluster size of the current frame

    }

    return mean_of_array(average_cluster_size,number_of_census);
}

float average_cluster_size_dp(int grid_size, float birth_probability, int number_of_census, int lag){

	int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site;


	if (birth_probability>0.65){
		updates_per_site = 25000; // DP reaches steady state by 25000 for all values of p > 0.65
	}
	else{
		updates_per_site = 100000; // close to the critical point, dynamics are slow.
	}

	float average_cluster_size[number_of_census];


	simulate_dp(frame,grid_size,birth_probability,updates_per_site);  // simulate NP till it reaches steady state


    for (int i=0; i<number_of_census; i++) {

     	simulate_dp(frame,grid_size,birth_probability,lag);  // simulate NP for an average of (lag) number of update per site
     	average_cluster_size[i] = find_average_cluster_size(frame, grid_size); // find the average cluster size of the current frame

    }

    return mean_of_array(average_cluster_size,number_of_census);
}

//---------------------------Correlation Function-----------------------------------------------//

zd_coordinates crosscorrelation_2D(int frame1[], int frame2[], int grid_size)
{
    float croscolmatr[2*grid_size -1][2*grid_size -1];

    //float rho1 = calculate_density(frame1, grid_size);
    //float rho2 = calculate_density(frame2, grid_size);
    f_coordinates rho_sig1 = calculate_SD(frame1, grid_size);
    f_coordinates rho_sig2 = calculate_SD(frame2, grid_size);

    int k_min, k_max, l_min, l_max;
    k_min = l_min = -(grid_size -1);
    k_max = l_max = grid_size -1;

    for(int k= k_min; k <= k_max; k++)
    {
      for(int l = l_min; l <= l_max; l++)
      {
        float sum =0.0;
        for(int m=0; m < grid_size; m++)
        {
          if( m - k >= grid_size || m - k < 0)
          {   continue;   }
          for(int n=0; n < grid_size; n++)
          {
            if( n - l >= grid_size || n - l < 0)
            {   continue;   }
            sum += (frame1[m*grid_size + n] - rho_sig1.x)*(frame2[(m - k)*grid_size + n -l] - rho_sig2.x);
          }
        }
        croscolmatr[k + grid_size -1][l + grid_size -1] =sum/(rho_sig1.y*rho_sig2.y*grid_size*grid_size);
        // Denomiator is the normalisation factor.
      }
    }
    zd_coordinates max ={-1000.0, -1.0, -1.0};
    for(int m=0; m < 2*grid_size-1; m++)
    {
      for(int n=0; n < 2*grid_size-1; n++)
      {
        if(croscolmatr[m][n] > max.x)
        {
          max.x =croscolmatr[m][n]; max.y = m -(grid_size -1); max.z = n - (grid_size -1);
        }
        //cout << setprecision(4) << croscolmatr[m][n] << "|";
      }
      //cout << endl;
    }
    stringstream msg;
    msg <<"\n Max Coeff:\t" <<setprecision(5) << max.x << " at: (" << max.y << "," << max.z << ")\n";
    cout << msg.str();
    return max;



}

void CrossCol_DP(int grid_size, float p, int number_of_census, int r_init, int lag)
{
  lag=1;
  stringstream peon, p_en, rini;

  peon << setprecision(3) << p;
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;

  ofstream outputACF; //Creating a file stream for storing the output in the form of  a CSV.

  outputACF.open("ACF/DP_L_"+ std::to_string(grid_size) + "_p_" + peon.str() + "_r_" + rini.str() + "_Cen_"+ std::to_string(number_of_census) + ".csv");
  // Creating CSV file in "ACF" sub-directory to store output data
  std::vector<std::vector<double>> output;
  //Will store ACF(t) data for every t <length for every division (trial) in order.

  //The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
  #pragma omp parallel
  {

      std::vector<std::vector<double>> kek;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

      #pragma omp for nowait schedule(static)
      for (int r=0; r< r_init; r++)
      {
        int seed = std::random_device{}();
        rng.seed(seed);
        int frame[grid_size*grid_size];
	      random_frame(frame, grid_size); // Initialize a random frame
	      int updates_per_site = 50000; // NP reaches steady state by 10000 for all values of p

        simulate_dp(frame, grid_size, p ,updates_per_site); // simulate NP till it reaches steady state

        int f1[grid_size*grid_size]; int f2[grid_size*grid_size];
        //copy(begin(frame), end(frame), begin(f1)); // Copying frame[] to f1[]

        for(int i=0; i<grid_size*grid_size; i++)
        {
          f1[i] = frame[i];
        }

        //std::vector<zd_coordinates> kek;

        for(int i = 0; i < number_of_census; i++)
        {

          for(int j=0; j<grid_size*grid_size; j++)
          {
            f2[j] = frame[j];
          }

          zd_coordinates kekistan = crosscorrelation_2D(f1, f2, grid_size);
          kek.push_back({i*grid_size*grid_size, kekistan.x, kekistan.y, kekistan.z});

          for(int n =0; n < grid_size*grid_size; n++)
            {
              dp_update(frame, grid_size, p); //Single update.
            }
            //copy(begin(frame), end(frame), begin(f2)); // Copying frame[] to f1[]
        }
      }
      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          output.insert(output.end(), kek.begin(), kek.end());
          // Inserting ACF data for each trial in order.
      }
  }
  /*std::vector<std::vector<double>> output;
  for(int i = 0; i < number_of_census; i++)
  {
    output.push_back({i, vec[i].x, vec[i].y, vec[i].z});
  } */
  cout << "The vector elements are: " << endl;
  cout << "Time step , Max Correl Coeff , k ,  l\n";
  for (int i = 0; i < output.size(); i++)
  {
    cout << setprecision(8) << output[i][0] << "  " << setprecision(6) << output[i][1] << "\t" << setprecision(5) << output[i][2] << "  " << setprecision(5) << output[i][3] << endl;
  }

  // Saving to aforementioned CSV

  outputACF << "Time step , Max Correl Coeff , k ,  l\n";
  for (int i = 0; i < output.size(); i++)
  {
    outputACF << setprecision(12) << output[i][0] << "," << setprecision(6) << output[i][1] << "," << setprecision(8) << output[i][2] << "," << setprecision(8) << output[i][3] << endl;
  }
  outputACF.close();

}
