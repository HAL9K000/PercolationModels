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

void solitary_droplet(int frame[], int grid_size)
{
  //A solitary droplet, lonesome tear crystallised in a lattice, spreading regret and loss.

  for (int i = 0; i < grid_size; ++i)
  {
        for (int j = 0; j < grid_size; ++j)
        {
            frame[i*grid_size + j] = 0;
        }
  }
  int numero_uno; //Stores index location of central position in frame

  if(grid_size%2 == 0)
  {
    //Grid size is even.
    numero_uno = (grid_size + 1)*grid_size/2;
  }
  else
  {
    //We can find an exact centre.
    numero_uno = (grid_size + 1)*(grid_size-1)/2;
  }
  frame[numero_uno] = 1;

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

int number_of_elem_of_array(int frame[],int grid_size)
{
  float sum = 0.0;

	for (int i =0; i<grid_size*grid_size; ++i)
  {
		sum += frame[i];
	}
  return sum;
}

void increase_stack_limit(int stack_size){

	// This becomes necessary for recursive depth first search when finding clusters.

	const rlim_t kStackSize = 1024 * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    cout << rl.rlim_max << endl;
    cout << rl.rlim_cur << endl;
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

	//https://www.faceprep.in/c/check-if-the-given-arrays-are-disjoint-in-c-c-java-and-python-faceprep/

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


void tp_update(int frame[], int grid_size, float birth_probability, float feedback_strength)
{

	/* Performs single asynchronous update for TP at a randomly selected site on the lattice

  Transitions take place as per the procedure delineated in LUBECK 2006.

  A ---> 0     with probability (1-q)(1-p)
  AA ---> 0A      with probability (1-q)(1-p) [Subcase of above]
  0A ---> 00      with probability q(1-p)
  0A ---> AA      with probability (1-q)(p)
  0AA ---> AAA      with probability q

  */

	int x = random_int(0, grid_size-1);
	int y = random_int(0, grid_size-1);

	coordinates site;
	site.x = x;
	site.y = y;

	if (frame[x*grid_size+y]==1)
  { // selected site is occupied

		coordinates neighbour = select_neighbor_of_site(site,grid_size);

		if (frame[neighbour.x*grid_size+neighbour.y]==1)
    { // neighbour of focal site is occupied. TCP like updates in this case

			double chance = random_real(0, 1);
      double another_chance = random_real(0, 1);

			if (chance < feedback_strength)
      {
				coordinates neighbour_of_pair = select_neighbor_of_pair(site,neighbour,grid_size);
				frame[neighbour_of_pair.x*grid_size+neighbour_of_pair.y] = 1; // birth at neighbour of pair with enhanced probability q
			}
			else if (another_chance < 1-birth_probability)
      {

				frame[x*grid_size+y]=0; // death at focal site with diminished death probability (1-q)(1-p)
        /** AA ---> 0A      (1-q)(1-p) */
			}
    }
    else
    { // neighbour of focal site is unoccupied. We have the plain old DP, but with a fresh twist.


			double chance = random_real(0, 1);
      double another_chance = random_real(0, 1);

			if (chance < 1 - feedback_strength)
      {
        if( another_chance < birth_probability)
        {
				      frame[neighbour.x*grid_size+neighbour.y] = 1; //birth at neighbour of occupied focal site with probability (1-q)*p
              /** 0A ---> AA      (1-q)(p) */
        }
        else
        {
          frame[x*grid_size+y]=0; //Death at focal site with probability (1-q)(1-p)
          /** A0 ---> 00      (1-q)(1-p)

          When taken together with AA ---> 0A transition above, we collectively have the transition:
           A ---> 0     with prob (1-q)(1-p)
           This is the probability of the focal site becomes empty, irrespective of the status of it's neighbouring site.
          */
        }
			}
			else if( another_chance < 1 - birth_probability) //Chance < q
      {
				  frame[x*grid_size+y]=0; //Death at focal site with probability q(1-p)
          /** 0A ---> 00      q(1-p) */
			}
		}

	}
}

/** UPDATE RULES PROVIDED BY AYAN BEFORE:


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

*/

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

float equilibrium_density_dp(int grid_size, float birth_probability, int r_init, int number_of_census, int lag, int updates_per_site, int collect_frames){

	// Simulates DP with specified parameters and prints the mean and standard deviation of vegetation cover to the terminal. The frames from
	// which vegetation cover is calculated can also be captured. Finally it returns the mean vegetation cover.
  long limit = r_init*number_of_census;
	float densities[limit];
	int frame[grid_size*grid_size];

  for(int r=0; r< r_init; r++)
  {

    random_frame(frame, grid_size);

    simulate_dp(frame,grid_size,birth_probability,updates_per_site);

    for (int i=0; i<number_of_census; i++)
    {

    	simulate_dp(frame,grid_size,birth_probability,lag);
    	densities[r*r_init + i] = calculate_density(frame,grid_size);

    	if (collect_frames ==1){ // conditional for collecting frames

    		ofstream outdata;

			stringstream p, rini, census, Lag;
			p << setprecision(3) << birth_probability;
      rini << r;
			census << i;
			Lag << lag;

			outdata.open("dump/dp_frame_"+std::to_string(grid_size)+"_"+p.str()+"_"+Lag.str()+'_'+census.str()+'_'+rini.str()+".txt");
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

    float mean_density = mean_of_array(densities,limit);
    float standard_deviation_density = standard_deviation_of_array(densities,limit);

    cout << endl;
    stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
    message << "p: " << setprecision(4) << birth_probability << " Mean: " << setprecision(5) << mean_density << " Standard Deviation: " << setprecision(7) << standard_deviation_density << endl;
    cout << message.str();
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

bool sort_s_cluster(coordinates a, coordinates b)
{ return a.y < b.y; // Returning smaller s value
}

vector<zd_coordinates> gen_ns_data(vector<coordinates>& cluster_details, double p)
{
  /*This function returns raw collected data (refer to "cluster_details" variable in the method below)
  and sorts and returns that data in the format of (p, s, n_s(p)) */

  vector<zd_coordinates> ns_data;
  // Vector to return sorted data in the format of (p, s, n_s(p)).

  //Sorting the cluster_details variable in ascending order of s.

  sort(cluster_details.begin(), cluster_details.end(), sort_s_cluster);

  //Finding max value of s (that is NOT spanning).

  int s_max = cluster_details[cluster_details.size() -1].y;

  int n_s=0; int s=1;
  for(int i=0; i<cluster_details.size();i++)
  {

    if(cluster_details[i].y == s)
    {
      n_s+=1;
    }
    else if(cluster_details[i].y > s)
    {
      // Size has increased. Store n_s data.
      zd_coordinates temp; //Creating temporary dummy variable to store present data.
      temp.x = p; temp.y = static_cast<double>(s); temp.z = static_cast<double>(n_s);
      ns_data.push_back(temp);

      n_s=1; //Resetting counter to 1.
      s= cluster_details[i].y;
    }
  }
  //For largest s (that is not spanning cluster).

  zd_coordinates temp; //Creating temporary dummy variable to store present data.
  temp.x = p; temp.y = static_cast<double>(s); temp.z = static_cast<double>(n_s);
  ns_data.push_back(temp);

  return ns_data;

}

void tau_patch_size_find_np(int grid_size, vector<zd_coordinates>& tau_data, double p, int r_init, int number_of_census, int lag)
{
  // For a given value of p, generate "number_of_census" updates of static frames, each lag distance apart.

  int frame[grid_size*grid_size];
  //int r_init = 5; //Five random frames will be generated one after another.
	//random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 8000; // NP reaches steady state by 10000 for all values of p
  //DP reaches steady state by 50000 for all values of p.
  long limit = r_init*number_of_census;     //Stores the length of percolation_probabilities array necessary.
	double percolation_probabilities[limit];

	//simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state

  //simulate_dp(frame,grid_size,birth_probability,updates_per_site);

  for(int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Assign a random frame
    simulate_dp(frame,grid_size,p,updates_per_site); // simulate NP till it reaches steady state

    //simulate_dp(frame, grid_size, p, updates_per_site);

    /*stringstream msg2;
    msg2 << "Nuka\t" << i << "\t" << lag << "\t" << number_of_census << "\n";
    cout << msg2.str(); */

    // random_frame_of_density(p, frame, grid_size);
    // Static percolation.

    for (int j = 0; j < number_of_census; ++j)
    {
      float k=0; // Stores the percolation strength, if applicable.

      //simulate_np(frame,grid_size,p,lag); // simulate NP for an average of (lag) number of update per site

      simulate_dp(frame,grid_size,p,lag);

      int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

      vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

      find_clusters_free_boundary(frame, labels, clusters, grid_size);
      // Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

      vector<coordinates> cluster_details;
      // Stores label id in x attribute, cluster size associated with label in 2nd coordinate.

      vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
      //Returns labels of spanning cluster(s) if present, -1 otherwise.

      /*stringstream msg5;
      msg5 << "Nunia\t" << j << "\t" << lag << "\n";
      cout << msg5.str(); */

      //Purpose of following nested loops is to populate cluster_details following data structure noted above,
      //minus the details of the spanning cluster(s).
    	for (int a=0; a<clusters.size(); a++)
      {
          // There may exist a spanning cluster.
          int flag=0; //Flag variable used to detect match with spanning cluster(s).
          if(spanning_cluster_labels[0]> -1)
          { //There exists a spanning cluster.
            for (int b=0; b<spanning_cluster_labels.size(); b++)
            {
              //Iterating over spanning cluster(s) indices.
              if(spanning_cluster_labels[b] == clusters[a].label)
              { flag=1;}
            }
          }
          if(flag == 0)
          {
            //No match of given cluster with spanning cluster.
            coordinates temp; //Temporary x, y variable declared.
            temp.x = clusters[a].label;
            temp.y = clusters[a].coords.size();
            cluster_details.push_back(temp);
          }
      }

      double bok= i*number_of_census + j + 1; //Stores current trial number

      if (static_cast<int>(bok)% 40 == 1)
      {
        //Every 25th term, the following  will be outputted to help with debugging.
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        if(spanning_cluster_labels[0]> -1)
        { message << "1st spanning clust size:\t" << clusters[spanning_cluster_labels[0]-1].coords.size() << "\t Num of spn clust:\t" << spanning_cluster_labels.size() << "\n";  }
        else
        { message << "No Spanning Cluster Present. \n"; }
        cout << message.str();

      }
      zd_coordinates bokaro = {bok, p, p};
      tau_data.push_back(bokaro);
      //All "headers" seperating different trials have (tr #, p, p) values in place of (p, s, n_s)

      vector<zd_coordinates> trail_data = gen_ns_data(cluster_details, p);
      //Sorting cluster_details in ascending order of s, and outputting (p,s, n_s(p,s)) data from it.

      tau_data.insert(tau_data.end(), trail_data.begin(), trail_data.end());

      if (spanning_cluster_labels[0] != -1)
      {
        //Spanning cluster(s) present
        for (int r =0 ; r < spanning_cluster_labels.size(); r++)
        {
          k += (float) clusters[spanning_cluster_labels[r]-1].coords.size();
        }
      }
      bokaro.x = p; bokaro.y = -10.0; bokaro.z = k;
      //Finally, storing details of spanning cluster at the very end as (p, -10, size_spn_cls)
      // Here s= -10 is the marker for spanning cluster. If absent, k =0
      tau_data.push_back(bokaro);

      trail_data.clear(); clusters.clear(); spanning_cluster_labels.clear(); cluster_details.clear();


    } //End of census for loop.

  } //End of outer loop.

  stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
  message << "Length of tau_data for all trials at this given p:\t" << tau_data.size() << "\n";
  cout << message.str();


}

void iter_patch_sizes_np(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag)
{
  /* This method will iteratively generate static cluster size (n_s(p)) distributions for all specified values
  of p. The ultimate goal is to estimate the value of tau. */

  //r_init =5;
  // Number of random trials to be initiated for each p value

  ofstream outputtau;
  // Creating a file instance called output to store output data (for Tau exponent) as CSV.

  std::vector<double> birth_probabilities = linspace(p_start, p_end, divisions);

  stringstream peon, p_en, rini;

  peon << setprecision(3) << p_start;
  p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;
  outputtau.open("CrtExp/Tau_SP_L_"+ std::to_string(grid_size) + "_p1_" + peon.str() + "_p2_" + p_en.str() + "_r_" + rini.str() + "_Cen_"+ std::to_string(number_of_census) + ".csv");
  // Creating CSV file in "ACF" sub-directory to store output data

  std::vector<zd_coordinates> vec;
  //Will store ACF(t) data for every t <length for every division (trial) in order.
  // First col stores p value, second the s value and the third the n_s value.

  //The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
  #pragma omp parallel
  {
      std::vector<zd_coordinates> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

      #pragma omp for nowait schedule(static)
      for (int i=0; i < divisions; i++)
      {

        int seed = std::random_device{}();
        rng.seed(seed);

        std::vector<zd_coordinates> tau_data;

        tau_patch_size_find_np(grid_size, tau_data, birth_probabilities[i], r_init, number_of_census, lag);
        //Finds and stores ns(p) data in tau_data variable.

        vec_private.insert(vec_private.end(), tau_data.begin(), tau_data.end());

      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting ACF data for each trial in order.
      }
  }
  vector <vector<double>> output;
  //Creating 2D vector.

  double trialno =1;
  for(int i=0; i <vec.size(); i++)
  {
    if(vec[i].x >= 1 && vec[i].z < 1)
    {
      //The above condition is only satisfied if vec[i] is a "header" row denoting separation b/w trials.
      trialno = vec[i].x; //Setting to current trial number.
    }
    else
    {
      //Otherwise proceed to fill entry.
      output.push_back({trialno, vec[i].x, vec[i].y, vec[i].z});
    }
  }

  //Printing out obtained results.
  cout << "The vector elements are: "<< endl;
  cout << "# Tr No , p , s ,  ns(p)\n";
  for (int i = 0; i < output.size(); i++)
  {
    cout << setprecision(3) << output[i][0] << "  " << setprecision(6) << output[i][1] << "  " << setprecision(8) << output[i][2] << "  " << setprecision(8) << output[i][3] << endl;
  }

  // Saving to aforementioned CSV

  outputtau << "# Tr No , p , s , ns(p)\n";
  for (int i = 0; i < output.size(); i++)
  {
    outputtau << setprecision(3) << output[i][0] << "," << setprecision(6) << output[i][1] << "," << setprecision(8) << output[i][2] << "," << setprecision(8) << output[i][3] << endl;
  }
  outputtau.close();

}

void tau_patch_size_find_tcp(int grid_size, vector<zd_coordinates>& tau_data, double p, double q, int r_init, int number_of_census, int lag)
{
  // For a given value of p, generate "number_of_census" updates of static frames, each lag distance apart.

  int frame[grid_size*grid_size];
	int updates_per_site = 100000; // NP reaches steady state by 100,000 for all values of p
  //DP reaches steady state by 50000 for all values of p.
  long limit = r_init*number_of_census;     //Stores the length of percolation_probabilities array necessary.
	double percolation_probabilities[limit];

	//simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state

  //simulate_dp(frame,grid_size,birth_probability,updates_per_site);

  for(int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Assign a random frame
    simulate_tp(frame,grid_size,p, q, updates_per_site); // simulate TCP till it reaches steady state

    for (int j = 0; j < number_of_census; ++j)
    {
      float k=0; // Stores the percolation strength, if applicable.

      simulate_tp(frame,grid_size,p,q, lag);

      int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

      vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

      find_clusters_free_boundary(frame, labels, clusters, grid_size);
      // Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

      vector<coordinates> cluster_details;
      // Stores label id in x attribute, cluster size associated with label in 2nd coordinate.

      vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
      //Returns labels of spanning cluster(s) if present, -1 otherwise.

      //Purpose of following nested loops is to populate cluster_details following data structure noted above,
      //minus the details of the spanning cluster(s).
    	for (int a=0; a<clusters.size(); a++)
      {
          // There may exist a spanning cluster.
          int flag=0; //Flag variable used to detect match with spanning cluster(s).
          if(spanning_cluster_labels[0]> -1)
          { //There exists a spanning cluster.
            for (int b=0; b<spanning_cluster_labels.size(); b++)
            {
              //Iterating over spanning cluster(s) indices.
              if(spanning_cluster_labels[b] == clusters[a].label)
              { flag=1;}
            }
          }
          if(flag == 0)
          {
            //No match of given cluster with spanning cluster.
            coordinates temp; //Temporary x, y variable declared.
            temp.x = clusters[a].label;
            temp.y = clusters[a].coords.size();
            cluster_details.push_back(temp);
          }
      }

      double bok= i*number_of_census + j + 1; //Stores current trial number

      if (static_cast<int>(bok)% 40 == 1)
      {
        //Every 25th term, the following  will be outputted to help with debugging.
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        if(spanning_cluster_labels[0]> -1)
        { message << "1st spanning clust size:\t" << clusters[spanning_cluster_labels[0]-1].coords.size() << "\t Num of spn clust:\t" << spanning_cluster_labels.size() << "\n";  }
        else
        { message << "No Spanning Cluster Present. \n"; }
        cout << message.str();

      }
      zd_coordinates bokaro = {bok, p, q};
      tau_data.push_back(bokaro);
      //All "headers" seperating different trials have (tr #, p, q) values in place of (p, s, n_s)

      vector<zd_coordinates> trail_data = gen_ns_data(cluster_details, p);
      //Sorting cluster_details in ascending order of s, and outputting (p,s, n_s(p,s)) data from it.

      tau_data.insert(tau_data.end(), trail_data.begin(), trail_data.end());

      if (spanning_cluster_labels[0] != -1)
      {
        //Spanning cluster(s) present
        for (int r =0 ; r < spanning_cluster_labels.size(); r++)
        {
          k += (float) clusters[spanning_cluster_labels[r]-1].coords.size();
        }
      }
      bokaro.x = p; bokaro.y = -10.0; bokaro.z = k;
      //Finally, storing details of spanning cluster at the very end as (p, -10, size_spn_cls)
      // Here s= -10 is the marker for spanning cluster. If absent, k =0
      tau_data.push_back(bokaro);

      trail_data.clear(); clusters.clear(); spanning_cluster_labels.clear(); cluster_details.clear();


    } //End of census for loop.

  } //End of outer loop.

  stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
  message << "Length of tau_data for all trials at this given p and q:\t" << tau_data.size() << "\n";
  cout << message.str();


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

void find_equilibrium_single_shot_transformations_tp(int grid_size, int divisions, int time_to_equilibrium,
	int how_many, vector<transformation>& transformations, vector<f_coordinates>& pq)
{

	// Simulates TP for the specified parameters and collects the specified number of transformations from the steady state

  std::vector<vector<double>> vec; //Stores final output.

  ofstream output_dels;
  // Creating a file instance called output to store output data as CSV.

  stringstream peon, peoff, div, quint, quaffle, gangsta, tallaq;

  peon << setprecision(4) << pq[0].x;
  quint << setprecision(3) << pq[0].y;
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  div << divisions;
  peoff << pq[divisions-1].x;
  quaffle << pq[divisions-1].y;
  gangsta << grid_size; tallaq << how_many;

  output_dels.open("dump/TCP_delS_G_" + gangsta.str() + "_N_" + tallaq.str() + "_p1_" +  peon.str() + "_q1_" + quint.str() + "_Div_" + div.str() + "_p2_"+ peoff.str() + "_q2_"+ quaffle.str() + ".csv");


  #pragma omp parallel
  {
    vector<vector<double>> collateral;

    int previous_frame[grid_size*grid_size];
  	int previous_frame_labels[grid_size*grid_size];
  	vector<cluster> previous_frame_clusters;

  	int current_frame[grid_size*grid_size];
  	int current_frame_labels[grid_size*grid_size];
  	vector<cluster> current_frame_clusters;

    #pragma omp for nowait schedule(static)
    for(int i= 0; i < divisions; i++)
    {


      int seed = std::random_device{}();
  	  rng.seed(seed);

      random_frame(current_frame, grid_size); // initialize the current frame
    	zeros(current_frame_labels,grid_size); // initialize labels for the current frame

    	simulate_tp(current_frame, grid_size, pq[i].x, pq[i].y, time_to_equilibrium); // simulate TP to reach steady state

      find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size); // find clusters from the current frame

      vector<transformation> trans; //To be used in collecting |del s| data.
      TCP_scattershot(trans, previous_frame, previous_frame_labels, previous_frame_clusters,
  		current_frame, current_frame_labels, current_frame_clusters, grid_size, how_many, pq[i].x, pq[i].y);
      // TCP_scattershot() is invaluable to completelting these simulations.

      for(int j=0; j<trans.size(); j++)
      {
        //Formatting trans data into a more appropriate format.
        collateral.push_back({pq[i].x, pq[i].y, j, trans[j].before, trans[j].after});
        //Filled in as p, q, #, s, s + del(s).
      }

    }
    #pragma omp for schedule(static) ordered
    for(int i=0; i< omp_get_num_threads(); i++)
    {
      #pragma omp ordered
        vec.insert(vec.end(), collateral.begin(), collateral.end());
        // Inserting critical exponent data for each grid size in order.
        stringstream message3;
        message3 << "Is this happening?\n";
        cout << message3.str();
    }

  }



	/* find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size); // find clusters from the current frame

	while (transformations.size() < how_many){ // loop till we have the specified number of transformations

		for (int i=0; i<grid_size; i++){
			for (int j=0; j<grid_size; j++){ // scan the current frame and duplicate it into the previous frame
				previous_frame[grid_size*i+j] = current_frame[grid_size*i+j];
				previous_frame_labels[grid_size*i+j] = current_frame_labels[grid_size*i+j];
			}
		}

		previous_frame_clusters = current_frame_clusters; // duplicate the clusters too

		tp_update(current_frame, grid_size, birth_probability, feedback_strength); // make a single update to the current frame

		current_frame_clusters.clear(); // wash the stale clusters off the current frame (since it has been updated now)
		zeros(current_frame_labels,grid_size); // wash the labels too

		find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size);
		// find the clusters and labels of the current frame (after update) again

		find_transformations_single_shot(transformations,
			previous_frame, previous_frame_labels, previous_frame_clusters,
			current_frame, current_frame_labels, current_frame_clusters, grid_size); // find the transformation and accumulate into the transformations vector
	} */

  cout << "| p , q, # Tr No , s ,  s + del(s) |\n";
  output_dels << " p , q, # Tr No, s,  s + del(s) \n";

  for (int i = 0; i < vec.size(); i++)
  {
    if(i%3333 ==1)
    {
    cout << setprecision(5) << vec[i][0] << "  " << setprecision(5) << vec[i][1] << "  " << setprecision(12) << vec[i][2] << "  " << setprecision(10) << vec[i][3] << "  " << setprecision(10) << vec[i][4] << endl;
    }
  }
  // Saving to aforementioned CSV

  for (int i = 0; i < vec.size(); i++)
  {
    output_dels << setprecision(8) << vec[i][0] << "," << setprecision(5) << vec[i][1] << "," << setprecision(12) << vec[i][2] << "," << setprecision(11) << vec[i][3] << "," << setprecision(11) << vec[i][4]  << endl;
  }
  output_dels.close();


}

void TCP_scattershot(vector <transformation>& trans, int previous_frame[], int previous_frame_labels[], vector<cluster>& previous_frame_clusters,
int current_frame[], int current_frame_labels[], vector<cluster>& current_frame_clusters, int grid_size, int how_many, double p, double q)
{

  auto start = high_resolution_clock::now();

	while (trans.size() < how_many){ // loop till we have the specified number of transformations

		for (int i=0; i<grid_size; i++){
			for (int j=0; j<grid_size; j++){
        // scan the current frame and duplicate it into the previous frame
				previous_frame[grid_size*i+j] = current_frame[grid_size*i+j];
				previous_frame_labels[grid_size*i+j] = current_frame_labels[grid_size*i+j];
			}
		}

		previous_frame_clusters = current_frame_clusters; // duplicate the clusters too

		tp_update(current_frame, grid_size, p, q); // make a single update to the current frame

		current_frame_clusters.clear(); // wash the stale clusters off the current frame (since it has been updated now)
		zeros(current_frame_labels,grid_size); // wash the labels too

		find_clusters(current_frame,current_frame_labels,current_frame_clusters,grid_size);
		// find the clusters and labels of the current frame (after update) again

		find_transformations_single_shot(trans,
			previous_frame, previous_frame_labels, previous_frame_clusters,
			current_frame, current_frame_labels, current_frame_clusters, grid_size); // find the transformation and accumulate into the transformations vector
	}

  auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

  stringstream msg; msg << "For (p,q) pair values as:\t(" << setprecision(6) << p << setprecision(5) << q << ")\t Transformation Finding Time: " << duration.count() << " seconds" << endl;

  cout << msg.str();

  return; //Exit method once desired number of |del s| values have been collected.
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

	if (is_spanning_horizontal(frame,grid_size))
  { // conditional to check if the frame has a cluster that spans the lattice horizontally

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

	// Returns the average cluster size of the frame (minus size of percolating cluster if any).

	int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0
	vector <cluster> clusters;

	find_clusters_free_boundary(frame, labels,clusters, grid_size);
	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

	vector<float> cluster_sizes; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

  vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
  //Stores list of spanning cluster IDs, -1 otherwise (if no spanning clusterfuck present).

  // The following loop will populate "cluster_sizes" with all cluster sizes other than the spanning cluster.
	for (int i=0; i<clusters.size(); i++)
  {

    // There may exist a spanning cluster.
    int flag=0; //Flag variable used to detect match with spanning cluster(s).
    if(spanning_cluster_labels[0]> -1)
    { //There exists a spanning cluster.
      for (int j=0; j<spanning_cluster_labels.size(); j++)
      {
        //Iterating over spanning cluster(s) indices.
        if(spanning_cluster_labels[j] == clusters[i].label)
        { flag=1;}
      }
    }
    if(flag == 0)
    {
      //No match of given cluster with spanning cluster.
      cluster_sizes.push_back((float)clusters[i].coords.size()); // accumulate the cluster sizes
    }
	}
  /* for (int i=0; i<cluster_sizes.size(); i++)
  {
    stringstream msg3;
    msg3 << setprecision(4) << cluster_sizes[i] << "\t";
    cout << msg3.str();
    if (i == cluster_sizes.size() -1)
    {
      stringstream msg4;
      msg4 << "\n";
      cout << msg4.str();
    }
  } */

	float average_cluster_size = mean_of_vector(cluster_sizes,cluster_sizes.size()); // find the average of the cluster sizes

  /*stringstream msg5;
  msg5 << "Average Cluster Size For This Census: " << setprecision(6) << average_cluster_size << "\n";
  cout << msg5.str();*/
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

	// Find percolation probabilities for NP over a range of value of p by running parallel simulations

	vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

	float percolation_probabilities[divisions] = {0}; // initialize percolation probabilities

	#pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
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

void average_cluster_size_np(int grid_size, vector<zd_coordinates> &avg_clust_size, float birth_probability,int r_init, int number_of_census, int lag){

	int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 10000; // NP reaches steady state by 10000 for all values of p
  int limit = r_init*number_of_census;
	float average_cluster_size[limit];

  for (int i=0; i< r_init; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Initialize a random frame
    //simulate_np(frame,grid_size,birth_probability,updates_per_site);
    // simulate NP till it reaches steady state
    random_frame_of_density(birth_probability, frame, grid_size);
    for (int j=0; j<number_of_census; j++)
    {

     	//simulate_np(frame,grid_size,birth_probability,lag); // simulate NP for an average of (lag) number of update per site
     	average_cluster_size[i*number_of_census+ j] = find_average_cluster_size(frame, grid_size); // find the average cluster size of the current frame

    }
  }
  float k = mean_of_array(average_cluster_size,limit);
  stringstream msg;
  msg << "For p:  " << setprecision(4) << birth_probability << " S(p):  " << k << "\n";
  cout << msg.str();

  /*for (int i=0; i< limit; i++)
  {
    stringstream msg2;
    msg2 << setprecision(4) << average_cluster_size[i] << "\t";
    cout << msg2.str();
    if (i == limit -1)
    {
      stringstream msg4;
      msg4 << "\n";
      cout << msg4.str();
    }
  }*/

  for (int i=0; i< limit; i++)
  {
    zd_coordinates temp; //Creating a temporary variable.
    int trl_no = i % number_of_census + 1; int r_no = int(i/number_of_census);
    temp.x = r_no*number_of_census + trl_no;
    temp.y = birth_probability;
    temp.z = average_cluster_size[i];
    avg_clust_size.push_back(temp);

  }

  //return avg_clust_size;

  //return mean_of_array(average_cluster_size,limit);
}

float average_cluster_size_dp(int grid_size, float birth_probability, int number_of_census, int lag){

	int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site;


	if (birth_probability>0.65){
		updates_per_site = 25000; // DP reaches steady state by 25000 for all values of p > 0.65
	}
	else
  {
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



// ----------------- Methods unique to theoretical calculation of Percolation Point For NP Model---------------------------------


int size_spanning_vertical(int frame[], int grid_size)
{
	// Returns size (i.e. number of nodes) of all the clusters spans the lattice vertically, 0 otherwise.

  vector<int> up;  // empty vector to store upper border
  vector<int> down;  // empty vector to store bottom border
  vector<int> commonality; //empty vector to store the common cluster labels present on both the upper and bottom margins

  int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

  int spanning_cluster_size=0;

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

  commonality = common_elements(up, down); // Returns a unique list of the common cluster labels on both the upper and lower borders.

  if (commonality.size() > 0)
  {
    //Vertically spanning clusters exist.

    for (int i=0; i < commonality.size(); i++)
    {
      // Looping through common cluster labels.

      spanning_cluster_size += clusters[commonality[i]-1].coords.size();

      // Storing values of number of sites belonging to spanning cluster.

    }

    return spanning_cluster_size;

  }
  else
  {
    return 0;
  }
}

int size_spanning_horizontal(int frame[], int grid_size)
{
	// Returns size (i.e. number of nodes) of all the clusters spans the lattice horizontally, 0 otherwise.

  vector<int> left;  // empty vector to store left border
  vector<int> right;  // empty vector to store right border
  vector<int> commonality; //empty vector to store the common cluster labels present on both the upper and bottom margins

  int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

  int spanning_cluster_size=0;

  vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

  find_clusters_free_boundary(frame, labels, clusters, grid_size);
  // Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions


  for (int i = 0; i<grid_size; i++){ // loop to populate up and down

    if (labels[grid_size*i] != 0){
      left.push_back(labels[grid_size*i]);
    }

    if (labels[grid_size*(i) + grid_size-1] != 0){
      right.push_back(labels[grid_size*(i) + grid_size-1]);
    }
  }

  commonality = common_elements(left, right); // Returns a unique list of the common cluster labels on both the upper and lower borders.

  if (commonality.size() > 0)
  {
    //Vertically spanning clusters exist.

    for (int i=0; i < commonality.size(); i++)
    {
      // Looping through common cluster labels.

      spanning_cluster_size += clusters[commonality[i]-1].coords.size();

      // Storing values of number of sites belonging to spanning cluster.

    }

    return spanning_cluster_size;

  }
  else
  {
    return 0;
  }
}

int size_spanning_2D(int frame[], int grid_size)
{
	// Returns size (i.e. number of nodes) of all the clusters that span the lattice both horizontally & Vertically,
  // 0 otherwise.

  vector<int> left;  // empty vector to store left border
  vector<int> right;  // empty vector to store right border
  vector<int> up;  // empty vector to store upper border
  vector<int> down;  // empty vector to store bottom border

  vector<int> commonality; //empty vector to store the common cluster labels present on both the upper and bottom margins

  int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

  int spanning_cluster_size=0;

  vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

  find_clusters_free_boundary(frame, labels, clusters, grid_size);
  // Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions


  for (int i = 0; i<grid_size; i++){ // loop to populate up, down, left & right.

    if (labels[grid_size*i] != 0){
      left.push_back(labels[grid_size*i]);
    }

    if (labels[grid_size*(i) + grid_size-1] != 0){
      right.push_back(labels[grid_size*(i) + grid_size-1]);
    }

    if (labels[i] != 0){
      up.push_back(labels[i]);
    }

    if (labels[grid_size*(grid_size-1)+i] != 0){
      down.push_back(labels[grid_size*(grid_size-1)+i]);
    }

  }

  vector<int> common1= common_elements(up, down);
  vector<int> common2= common_elements(right, common1);

  commonality = common_elements(left, common2);
  // Returns a unique list of the common cluster labels on both the upper and lower borders.

  if (commonality.size() > 0)
  {
    //2D spanning clusters exist.

    for (int i=0; i < commonality.size(); i++)
    {
      // Looping through common cluster labels.

      spanning_cluster_size += clusters[commonality[i]-1].coords.size();

      // Storing values of number of sites belonging to the 2D spanning cluster.

    }

    return spanning_cluster_size;

  }
  else
  {
    return 0;
  }
}

float theoretical_percolation_probability_np(int grid_size, float birth_probability, int r_init, int number_of_census, int lag)
{
  //Finds a Theoretical Estimate Of Percolation Probability.

  int frame[grid_size*grid_size];
	//random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 10000; // NP reaches steady state by 10000 for all values of p
  long limit = r_init*number_of_census;     //Stores the length of percolation_probabilities array necessary.
	float percolation_probabilities[limit];


	//simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state

  for(int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Assign a random frame
    //simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state
    random_frame_of_density(birth_probability, frame, grid_size); // Assign a random frame of density p.

    for (int j = 0; j < number_of_census; ++j)
    {
      //simulate_np(frame,grid_size,birth_probability,lag); // simulate NP for an average of (lag) number of update per site

      float k= (float) size_spanning_2D(frame, grid_size);
      //Stores the number of vertices that belong to 2D spanning tree.

      percolation_probabilities[i*number_of_census + j] = k/(grid_size*grid_size);
    }

  }

  float percolation_probability = mean_of_array(percolation_probabilities,limit);
	// Theoretical Percolation Probability averaged over censuses

	return percolation_probability;
}

void theoretical_percolation_probabilities_np(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag){

	// Find "THEORETICAL" percolation probabilities for NP over a range of value of p by running parallel simulations

	vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

	float percolation_probabilities[divisions] = {0}; // initialize percolation probabilities

  ofstream output;
  // Creating a file instance called output to store output data as CSV.
  ofstream outputbeta;

  stringstream p_st, p_en, rini, lagger;

  p_st << setprecision(3) << p_start;
  p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;
  lagger << lag;
  output.open("Theoretical_Percol/SP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + "_Lag_" + lagger.str() + ".csv");
  // Creating CSV file in "Theoretical_Percol" sub-directory to store output data.

  outputbeta.open("CrtExp/Beta_SP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");

	#pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
	{
	  float percolation_probabilities_private[divisions] = {0};

	  #pragma omp for
		for (int i=0; i < divisions; i++){

			int seed = std::random_device{}();
			rng.seed(seed);

			percolation_probabilities_private[i] = theoretical_percolation_probability_np(grid_size, birth_probabilities[i], r_init, number_of_census, lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            percolation_probabilities[n] += percolation_probabilities_private[n];
	        }
	    }
	}

  double p_c = 0.592746;
	//Percolation point for 2D Grid Network.
	std::vector<double> lnbase; 	// Stores ln(p - pc) values for various values of p.

	// Storing ln| p - p_c | for p ----> p_c+
	for(int i=0; i<divisions; i++)
	{
		double x =0.0;
		if (birth_probabilities[i] < p_c)
		{
			x= 0.0; //Only values greater than p_c are stored.
		}
		else
		{
			x= log(birth_probabilities[i] - p_c);
		}

		lnbase.push_back(x);
	}

	for (int i=0; i< divisions; i++){
		cout << "p " << setprecision(7) << birth_probabilities[i] << " Strength Of Percolation (P_infty) (Theoretical) " << setprecision(4) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}

  // Writing results to CSV File.

  output << "# Birth Probability (p) , Strength Of Percolation (P_infty) (Theoretical)\n";
  for (int i=0; i< divisions; i++){
		output << setprecision(8) << birth_probabilities[i] << "," << setprecision(7) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}
  output.close();


  outputbeta << "# Birth Probability (p) , ln|p - p_c|, Strength Of Percolation (P_infty) (Theoretical)\n";
  for (int i=0; i< divisions; i++)
  {
      if (lnbase[i] != 0)
      {
        // In other words, for all p > p_c
        outputbeta << setprecision(7) << birth_probabilities[i] << ","<< setprecision(10) << lnbase[i] << "," << setprecision(10) << percolation_probabilities[i] << endl;
    		/* Prints parameter value, ln|p - p_c|. and percolation probabilities to the terminal.
        The ordering is the same as the order of p in birth_probabilities */
      }

	}
  outputbeta.close();

}

float theoretical_percolation_probability_dp(int grid_size, float birth_probability, int r_init, int number_of_census, int lag)
{
  //Finds a Theoretical Estimate Of Percolation Probability For DP.

  int frame[grid_size*grid_size];
	//random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site =30000; // DP reaches steady state by 30000 for all values of p
  long limit = r_init*number_of_census;     //Stores the length of percolation_probabilities array necessary.
	float percolation_probabilities[limit];


	//simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state

  for(int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Assign a random frame
    simulate_dp(frame,grid_size,birth_probability,updates_per_site); // simulate DP till it reaches steady state
    //random_frame_of_density(birth_probability, frame, grid_size); // Assign a random frame of density p.

    for (int j = 0; j < number_of_census; ++j)
    {
      simulate_dp(frame,grid_size,birth_probability,lag); // simulate DP for an average of (lag) number of update per site

      float k= (float) size_spanning_2D(frame, grid_size);
      //Stores the number of vertices that belong to 2D spanning tree.

      percolation_probabilities[i*number_of_census + j] = k/(grid_size*grid_size);
    }

  }

  float percolation_probability = mean_of_array(percolation_probabilities,limit);
	// Theoretical Percolation Probability averaged over censuses

	return percolation_probability;
}

void theoretical_percolation_probabilities_dp(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag){

	// Find "THEORETICAL" percolation probabilities for NP over a range of value of p by running parallel DP simulations

	vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

	float percolation_probabilities[divisions] = {0}; // initialize percolation probabilities

  ofstream output;
  // Creating a file instance called output to store output data as CSV.
  ofstream outputbeta;

  stringstream p_st, p_en, rini, lagger;

  p_st << setprecision(3) << p_start;
  p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;
  lagger << lag;
  output.open("Theoretical_Percol/DP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + "_Lag_" + lagger.str() + ".csv");
  // Creating CSV file in "Theoretical_Percol" sub-directory to store output data.

  //outputbeta.open("CrtExp/Beta_DP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");

	#pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
	{
	  float percolation_probabilities_private[divisions] = {0};

	  #pragma omp for
		for (int i=0; i < divisions; i++){

			int seed = std::random_device{}();
			rng.seed(seed);

			percolation_probabilities_private[i] = theoretical_percolation_probability_dp(grid_size, birth_probabilities[i], r_init, number_of_census, lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            percolation_probabilities[n] += percolation_probabilities_private[n];
	        }
	    }
	}

  /*double p_c = 0.592746;
	//Percolation point for 2D Grid Network.
	std::vector<double> lnbase; 	// Stores ln(p - pc) values for various values of p.

	// Storing ln| p - p_c | for p ----> p_c+
	for(int i=0; i<divisions; i++)
	{
		double x =0.0;
		if (birth_probabilities[i] < p_c)
		{
			x= 0.0; //Only values greater than p_c are stored.
		}
		else
		{
			x= log(birth_probabilities[i] - p_c);
		}

		lnbase.push_back(x);
	}*/

	for (int i=0; i< divisions; i++){
		cout << "p " << setprecision(7) << birth_probabilities[i] << " Strength Of Percolation (P_infty) (Theoretical) " << setprecision(4) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}

  // Writing results to CSV File.

  output << "# Birth Probability (p) , Strength Of Percolation (P_infty) (Theoretical)\n";
  for (int i=0; i< divisions; i++){
		output << setprecision(8) << birth_probabilities[i] << "," << setprecision(7) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}
  output.close();


  /*outputbeta << "# Birth Probability (p) , ln|p - p_c|, Strength Of Percolation (P_infty) (Theoretical)\n";
  for (int i=0; i< divisions; i++)
  {
      if (lnbase[i] != 0)
      {
        // In other words, for all p > p_c
        outputbeta << setprecision(7) << birth_probabilities[i] << ","<< setprecision(10) << lnbase[i] << "," << setprecision(10) << percolation_probabilities[i] << endl;
    		/* Prints parameter value, ln|p - p_c|. and percolation probabilities to the terminal.
        The ordering is the same as the order of p in birth_probabilities
      }

	}
  outputbeta.close(); */

}

double theoret_percol_prob_denovo_dp(int grid_size, float birth_probability,  int r_init, int number_of_census, int lag)
{
  //Finds a Theoretical Estimate Of Percolation Probability For DP.

  int frame[grid_size*grid_size];
	//random_frame(frame, grid_size); // Initialize a random frame
  int updates_per_site =8000; // DP reaches steady state by 8000 for most values of p.
  if( birth_probability >= 0.4 && birth_probability >= 0.4)
  {
    //There is critical slowing down.
    updates_per_site =50000; // DP reaches steady state by 50000 for critical slowing down.
  }
  long limit = r_init*number_of_census;     //Stores the length of percolation_probabilities array necessary.
	float percolation_probabilities[limit];


	//simulate_np(frame,grid_size,birth_probability,updates_per_site); // simulate NP till it reaches steady state

  for(int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Assign a random frame
    simulate_dp(frame,grid_size,birth_probability,updates_per_site); // simulate DP till it reaches steady state
    //random_frame_of_density(birth_probability, frame, grid_size); // Assign a random frame of density p.

    for (int j = 0; j < number_of_census; ++j)
    {
      simulate_dp(frame,grid_size,birth_probability,lag); // simulate DP for an average of (lag) number of update per site

      int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0

      vector<cluster> clusters; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

      find_clusters_free_boundary(frame, labels, clusters, grid_size);
      // Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

      vector<coordinates> cluster_details;
      // Stores label id in x attribute, cluster size associated with label in 2nd coordinate.

      vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
      //Returns labels of spanning cluster(s) if present, -1 otherwise.
      float k =0.0; //If a spanning cluster is absent.
      if( spanning_cluster_labels[0] > -1)
      {
        // Spanning cluster in either direction present.
        k=1.0;
      }
      percolation_probabilities[i*number_of_census + j] = k;
    }

  }

  double percolation_probability = mean_of_array(percolation_probabilities,limit);
	// Pi[p] averaged over censuses

	return percolation_probability;
}


void calculate_pc_dp(int grid_size, float p_start, float p_end, int divisions, int r_init, int number_of_census, int lag){

	// Find "THEORETICAL" percolation probabilities for NP over a range of value of p by running parallel DP simulations

	vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

	double percolation_probabilities[divisions] = {0}; // initialize percolation probabilities

  ofstream output;
  // Creating a file instance called output to store output data as CSV.
  ofstream outputpc;

  std::vector<zd_coordinates> vec;
  // Stores collated output from parallel method calls in proper scending order of grid sizes.

  stringstream p_st, p_en, rini, lagger;

  p_st << setprecision(3) << p_start;
  p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;
  lagger << lag;
  output.open("Theoretical_Percol/Pc_DP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + "_Lag_" + lagger.str() + ".csv");
  // Creating CSV file in "Theoretical_Percol" sub-directory to store output data.

  //outputbeta.open("CrtExp/Beta_DP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");

	#pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
	{
	  double percolation_probabilities_private[divisions] = {0};

	  #pragma omp for
		for (int i=0; i < divisions; i++){

      stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
      message << "We are working on Occupy WS Prob:\t" << birth_probabilities[i] <<endl;
      cout << message.str();

			int seed = std::random_device{}();
			rng.seed(seed);

			percolation_probabilities_private[i] = theoret_percol_prob_denovo_dp(grid_size, birth_probabilities[i], r_init, number_of_census, lag);
		}
	    #pragma omp critical
	    {
	        for(int n=0; n < divisions; ++n) {
	            percolation_probabilities[n] += percolation_probabilities_private[n];
	        }
	    }
	}

	for (int i=0; i< divisions; i++){
		cout << "p " << setprecision(7) << birth_probabilities[i] << " Prob Of Existence Of Spanning Cluster (Pi(p)) " << setprecision(5) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}

  // Writing results to CSV File.

  output << "# Birth Probability (p) ,  Prob Of Existence Of Spanning Cluster (Pi(p)) \n";
  for (int i=0; i< divisions; i++){
		output << setprecision(8) << birth_probabilities[i] << "," << setprecision(10) << percolation_probabilities[i] << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}
  output.close();

}


void pavg_map_pc_dp(int grid_size, int r_init, int number_of_census, int lag)
{

  //Using Stauffer-Aharony to find infinite percolation threshold.
  int updates_per_site=50000;
  //vector<double> birth_probabilities = linspace(p_start, p_end, divisions); // create a vector of specified values of p

  //double rho_stat[divisions] = {0}; // initialize percolation probabilities

  ofstream output;
  // Creating a file instance called output to store output data as CSV.
  ofstream outputpc;

  std::vector<zd_coordinates> vec;
  // Stores collated output from parallel method calls in proper scending order of grid sizes.

  stringstream g, rini;

  /**p_st << setprecision(3) << p_start;
  p_en << setprecision(3) << p_end; */
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;
  g << grid_size;

  output.open("Theoretical_Percol/DP_L_"+ g.str() +  "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");
  // Creating CSV file in "Theoretical_Percol" sub-directory to store output data.

  //outputbeta.open("CrtExp/Beta_DP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");

  #pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
  {
    std::vector<zd_coordinates> vec_private;

    //Grants a static schedule with a chunk size of 1.
    /* Based on procedure suggested in:
    https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

    #pragma omp for nowait schedule(static)
    for (int i=0; i < r_init; i++)
    {

      stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
      message << "We are working on Occupy WS Trial:\t" << i <<endl;
      cout << message.str();

      int frame[grid_size*grid_size];
      int seed = std::random_device{}();
      rng.seed(seed);
      random_frame(frame, grid_size); // Assign a random frame
      simulate_dp(frame,grid_size,0.5,updates_per_site);

      zd_coordinates temp = binsearch_p_c(0.5, frame, grid_size, number_of_census, seed);
      //Stores [L,p, p^2] value for a particular trial.

      vec_private.push_back(temp);
      //rho_stat_private[i] = equilibrium_density_dp(grid_size, birth_probabilities[i], r_init, number_of_census, lag, updates_per_site, 0);
      // No collection of frames
    }
    #pragma omp for schedule(static) ordered
    for(int i=0; i< omp_get_num_threads(); i++)
    {
      #pragma omp ordered
        vec.insert(vec.end(), vec_private.begin(), vec_private.end());
        // Inserting p-values for each trial in order.
        stringstream message3;
        message3 << "Is this happening?\n";
        cout << message3.str();

    }
  }
  cout << " L, # Tr No , p\n";
  for (int i=0; i<vec.size(); i++){
    cout << vec[i].x << " " << setprecision(5) << i << " " << setprecision(10) << vec[i].y << endl;
    // Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
  }

  // Writing results to CSV File.

  output << "# L, # Tr No , p\n";
  for (int i=0; i<vec.size(); i++){
    output << vec[i].x << "," << setprecision(5) << i << ","  << setprecision(10) << vec[i].y << endl;
    // Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
  }
  output.close();
}


void pavg_map_pc_tcp(double q, int grid_size, int r_init, int number_of_census, int lag)
{

  //Using Stauffer-Aharony to find infinite percolation threshold in TCP set-up.
  int updates_per_site=50000;

  ofstream output;
  // Creating a file instance called output to store output data as CSV.
  ofstream outputpc;

  std::vector<zd_coordinates> vec;
  // Stores collated output from parallel method calls in proper scending order of grid sizes.

  stringstream g, q_no, rini;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  rini << r_init;
  q_no << q;
  g << grid_size;

    output.open("Theoretical_Percol/TCP_L_"+ g.str() + "_Q_"+ q_no.str() +  "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");
  // Creating CSV file in "Theoretical_Percol" sub-directory to store output data.

  //outputbeta.open("CrtExp/Beta_DP_L_"+ std::to_string(grid_size) + "_p1_" + p_st.str() + "_p2_" + p_en.str() + "_Cen_"+ std::to_string(number_of_census) + "_R_"+ rini.str() + ".csv");

  #pragma omp parallel // The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
  {
    std::vector<zd_coordinates> vec_private;

    //Grants a static schedule with a chunk size of 1.
    /* Based on procedure suggested in:
    https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

    #pragma omp for nowait schedule(static)
    for (int i=0; i < r_init; i++)
    {

      stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
      message << "We are working on Occupy WS Trial:\t" << i <<endl;
      cout << message.str();

      int frame[grid_size*grid_size];
      int seed = std::random_device{}();
      rng.seed(seed);
      random_frame(frame, grid_size); // Assign a random frame
      simulate_tp(frame,grid_size,0.5, q, updates_per_site);

      zd_coordinates temp = binsearch_p_c_TCP(0.5, q, frame, grid_size, number_of_census, seed);
      //Stores [L,p, p^2] value for a particular trial.

      vec_private.push_back(temp);
      //rho_stat_private[i] = equilibrium_density_dp(grid_size, birth_probabilities[i], r_init, number_of_census, lag, updates_per_site, 0);
      // No collection of frames
    }
    #pragma omp for schedule(static) ordered
    for(int i=0; i< omp_get_num_threads(); i++)
    {
      #pragma omp ordered
        vec.insert(vec.end(), vec_private.begin(), vec_private.end());
        // Inserting p-values for each trial in order.
        stringstream message3;
        message3 << "Is this happening?\n";
        cout << message3.str();

    }
  }
  cout << "q,  L, # Tr No , p\n";
  for (int i=0; i<vec.size(); i++){
    cout << q << " " << vec[i].x << " " << setprecision(5) << i << " " << setprecision(10) << vec[i].y << endl;
    // Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
  }

  // Writing results to CSV File.

  output << "#q, L, # Tr No , p\n";
  for (int i=0; i<vec.size(); i++){
    output << q << "," << vec[i].x << "," << setprecision(5) << i << ","  << setprecision(10) << vec[i].y << endl;
    // Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
  }
  output.close();
}

//----------------------------- Function For Calculating Custom ACF Data----------------------------------------//



float xdynamic(int frame[], int grid_size)
{
  //Calculates & Returns X dynamic variable for a frame instance.

  float x=0.0;
  for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
          {
            x += (float)(i*grid_size+j)*frame[i*grid_size + j];
          }
    }

    return x;

}



void acf_np(int grid_size, vector<f_coordinates>& acf_data, float p, int length, int lag)
{
  lag=1;
  int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 10000; // NP reaches steady state by 10000 for all values of p

  simulate_np(frame, grid_size, p ,updates_per_site); // simulate NP till it reaches steady state

  float r = (float) grid_size*grid_size;

  float mean= p*r*(r + 1)/2.0;
  float var = p*(1-p)*r*(r + 1 )*(2*r + 1)/6.0;

  stringstream msg;     //To make cout thread-safe as well as non-garbled due to race conditions.
  msg << "Mean & Variance this particular trial:\t" << setprecision(3) << mean << "  " << setprecision(3) << var << "\n";
  cout << msg.str();


  float X_tk[length];    // Stores value of X(t+k) for different values of 0<= k < length
  /* X is a dynamic variable (check notebook for details) based on frame (NOTE: NOT INJECTIVE), with mean and variance
  as above. */

  for(int i=0; i<length; i++)
  {
    X_tk[i] = xdynamic(frame, grid_size);  // Stores value of X_t (dynamic variable) for current frame (t).

    simulate_np(frame,grid_size, p , 1); // simulate NP for an average of (lag = 1) number of update per site
  }

  for(int i=0; i<length; i++)
  {
    int m = length - i;
    float sum =0.0;
    for(int j=0; j < m; j++)
    {
      sum += (X_tk[j] - mean)*(X_tk[j+i] - mean);
    }
    f_coordinates absalom;
    absalom.x = (float) i;
    absalom.y = (float) sum/(m*var);
    acf_data.push_back(absalom);

  }
  //Using estimator as per: https://en.wikipedia.org/wiki/Autocorrelation

  stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
  message << "Length of acf_data for this particular trial:\t" << acf_data.size() << "\n";
  cout << message.str();


}


void acf_np_custom(int grid_size, float p, int divisions, int length, int lag)
{
  ofstream output;
  // Creating a file instance called output to store output data as CSV.

  stringstream peon, div;

  peon << setprecision(3) << p;
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  div << divisions;
  output.open("ACF/NP_L_"+ std::to_string(grid_size) + "_p_" + peon.str() + "_Div_" + div.str() + "_Len_"+ std::to_string(length) + ".csv");
  // Creating CSV file in "ACF" sub-directory to store output data

  std::vector<f_coordinates> vec;
  //Will store ACF(t) data for every t <length for every division (trial) in order.

  //The implementation below is obtain a order of percolation probabilities that shadows the order of birth_probabilities
  #pragma omp parallel
  {
      std::vector<f_coordinates> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

      #pragma omp for nowait schedule(static)
      for (int i=0; i < divisions; i++)
      {

        int seed = std::random_device{}();
        rng.seed(seed);

        std::vector<f_coordinates> acf_data;

        acf_np(grid_size, acf_data, p, length, lag);
        //Finds and stores ACF functional form based on custom (NON-INJECTIVE) dynamic variable defined in xdynamic()

        vec_private.insert(vec_private.end(), acf_data.begin(), acf_data.end());

      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting ACF data for each trial in order.
      }
  }


  //Printing out obtained results.
  cout << "The vector elements are: "<< endl;
  cout << "# Tau (T) , ACF(T)\n";
  for (int i = 0; i < vec.size(); i++)
  {
    cout << vec[i].x << "  " << setprecision(3) << vec[i].y << endl;
  }

  // Writing results to CSV File.

  output << "# Tau (T) , ACF(T)\n";
  for (int i=0; i< vec.size(); i++){
		output << vec[i].x << "," << setprecision(4) << vec[i].y << endl;
		// Prints parameter value and percolation probabilities to the terminal. The ordering is the same as the order of p in birth_probabilities
	}
  output.close();


}



//------------------------------------ Finite Scaling Gimmickery (TCP) ---------------------------------

void output_dump(int frame[], int grid_size, int i, double p, int labels[], vector <cluster>& clusters, vector<int>& spanning_cluster_labels)
{
  ofstream ofputter;
  stringstream peon, div ,g;

  peon << setprecision(8) << p;
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  div << i;
  g << grid_size;
  ofputter.open("dump/SP_p_" + peon.str() + "_G_"+ g.str() + "_i_" + div.str() + ".csv");
  // Creating CSV file in "dump" sub-directory to store frame data

  for (int j = 0; j < grid_size; ++j)
    {
        for (int k = 0; k < grid_size; ++k)
        {
            int flag=0;
            for(int r =0 ; r < spanning_cluster_labels.size(); r++)
            {
              if(labels[j*grid_size + k] == spanning_cluster_labels[r])
              {
                //If current node is part of spanning cluster, trip up the flag variable
                flag=1;
              }
            }
            if ( flag == 1)
            {
              //Node is part of spanning cluster, outdata as 2.
              ofputter << -1 << ',';
            }
            else
            { ofputter << labels[j*grid_size + k] << ','; }
        }
        ofputter << '\n';
    }
  ofputter.close();
  }

zd_coordinates binsearch_p_c_TCP(double p, double q, int frame[], int grid_size, int num, int seed)
  {
    int updates_per_site =50000;
    // Using binary search to find first instance of grid percolation.
    double p_old= 2*p;
    for( int i=0; i < num; i++)
    {
      int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0
      vector <cluster> clusters;
      find_clusters_free_boundary(frame, labels, clusters, grid_size);
      // Performing DFS, finding clusters.

      vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
      if (spanning_cluster_labels[0] != -1)
      {
        //There exists a spanning cluster, search for cluster in lower half-interval.
        double temp = p;
        p -= fabs( p_old - temp )/2;
        p_old = temp;
      }
      else
      {
        // No spanning cluster found, search for spanning cluster in upper half-interval.
        double temp = p;
        p += fabs( p_old - temp )/2;
        p_old = temp;
      }
      rng.seed(seed);
      //random_frame_of_density(p, frame, grid_size);
      // Generating new random frame (with same seed) in half-interval

      random_frame(frame, grid_size); // Assign a random frame (with same seed) in half-interval
      simulate_tp(frame,grid_size,p, q, updates_per_site);
      // Simulate DP to equilibrium with p
    }

    zd_coordinates arkham; //Dummy variable.
    arkham.x = grid_size;
    arkham.y = p;
    arkham.z = (double) p*p;

    return arkham;

  }

void crtexp_gamma_TCP(int grid_size,vector<zd_coordinates> &comp_data, double p, double q, int r_init, int number_of_census, int lag)
  {
    // Using the formal definition of the average cluster size.
    //int lag=5000;
    tau_patch_size_find_tcp(grid_size, comp_data, p, q, r_init, number_of_census, lag);
    // Returns comp_data in the form of (p, s, n_s(p,q)).
    // We need to modify p with L.

    for(int i=0; i < comp_data.size(); i++)
    {
      if(comp_data[i].z >= 1)
      {
        // i.e Not A Header File [{#, p, q}] OR row denoting 0 P(p) [{p, -10, 0}]
        comp_data[i].x = grid_size;
        comp_data[i].z /= (grid_size*grid_size);
      }
      if(comp_data[i].y < 0)
      {
        // We have ourselves a null ('0') percolation strength newline [{p, -10, 0}]
        comp_data[i].x = grid_size;
      }
    }

    stringstream msg;
    msg << "Hola\n";
    cout << msg.str();
  }

void finite_scaling_crtexp_TCP(int grid_sizes[], double p, double q, string type, int divisions, int r_init, int number_of_census, int lag)
{
    //This method will find finite scaling relations for a given p ----------------> p_c.

    //int r_init =25;
    //number_of_census= 1;

    type= "Gam";

    std::vector<zd_coordinates> vec;
    // Stores collated output from parallel method calls in proper scending order of grid sizes.
    ofstream outputfinsc;
    // Creating a file instance called output to store output data as CSV.

    stringstream peon, quint, div ,g1, g2;

    peon << setprecision(4) << p;
    quint << setprecision(3) << q;
    //p_en << setprecision(3) << p_end;
    // setprecision() is a stream manipulator that sets the decimal precision of a variable.
    div << divisions;
    g1 << grid_sizes[0];
    g2 << grid_sizes[divisions-1];
    if (type == "Nu" || type == "Gam")
      {
        outputfinsc.open("CrtExp/FinSc"+ type +"_TCP_p_" + peon.str() + "_q_" + quint.str() + "_Div_" + div.str() + "_G1_"+ g1.str() + "_G2_"+ g2.str() + ".csv");
      }
    else
      {
        //Default is Beta/Gamma Computation.
        outputfinsc.open("CrtExp/FinSc_TCP_p_" + peon.str() + "_q_" + quint.str() + "_Div_" + div.str() + "_G1_"+ g1.str() + "_G2_"+ g2.str() + ".csv");
      }

    // Creating CSV file in "ACF" sub-directory to store output data

    #pragma omp parallel
    {
        std::vector<zd_coordinates> vec_private;

        //Grants a static schedule with a chunk size of 1.
        /* Based on procedure suggested in:
        https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

        #pragma omp for nowait schedule(static)
        for (int i=0; i < divisions; i++)
        {
          //type="Gam";
          stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
          message << "We are working on Grid Size:\t" << grid_sizes[i] <<endl;
          cout << message.str();
          int seed = std::random_device{}();
          rng.seed(seed);

          std::vector<zd_coordinates> comp_data;

          if(type == "Nu")
          {
            //crtexp_nu(grid_sizes[i], comp_data, r_init, number_of_census);
            continue;
            //Finds and returns nu (critical exponent) related data for a given grid_size
            // Returns {L, p, p^2}
          }
          else if(type =="Gam")
          {
            crtexp_gamma_TCP(grid_sizes[i], comp_data, p, q, r_init, number_of_census, lag);
            //Finds and returns gamma (critical exponent) related data for a given grid_size
            // Returns {L, p, p^2}
          }

          vec_private.insert(vec_private.end(), comp_data.begin(), comp_data.end());

        }

        #pragma omp for schedule(static) ordered
        for(int i=0; i< omp_get_num_threads(); i++)
        {
          #pragma omp ordered
            vec.insert(vec.end(), vec_private.begin(), vec_private.end());
            // Inserting critical exponent data for each grid size in order.
            stringstream message3;
            message3 << "Is this happening?\n";
            cout << message3.str();

        }
    }
    cout << "Game Over" << endl;
    vector <vector<double>> output;
    //Creating 2D vector to store final output

    if(type == "Gam")
    {
      double trialno =1; double q_anon=q;
      for(int i=0; i <vec.size(); i++)
      {
        if(vec[i].z < 1 && vec[i].y < 1 && vec[i].y >=0)
        {
          //Header line represented by row vector [# No, p, q]
          trialno = vec[i].x; q_anon= vec[i].z;
          continue;
        }
        output.push_back({p, q_anon, vec[i].x, trialno, vec[i].y, vec[i].z});
        //Filled as p,q, L, # No, s, n_s(p).
      }
    }
    else
    {
      double trialno =1;
      for(int i=0; i <vec.size(); i++)
      {
        int j = i % r_init;
        int trl_no = j % number_of_census + 1; int r_no = int(j/number_of_census);
        trialno = r_no*number_of_census + trl_no;

        output.push_back({p, q, vec[i].x, trialno, vec[i].y, vec[i].z});
        //Filled as p, q, L, # No, P(p), S(p)
      }
    }


    //Printing out obtained results.
    cout << "The vector elements are: "<< endl;
    //cout << " p , L, # Tr No , P[p] ,  S[p]\n";
    if(type == "Nu")
    {
      cout << " p , q, L, # Tr No , p ,  p^2 \n";
      outputfinsc << " p , q, L, # Tr No , p ,  p^2 \n";
    }
    else if(type == "Gam")
    {
      cout << " p , q, L, # Tr No , s ,  n_s(p) \n";
      outputfinsc << " p , q, L, # Tr No , s ,  n_s(p) \n";
    }
    else
    {
      cout << " p ,q, L, # Tr No , s ,  ns_(p)\n";
      outputfinsc << " p ,q, L, # Tr No , P[p] ,  S[p]\n";
    }
    for (int i = 0; i < output.size(); i++)
    {
      cout << setprecision(8) << output[i][0] << "  " << setprecision(3) << output[i][1] << "  " << setprecision(5) << output[i][2] << "  " << setprecision(3) << output[i][3] << "  " << setprecision(16) << output[i][4] << "  " << setprecision(16) << output[i][5] << endl;
    }
    // Saving to aforementioned CSV

    for (int i = 0; i < output.size(); i++)
    {
      outputfinsc << setprecision(8) << output[i][0] << "," << setprecision(5) << output[i][1] << "," << setprecision(8) << output[i][2] << "," << setprecision(3) << output[i][3] << "," << setprecision(16) << output[i][4] << "," << setprecision(16) << output[i][5] << endl;
    }
    outputfinsc.close();


  }

//--------------------------------- Crtical Exponents (Beta, Gamma etc) [Finite Scaling] [NP/DP]------------------------//

zd_coordinates binsearch_p_c(double p, int frame[], int grid_size, int num, int seed)
{
  int updates_per_site =50000;
  // Using binary search to find first instance of grid percolation.
  double p_old= 2*p;
  for( int i=0; i < num; i++)
  {
    int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0
    vector <cluster> clusters;
    find_clusters_free_boundary(frame, labels, clusters, grid_size);
    // Performing DFS, finding clusters.

    vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
    if (spanning_cluster_labels[0] != -1)
    {
      //There exists a spanning cluster, search for cluster in lower half-interval.
      double temp = p;
      p -= fabs( p_old - temp )/2;
      p_old = temp;
    }
    else
    {
      // No spanning cluster found, search for spanning cluster in upper half-interval.
      double temp = p;
      p += fabs( p_old - temp )/2;
      p_old = temp;
    }
    rng.seed(seed);
    //random_frame_of_density(p, frame, grid_size);
    // Generating new random frame (with same seed) in half-interval

    random_frame(frame, grid_size); // Assign a random frame (with same seed) in half-interval
    simulate_dp(frame,grid_size,p,updates_per_site);
    // Simulate DP to equilibrium with p
  }

  zd_coordinates arkham; //Dummy variable.
  arkham.x = grid_size;
  arkham.y = p;
  arkham.z = (double) p*p;

  return arkham;

}

void crtexp_nu(int grid_size,vector<zd_coordinates> &comp_data, int r_init, int number_of_census)
{
  // Based on methods highlighted in pages 73-74 of Stauffer & Anthony. Finding <p> & <p^2>.

  number_of_census=12; int updates_per_site =8000; //DP reaches equilibrium by 50000.

  int frame[grid_size*grid_size];
  for (int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame(frame, grid_size); // Assign a random frame
    simulate_dp(frame,grid_size,0.5,updates_per_site);
    // Simulate DP to equilibrium with p=0.5
    //random_frame_of_density(0.5, frame, grid_size); // Initialize a random frame of density 0.5 .
    zd_coordinates temp = binsearch_p_c(0.5, frame, grid_size, number_of_census, seed);
    comp_data.push_back(temp); //Writing up the values for a single trial in a given grid.
    // Stores data in the form of {L, p, p^2}
  }

}

void crtexp_gamma(int grid_size,vector<zd_coordinates> &comp_data, double p, int r_init, int number_of_census)
{
  // Using the formal definition of the average cluster size.
  int lag=15;
  tau_patch_size_find_np(grid_size, comp_data, p, r_init, number_of_census, lag);
  // Returns comp_data in the form of (p, s, ns(p)).
  // We need to modify p with L.

  for(int i=0; i < comp_data.size(); i++)
  {
    if(comp_data[i].z >= 1)
    {
      // i.e Not A Header File [{#, p, p}] OR row denoting 0 P(p) [{p, -10, 0}]
      comp_data[i].x = grid_size;
      comp_data[i].z /= (grid_size*grid_size);
    }
    if(comp_data[i].y < 0)
    {
      // We have ourselves a null ('0') percolation strength newline [{p, -10, 0}]
      comp_data[i].x = grid_size;
    }
  }

  stringstream msg;
  msg << "Hola\n";
  cout << msg.str();
}

void crtexp_beta_gamma(int grid_size,vector<zd_coordinates> &comp_data, double p, int r_init, int number_of_census)
{

  number_of_census=1;

  int frame[grid_size*grid_size];
	random_frame(frame, grid_size); // Initialize a random frame
	int updates_per_site = 10000; // NP reaches steady state by 10000 for all values of p
  int limit = r_init*number_of_census;
	float average_cluster_size[limit]; double percolation_probabilities[limit];

  for (int i=0; i< r_init; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    random_frame_of_density(p, frame, grid_size); // Initialize a random frame of density p.
    //simulate_np(frame,grid_size,birth_probability,updates_per_site);
    // simulate NP till it reaches steady state
    for (int j=0; j<number_of_census; j++)
    {
      float k =0.0; //Stores size of spanning cluster(s).

     	//simulate_np(frame,grid_size,birth_probability,lag); // simulate NP for an average of (lag) number of update per site
     	average_cluster_size[i*number_of_census+ j] = find_average_cluster_size(frame, grid_size);
      // find the average cluster size of the current frame (minus spanning cluster if present)

      int labels[grid_size*grid_size] = {0}; // initialize all sites of label lattice to 0
    	vector <cluster> clusters;

    	find_clusters_free_boundary(frame, labels, clusters, grid_size);
    	// Segregates clusters, populates labels lattice and accumulates clusters with free boundary conditions

    	vector<float> cluster_sizes; // See cluster_dynamics.h for the data structure cluster. It has two attributes: label and coords.

      vector<int> spanning_cluster_labels = spanning_cluster_label_id(frame, grid_size, labels, clusters);
      //Stores list of spanning cluster IDs, -1 otherwise (if no spanning clusterfuck present).

      if (spanning_cluster_labels[0] != -1)
      {
        //Spanning cluster(s) present
        for (int r =0 ; r < spanning_cluster_labels.size(); r++)
        {
          k += (float) clusters[spanning_cluster_labels[r]-1].coords.size();
        }
      }
      //Stores the number of vertices that belong to 2D spanning tree.
      percolation_probabilities[i*number_of_census + j] = k/(grid_size*grid_size);

      /*if ( i % 50 == 0)
      {
        output_dump(frame, grid_size, i, p, labels, clusters, spanning_cluster_labels);
        // Dumps frame as CSV file
      } */

    }
  }
  for (int i=0; i< limit; i++)
  {
    zd_coordinates temp; //Creating a temporary variable.
    //int trl_no = i % number_of_census + 1; int r_no = int(i/number_of_census);
    temp.x = grid_size;
    temp.y = percolation_probabilities[i];
    temp.z = average_cluster_size[i];
    comp_data.push_back(temp);
  }
}

void finite_scaling_crtexp(int grid_sizes[], double p, string type, int divisions, int r_init, int number_of_census)
{
  //This method will find finite scaling relations for a given p ----------------> p_c.

  //int r_init =25;
  //number_of_census= 1;

  type= "Nu";

  std::vector<zd_coordinates> vec;
  // Stores collated output from parallel method calls in proper scending order of grid sizes.
  ofstream outputfinsc;
  // Creating a file instance called output to store output data as CSV.

  stringstream peon, div ,g1, g2;

  peon << setprecision(4) << p;
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  div << divisions;
  g1 << grid_sizes[0];
  g2 << grid_sizes[divisions-1];
  if (type == "Nu" || type == "Gam")
    {
      outputfinsc.open("CrtExp/FinSc"+ type +"_DP_p_" + peon.str() + "_Div_" + div.str() + "_G1_"+ g1.str() + "_G2_"+ g2.str() + ".csv");
    }
  else
    {
      //Default is Beta/Gamma Computation.
      outputfinsc.open("CrtExp/FinSc_DP_p_" + peon.str() + "_Div_" + div.str() + "_G1_"+ g1.str() + "_G2_"+ g2.str() + ".csv");
    }

  // Creating CSV file in "ACF" sub-directory to store output data

  #pragma omp parallel
  {
      std::vector<zd_coordinates> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

      #pragma omp for nowait schedule(static)
      for (int i=0; i < divisions; i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on Grid Size:\t" << grid_sizes[i] <<endl;
        cout << message.str();
        int seed = std::random_device{}();
        rng.seed(seed);

        std::vector<zd_coordinates> comp_data;

        if(type == "Nu")
        {
          crtexp_nu(grid_sizes[i], comp_data, r_init, number_of_census);
          //Finds and returns nu (critical exponent) related data for a given grid_size
          // Returns {L, p, p^2}
        }
        else if(type =="Gam")
        {
          crtexp_gamma(grid_sizes[i], comp_data, p, r_init, number_of_census);
          //Finds and returns gamma (critical exponent) related data for a given grid_size
          // Returns {L, p, p^2}
        }

        vec_private.insert(vec_private.end(), comp_data.begin(), comp_data.end());

      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
          message3 << "Is this happening?\n";
          cout << message3.str();

      }
  }
  cout << "Game Over" << endl;
  vector <vector<double>> output;
  //Creating 2D vector to store final output

  if(type == "Gam")
  {
    double trialno =1;
    for(int i=0; i <vec.size(); i++)
    {
      if(vec[i].z < 1 && vec[i].y < 1 && vec[i].y >=0)
      {
        //Header line represented by row vector [# No, p, p]
        trialno = vec[i].x;
        continue;
      }
      output.push_back({p, vec[i].x, trialno, vec[i].y, vec[i].z});
      //Filled as p, L, # No, s, n_s(p).
    }
  }
  else
  {
    double trialno =1;
    for(int i=0; i <vec.size(); i++)
    {
      int j = i % r_init;
      int trl_no = j % number_of_census + 1; int r_no = int(j/number_of_census);
      trialno = r_no*number_of_census + trl_no;

      output.push_back({p, vec[i].x, trialno, vec[i].y, vec[i].z});
      //Filled as p, L, # No, P(p), S(p)
    }
  }


  //Printing out obtained results.
  cout << "The vector elements are: "<< endl;
  //cout << " p , L, # Tr No , P[p] ,  S[p]\n";
  if(type == "Nu")
  {
    cout << " p , L, # Tr No , p ,  p^2 \n";
    outputfinsc << " p , L, # Tr No , p ,  p^2 \n";
  }
  else if(type == "Gam")
  {
    cout << " p , L, # Tr No , s ,  s*n_s(p) \n";
    outputfinsc << " p , L, # Tr No , s ,  s*n_s(p) \n";
  }
  else
  {
    cout << " p , L, # Tr No , s ,  ns_(p)\n";
    outputfinsc << " p , L, # Tr No , P[p] ,  S[p]\n";
  }
  for (int i = 0; i < output.size(); i++)
  {
    cout << setprecision(8) << output[i][0] << "  " << setprecision(5) << output[i][1] << "  " << setprecision(3) << output[i][2] << "  " << setprecision(16) << output[i][3] << "  " << setprecision(16) << output[i][4] << endl;
  }
  // Saving to aforementioned CSV

  for (int i = 0; i < output.size(); i++)
  {
    outputfinsc << setprecision(8) << output[i][0] << "," << setprecision(8) << output[i][1] << "," << setprecision(3) << output[i][2] << "," << setprecision(16) << output[i][3] << "," << setprecision(16) << output[i][4]  << endl;
  }
  outputfinsc.close();


}


//-------------------------------- Critical Exponent Calculation [DP Model]----------------------------


void crtexp_DP_Basic(int grid_size,vector<zd_coordinates> &comp_data, double p, int r_init, int length)
{

  /* We are following the Dynamic Monte Carlo alogrithms highlighted in
     SECTION 3.4.3 (PG-- 867), Hirinschen 2000*/

  int frame[grid_size*grid_size];
  for (int i = 0; i < r_init ; i++)
  {
    int seed = std::random_device{}();
    rng.seed(seed);
    solitary_droplet(frame, grid_size); // Creates a lone wolf tear in the fabric of the lattice.

    zd_coordinates header; header.x= i+1; header.y = p; header.z = -1;
    comp_data.push_back(header); //Serves as a header row, seperating any two random trials.

    for(int j=0; j < length; j++)
    {
      //Looping over all the time steps, one at a time.
      int n = number_of_elem_of_array(frame, grid_size);
      //Checking if we are creamed or not.

      zd_coordinates temp; //Creating a temporary variable.
      temp.x = j;
      temp.y = 0;
      temp.z = 0;

      if( n == 0)
      {
        //There are no active elements left, hence all P(t) and N(t) values for all
        // subsequent time steps will be 0.

        for(int k =j; k <length; k++)
        {
          temp.x = k;
          comp_data.push_back(temp);
        }
        break;

      }
      else
      {
        temp.y = 1; //P(t) =1
        temp.z = n; //N(t) > 0

        comp_data.push_back(temp);
      }

      simulate_dp(frame, grid_size, p, 1);

    }

  }



}

void crtexp_dynamo_dp(int grid_size, double p_start, double p_end, int divisions, int r_init, int length)
{
  /* THE POINT OF THIS FUNCTION IS TO FIRSTLY ASCERTAIN P_C, FOLLOWED BY DYNAMIC ASCERTATION OF
     CERTAIN CRITICAL EXPONENTS (SECTION 3.4.3 (PG-- 867), Hirinschen 2000)*/

  vector<double> p_space = linspace(p_start, p_end, divisions);
  // The pspace to iterate over.

  std::vector<zd_coordinates> vec;
  // Stores collated output from parallel method calls in proper ascending order of p values.


  stringstream g, div ,p1, p2, rini;

  g << grid_size;
  //p_en << setprecision(3) << p_end;
  // setprecision() is a stream manipulator that sets the decimal precision of a variable.
  div << divisions;
  p1 << p_start;
  p2 << p_end;
  rini << r_init;

  //output_dp.open("CrtExp/P_c_DP_G_" + g.str() + "_Div_" + div.str() + "_p1_"+ p1.str() + "_p2_"+ p2.str() + "_R_"+ rini.str() + ".csv");
  // Creating CSV file.

  #pragma omp parallel
  {
      std::vector<zd_coordinates> vec_private;

      //Grants a static schedule with a chunk size of 1.
      /* Based on procedure suggested in:
      https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector */

      #pragma omp for nowait schedule(static)
      for (int i=0; i < p_space.size(); i++)
      {
        //type="Gam";
        stringstream message;     //To make cout thread-safe as well as non-garbled due to race conditions.
        message << "We are working on P Value:\t" << p_space[i] <<endl;
        cout << message.str();
        int seed = std::random_device{}();
        rng.seed(seed);

        std::vector<zd_coordinates> comp_data;

        crtexp_DP_Basic(grid_size, comp_data, p_space[i], r_init, length);

        vec_private.insert(vec_private.end(), comp_data.begin(), comp_data.end());

      }

      #pragma omp for schedule(static) ordered
      for(int i=0; i< omp_get_num_threads(); i++)
      {
        #pragma omp ordered
          vec.insert(vec.end(), vec_private.begin(), vec_private.end());
          // Inserting critical exponent data for each grid size in order.
          stringstream message3;
          message3 << "Is this happening?\n";
          cout << message3.str();

      }
  }

  vector <vector<double>> output; //Creating 2-D vector to store final output.

  ofstream output_dp;
  // Creating a file instance called output to store output data as CSV.

  double trialno =1; double p= p_space[0];

  stringstream pugalo; pugalo << setprecision(2) << p;
  int chk= mkdir(("CrtExp/"+pugalo.str()).c_str(), 0777);
  output_dp.open("CrtExp/"+pugalo.str()+"/P_c_DP_G_" + g.str() + "_Div_" + div.str() + "_p1_"+ p1.str() + "_R_"+ rini.str() + ".csv");
  output_dp << " p , # Tr No , t ,  P(t) ,  N(t) \n";

  for(int i=0; i <vec.size(); i++)
  {
    if(vec[i].z < 0 && vec[i].y < 1 && vec[i].y >= 0)
    {
      //Header line represented by row vector [# No, p, -1]

      if(vec[i].y > p)
      {
        //New p value being accessed.
        output_dp.close();
        trialno = vec[i].x; p = vec[i].y;
        stringstream p_new, p_final;
        p_new << setprecision(2) << p;
        p_final << p;

        int chk= mkdir(("CrtExp/"+p_new.str()).c_str(), 0777);
        output_dp.open("CrtExp/"+p_new.str()+"/P_c_DP_G_" + g.str() + "_Div_" + div.str() + "_p1_"+ p_final.str() + "_R_"+ rini.str() + ".csv");
        output_dp << " p , # Tr No , t ,  P(t) ,  N(t) \n";
      }

      trialno = vec[i].x; p = vec[i].y;
      continue;
    }
    output.push_back({p, trialno, vec[i].x, vec[i].y, vec[i].z});
    //Filled as p, # No, t,  P(t), N(t).
    output_dp << setprecision(8) << p << "," << setprecision(3) << trialno << "," << setprecision(5) << vec[i].x << "," << setprecision(3) << vec[i].y << "," << setprecision(12) << vec[i].z << endl;
  }

  cout << "The vector elements are: "<< endl;
  //cout << " p , L, # Tr No , P[p] ,  S[p]\n";
  cout << " p , # Tr No , t ,  P(t) ,  N(t) \n";


 for(int i=0; i <output.size(); i++)
 {
   if( i%4800 ==1)
   {
     cout << setprecision(8) << output[i][0] << "  " << setprecision(3) << output[i][1] << "  " << setprecision(5) << output[i][2] << "  " << setprecision(2) << output[i][3] << "  " << setprecision(10) << output[i][4] << endl;
   }
   //output_dp << setprecision(8) << output[i][0] << "," << setprecision(3) << output[i][1] << "," << setprecision(5) << output[i][2] << "," << setprecision(3) << output[i][3] << "," << setprecision(12) << output[i][4] << endl;

 }
 output_dp.close();

}
