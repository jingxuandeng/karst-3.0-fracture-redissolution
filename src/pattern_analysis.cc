#include "network.h"
#include "printing.h"


/**
* This function generates pattern starting from root with condition given by threshold.
* Then the pattern is transformed to the minimal tree spanning the pattern.
* Pore weight are set to the min sum from root of 1/diameter (resistance is also an option), pure d or resistance is alo an option.
* After generating minimal tree nodes and pore has consecutive flags:
* node.x <= 0 - node does't belong to the patter/tree
* node.x == n  - nod is a n-th generation (for root n=1)
* pore.x = 2   - pore belongs to the tree
* pore.x = 1   - pore belongs to the pattern but not belongs to the tree
* pore.x = 0   - pore do not belong to the pattern
* node.tmp is set to actual weight of a nod ("distance from a root")
*
* @param root node to be a root for a tree
* @param threshold threshold for pattern generation; pattern contains from connected pores with d>=threshold*d0
* @author Agnieszka Budek
* @date 26/04/2020
*/
void Network::find_minimal_spanning_tree(double threshold, Node* root){

	cerr<<"Finding minimal spanning tree..."<<endl;

	double big_number = 10e100;    //numerical infinity for tree generation algorithm
	int nodes_left    = 0;   //nr of nodes belonging to the longest pattern

	//Preparing dissolution pattern
	for (int i =0;i<NN;i++)   n[i]->x=0;
	for (int i =0;i<NP;i++)   p[i]->x=0;
	root->check_diss_pattern(threshold*d0);
	for(int i=0; i<NN; i++){
		if(n[i]->x==1)  n[i]->x=0;
		if(n[i]->x==2){ n[i]->x=1; nodes_left++;}
	}

	//finding tree
	for(int i=0;i<NN;i++) if(n[i]->x==1) n[i]->tmp=big_number;    //setting nod key to inf
	double min_k_value = big_number;   //minimal nod key
	Node * min_nod     = root; //first nod is a max root

	while(nodes_left>0){

		if(min_nod == root) {
			min_nod->tmp = 0;               //adding node to the tree
			min_nod->x   = 1;               //distance from the root+1
		}
		else{
			int min_s    = abs(min_nod->x)-1;            //pore to connect to the tree
			min_nod->x   = min_nod->n[min_s]->x+1;       //distance from the root+
			min_nod->p[min_s]->x=2;	                     //adding pore to the tree;
		}
		//updating keys
		for (int s=0; s< min_nod->b; s++){
			Node * nod_tmp = min_nod->n[s];
			if(nod_tmp->tmp==0)  continue;  //node added to the pattern (root)
			if(nod_tmp->x>1)     continue;  //node already added to the tree
			if(nod_tmp->x==0)    continue;  //node dose't belong to the pattern
			//double dist    = min_nod->tmp + 1./min_nod->p[s]->perm(mu_0); //weight is a distance from root
			double dist    = min_nod->tmp + 1./min_nod->p[s]->d; //weight is a distance from root
			if(nod_tmp->tmp > dist){ nod_tmp->tmp=dist; nod_tmp->x = -nod_tmp->neighbor_number(min_nod)-1;}
		}
		nodes_left --;
		if(nodes_left==0) break;
		//find new min node
		min_k_value = big_number;
		for(int i=0;i<NN;i++) if(n[i]->x<0)
			if (min_k_value > n[i]->tmp){
				min_nod     = n[i];
				min_k_value = n[i]->tmp;}
	}
}

/**
* This function finds the longest pattern (the one that achieves the largest y-coordinates).
* Then the minimal spanning tree is find for this pattern.
*
* @param threshold threshold for pattern generation; pattern contains from connected pores with d>=threshold*d0
* @author Agnieszka Budek
* @date 26/04/2020
*/

void Network::find_the_largest_tree(double threshold){

	int max_root      = 0;
	Node *n_max_y     = n[0];      //the tip of the longest pattern/tree
	double max_y      = -10e-10;
	double big_number = 10e100;    //numerical infinity for tree generation algorithm

	//looking for the longest  pattern
	for (int j=0;j<N_wi;j++){
		for (int i =0;i<NN;i++)   n[i]->x=0;
		for (int i =0;i<NP;i++)   p[i]->x=0;
		wi[j]->check_diss_pattern(threshold*d0);
		for(int i=0;i<NN;i++) if(n[i]->x==2) if(n[i]->xy.y>max_y) {
			max_y=n[i]->xy.y; n_max_y = n[i]; max_root=j;}
	}

	//generating the longest pattern
	for (int i =0;i<NN;i++)   n[i]->x=0;
	for (int i =0;i<NP;i++)   p[i]->x=0;
	wi[max_root]->check_diss_pattern(threshold*d0);

	//finding the best root
	double min_r=big_number;
	for(int i=0;i<N_wi;i++) if(wi[i]->x==2){
		Node * n_tmp = wi[i];
		for(int s=0;s<n_tmp->b;s++) if(n_tmp->p[s]->x==1)
			if(1./n_tmp->p[s]->perm(mu_0)) {min_r = 1./n_tmp->p[s]->perm(mu_0); max_root = i;}
		}


	find_minimal_spanning_tree(threshold,wi[max_root]);
}

/**
* This function returns length of the tree spanning dissolution pattern.
* The assumption is that the tree has been already generated (x flags are set correctly).
*
* @author Agnieszka Budek
* @date 26/04/2020
*/
double Network::find_minimal_tree_length(){

	cerr<<"Finding tree length..."<<endl;

	double min_dist=10e10;

	//for each outlet node
	for(int i=0;i<N_wo;i++){
		double dist_tmp = distance_to_root(wo[i]);
		if(dist_tmp<min_dist && dist_tmp>=0) min_dist = dist_tmp;
	}
	// if tree not connected to the outlet we find the largest child generation
	if(min_dist==10e10){
		int x_max=0; Node *n_max=NULL;
		for (int i=0; i<NN; i++) if(n[i]->x > x_max)   {x_max = n[i]->x;     n_max = n[i];}
		double y_max=0;
		for (int i=0; i<NN; i++) if(n[i]->x == x_max && n[i]->xy.y >y_max) {y_max = n[i]->xy.y;  n_max = n[i];}
		if (n_max != NULL) min_dist = distance_to_root(n_max);
	}

	return min_dist;
}

/**
* This function returns distance form given node (n_tmp)
* The assumption is that the tree has been already generated (x flags are set correctly).
*
* @param n_tmp node we are looking distance for
* @author Agnieszka Budek
* @date 26/04/2020
*/
double Network::distance_to_root(Node* n_tmp){

	if(n_tmp->x<=0) return -1;

	double	dist = 0;
	while (n_tmp->x>1)
		for(int s=0;s<n_tmp->b;s++) if(n_tmp->p[s]->x==2 && n_tmp->n[s]->x== n_tmp->x-1){
			dist+= point_distance(n_tmp->p[s]->n[0]->xy, n_tmp->p[s]->n[1]->xy);
			n_tmp=n_tmp->n[s];
		}

	return dist;

}

/**
* This function returns child distribution: it gives nr of children in n-th generation.
* The assumption is that the tree has been already generated (x flags are set correctly).
*
* @author Agnieszka Budek
* @date 26/04/2020
*/
int* Network::find_child_distribution (){

	//find max generation
	int max_child=1;
	for(int i=0;i<NN;i++) if(max_child<n[i]->x) max_child=n[i]->x;
	int* child_dist = new int[max_child+1];
	for(int i=0; i<=max_child; i++) child_dist[i]=0;

	//generate histogram
	for(int i=0;i<NN;i++) child_dist[int(n[i]->x)]++;

	child_dist[0] = max_child;
	return child_dist;

}

/**
* This function finds forking in a tree (starting form the root).
* To estimate fractal dimension it is better to use revert forking.
* As a result nodes flags are changed (node.tmp).
* The assumption is that the tree has been already generated (x flags are set correctly).
* node.tmp == 0  - node does't belong to the patter/tree
* node.tmp == n  - nod is a n-th generation fork (for root n=1)
* node.tmp = -n   - node belongs to the tree but is not a fork (is a child of for of n-th generation)
* @author Agnieszka Budek
* @date 26/04/2020
*/
void Network::find_forks(){

	cerr<<"Finding forks..."<<endl;

	//set initial flags
	for( int i=0; i<NN; i++) n[i]->tmp=0;

	Node* root = NULL;
	//find root
	for( int i=0; i<NN; i++) if(n[i]->x==1) {root = n[i]; break;}

	if(root == NULL)   {cerr<<"Problem in Network::find_forks(), root is not to be found."<<endl;}
	else               root->find_forks();


}


/**
* This function finds forking in a tree (starting form the leafs).
* As a result nodes flags are changed (node.tmp).
* The assumption is that the tree has been already generated (x flags are set correctly).
* node.tmp == 0  - node does't belong to the patter/tree
* node.tmp == n  - nod is a n-th generation fork (for root n=1)
* node.tmp = -n   - node belongs to the tree but is not a fork (is a child of for of n-th generation)
* @author Agnieszka Budek
* @date 26/04/2020
*/
double Network::find_reverse_forks(){

	cerr<<"Finding reverse forks..."<<endl;

	double max_horton=0;

	//set initial flags
	for( int i=0; i<NN; i++) if(n[i]->x > 0) n[i]->tmp=1;   else n[i]->tmp=0;

	int max_child=0;
	//find maximum child number
	for( int i=0; i<NN; i++) if(n[i]->x>max_child) max_child = n[i]->x;

	//for all children generation (starting form the largest one)
	while(max_child>0){

		//for all nodes in this generation
		for(int i=0; i<NN; i++) if(n[i]->x==max_child){

			Node *nn = n[i];
			int children_of_proper_order=0;

			//find nr of proper children
			for(int s=0; s<nn->b;s++)if(nn->p[s]->x==2)
				if(nn->n[s]->x == nn->x+1 && fabs(nn->tmp) == fabs(nn->n[s]->tmp))    children_of_proper_order++;

			//update node's stream number
			if(children_of_proper_order>1) nn->tmp = fabs(nn->tmp)+1;
			else                           nn->tmp = fabs(nn->tmp); //optionally one can add here minus to distinguish forks (this is why all those fabs are here)
			//update node's parent stream number
			for(int s=0; s<nn->b;s++)if(nn->p[s]->x==2)
				if(nn->n[s]->x == nn->x-1 && fabs(nn->n[s]->tmp)<fabs(nn->tmp))       nn->n[s]->tmp = fabs(nn->tmp);

			//updating max horton
			if(nn->tmp > max_horton) max_horton = nn->tmp;
		}
		max_child--;
	}

	return max_horton;
}



/**
* This function returns fork distribution: it gives nr of forks in n-th generation.
* The assumption is that the tree has been already generated (x flags are set correctly).
*
* @author Agnieszka Budek
* @date 26/04/2020
*/
int* Network::find_fork_distribution(){

	//find max generation
	int max_fork=1;
	for(int i=0;i<NN;i++) if(max_fork<n[i]->tmp) max_fork=n[i]->tmp;
	int* fork_dist = new int[max_fork+1];
	for(int i=0; i<=max_fork; i++) fork_dist[i]=0;

	//generate histogram
	for(int i=0;i<NN;i++) if(n[i]->tmp>0) fork_dist[int(n[i]->tmp)]++;

	fork_dist[0] = max_fork;
	return fork_dist;

}


/**
* This function returns avarage size of a cluster (belonging to the tree) of the length n.
* Randomly clusters roots are selected and then the size of a cluster is calculated.
* The assumption is that the tree has been already generated (x flags are set correctly).
* Clusters that are shorter than n are not taken into account.
* @author Agnieszka Budek
* @date 07/05/2020
*/
double Network::find_cluster_size(int l){


	int x_tot = 0;
	for(int i=0;i<NN;i++) if(n[i]->x>0) x_tot++;
	if(x_tot<sqrt(NN)) return 0;

	int N   = 1000;
	int ile = 0;
	int sum = 0;
	for (int i=0;i<N;i++){
		Node *n_tmp = NULL;
		while(n_tmp == NULL) {int j = rand()%NN; if(n[j]->x>0) n_tmp = n[j];}
		int size = n_tmp->cluster_size(l);
		if(size>0){
			sum+=size;
			ile++;
		}
	}

	if(ile>=100) return double(sum)/ile;
	else       return 0;
}



