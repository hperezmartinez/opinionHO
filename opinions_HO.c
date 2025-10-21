// The following code corresponds to the opinion formation model presented in "Social polarization promoted by
// higher order interactions" (https://doi.org/10.48550/arXiv.2507.12325)
// The most updated version of this code can be found in the repository: https://github.com/hperezmartinez/opinionHO
//
// This code generates <<N_stats>> opinion equilibrium states starting from random initial conditions for a given
// combination of parameters <<lambda>> and <<beta>>, taking <<lambda_T>> = gamma-<<lambda>>, with gamma = 20 in
// our particular case.
//
// It measures the average opinion, standard deviation, and higher-order exposure, as defined in the paper.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

int **A; // List of neighbors
int **A_neigh; // Auxiliary list of neighbors
int **A_T1, **A_T2; // List of neighbors through triangles
double **weights; // weights for pairwise interactions
double **weights_T; // weights for triangle interactions
int *v_degree; // pairwise degrees
int *v_degree_T; // Higher-order degrees
int N, L, L_T; // number of nodes, links, and triangles
double total_change; // For the early stopping condition in case of convergence to a stable state

#define N_stats 100 // Number of simulations to run
#define dt 0.1 // time step
#define N_max_steps (1000/dt)

// The following lines mark the network used
// Only one line can be uncommented at any given time
#define RSC
//#define HGMRI
//#define HG // For inclusiveness


// Random number generator, obtained and adapted from:
// Teukolsky, S. A., Flannery, B. P., Press, W., & Vetterling, W. (1992). Numerical recipes in C. SMR, 693(1), 59-70.

#define MBIG 1000000000
#define MZ 0
#define FAC (1.0/MBIG)
int MSEED;
long idum = -1;
// Any large MBIG, and any smaller (but still large) MSEED can be substituted for the above values.

double ran3(){
// Set idum to any negative value to initialize or reinitialize the sequence.
    static int inext,inextp;
    static long ma[56]; //The value 56 (range ma[1..55]) is special and should not be modified; see Knuth.
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    
    if (idum < 0 || iff == 0) { // Initialization
        iff=1;
        mj=labs(MSEED-labs(idum)); //Initialize ma[55] using the seed idum and the large number MSEED.
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        
        for (i=1;i<=54;i++) { // Now initialize the rest of the table in a slightly random order with numbers that are not especially random.
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++) //We randomize them by “warming up the generator”
            for (i=1;i<=55;i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
        }
        
        inext=0; // Prepare indices for our first generated number.
        inextp=31; // The constant 31 is special; see Knuth.
        idum=1;
    }
    
    // Here is where we start, except on initialization.
    if (++inext == 56) inext=1; // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp]; //Generate a new random number subtractively.
    if (mj < MZ) mj += MBIG; // Be sure that it is in range.
    ma[inext]=mj;  // Store it and output the derived uniform deviate.
    return mj*FAC;
}


void ini_ran3(){
    MSEED = (int)(time(NULL));
    //MSEED = 1234567890;
    ran3();
}


void create_header(FILE *f, char *network_type, int HO_degree, double mean_opinion, double dev_opinion,
                   double gamma, double lmbda, double lmbda_T, double beta){
    time_t current_time = time(NULL);
    char *date;

    date = ctime(&current_time);
    printf("Creating header...\t");
    fprintf(f, "# Opinion formation model with higher-order interactions (opinions_HO.c)\n");
    fprintf(f, "# Original code contained in https://github.com/hperezmartinez/opinionHO\n");
    fprintf(f, "# Results are reported in https://doi.org/10.48550/arXiv.2507.12325\n");
    fprintf(f, "#\n");
    fprintf(f, "# date: %s", date);
    fprintf(f, "#\n");
    fprintf(f, "# Network data:\n");
    fprintf(f, "#\t %s, d=%d, nodes: %d, links: %d, triangles: %d\n", network_type, HO_degree, N, L, L_T);
    fprintf(f, "#\n");
    fprintf(f, "# Initial opinion: %lf +- %lf\n", mean_opinion, dev_opinion);
    fprintf(f, "#\n");
    fprintf(f, "# Variables:\n");
    fprintf(f, "#\t gamma: %.2lf, lambda: %.3lf, lambda_T: %.3lf, beta: %.3lf\n", gamma, lmbda, lmbda_T, beta);
    fprintf(f, "#\n");
    fprintf(f, "# Simulation parameters:\n");
    fprintf(f, "#\t Max. steps: %.0lf, dt: %.2lf, N_stats: %d\n", N_max_steps, dt, N_stats);
    fprintf(f, "#\n");
    fprintf(f, "#==================================================\n");
    fprintf(f, "#\n");
    printf("Header created.\n");

}


void readNetwork(char *network_type, int *HO_degree){
    // Networks are read from files, not generated. After generation, nodes that don't belong to 
    // at least one hyperlink are elliminated, in order to discount any effects that different
    // types of interactions could have.

    int i, j, n, l, m, *v_degree_aux, trash;
    char filename[200];
    FILE *fin;

    // PAIRWISE INTERACTIONS
    #ifdef RSC
    *HO_degree = 3;
    sprintf(filename, "networks/RSC1899k=10d=%d_pairs.txt", *HO_degree);
    #endif // RSC
    #ifdef HGMRI
    *HO_degree = 7; // 15
    sprintf(filename, "networks/HGMRI2000k=10d=%d_pairs.txt", *HO_degree);
    #endif // HGMRI
    #ifdef HG
    *HO_degree = 3;
    double inclusiveness = 0.5; // 0.
    sprintf(filename, "networks/HG1899k=10d=%dI=%.1lf_pairs.txt", *HO_degree, inclusiveness);
    #endif // HG

    fin=fopen(filename, "r");

    if(fin == NULL){
        printf("Error, input file for pairwise interactions not found:\n");
        printf("%s\n", filename);
        return;
    }

    fscanf(fin, "%d %d", &N, &L); // The first line contains the number of nodes and number of links
    printf("Reading network:\n");
    printf("\tNetwork data (parwise): %d %d\n", N, L);
    
    A = (int**)malloc(N*sizeof(int*)); 
    weights = (double**)malloc(N*sizeof(double*));
    v_degree = (int*)calloc(N, sizeof(int));
    v_degree_aux = (int*)calloc(N, sizeof(int));

    for(i=0; i<L; i++){ // Reading as many lines as links there are, and obtaining degrees
        fscanf(fin, "%d %d", &n, &l);
        //printf("%d: %d %d\n", i, n, l);
        v_degree[n] += 1;
        v_degree[l] += 1;
    }

    // Sizing list of neighbors after reading the network
    for(i=0; i<N; i++){
        A[i]=(int*)calloc(v_degree[i], sizeof(int));
        weights[i] = (double*)calloc(v_degree[i], sizeof(double));
    }

    rewind(fin);
    fscanf(fin, "%d %d", &trash, &trash);
    
    // A[i][j] indicates the j-th neighbor of node i
    // weights[i][j] indicates the weight of agent i with its j-th neighbor
    for(i=0; i<L; i++){
        fscanf(fin, "%d %d", &n, &l);
        A[n][v_degree_aux[n]]=l;
        A[l][v_degree_aux[l]]=n;
        weights[n][v_degree_aux[n]] = 1.0; 
        weights[l][v_degree_aux[l]] = 1.0;
        v_degree_aux[n] += 1;
        v_degree_aux[l] += 1;
    }

    fclose(fin);


    // A_neigh[i][j] indicates the coordinate of agent i in the list of neighbors of j
    // In other words, A[j][A_neigh[i][j]] == i
    // This is done to accelerate search in the computing of the weights
    int neigh;
    A_neigh=(int**)malloc(N*sizeof(int*));

    for(i=0; i<N; i++){
        A_neigh[i]=(int*)calloc(v_degree[i], sizeof(int));

        for(j=0; j<v_degree[i]; j++){
            neigh = A[i][j];
            for(l=0; l<v_degree[neigh]; l++)
                if(A[neigh][l] == i) break;
            A_neigh[i][j] = l;
        }
    }
    

    
    /// HIGHER ORDER INTERACTIONS
    #ifdef RSC
    sprintf(filename, "networks/RSC1899k=10d=%d_triangles.txt", *HO_degree);
    #endif // RSC
    #ifdef HGMRI
    sprintf(filename, "networks/HGMRI2000k=10d=%d_triangles.txt", *HO_degree);
    #endif // HGMRI
    #ifdef HG
    // For the case with vaying inclusiveness, the subjacent triangle network is the one
    // corresponding to the RSC of N=1899, k=10, d=3.
    sprintf(filename, "networks/RSC1899k=10d=%d_triangles.txt", *HO_degree);
    #endif // HG

    fin=fopen(filename, "r");

    if(fin == NULL){
        printf("Error, input file for higher-order interactions not found:\n");
        printf("%s\n", filename);
        return;
    }

    fscanf(fin, "%d %d", &trash, &L_T); // The first line contains the number of nodes and number of triangles
    printf("\tNetwork data (HO interactions): %d %d\n", N, L_T);
    printf("\n");
    
    v_degree_T = (int*)calloc(N, sizeof(int));
    weights_T = (double**)malloc(N*sizeof(double*));
    
    A_T1 = (int**)malloc(N*sizeof(int*));
    A_T2 = (int**)malloc(N*sizeof(int*));
    
    for(i=0; i<N; i++)
        v_degree_aux[i] = 0;

    for(i=0; i<L_T; i++){
        fscanf(fin, "%d %d %d", &n, &l, &m);
        v_degree_aux[n] += 1;
        v_degree_aux[l] += 1;
        v_degree_aux[m] += 1;
        
    }

    for(i=0; i<N; i++)
        v_degree_T[i] = v_degree_aux[i];
    
    for(i=0; i<N; i++){
        A_T1[i] = (int*)malloc(v_degree_T[i]*sizeof(int));
        A_T2[i] = (int*)malloc(v_degree_T[i]*sizeof(int));
        weights_T[i] = (double*)malloc(v_degree_T[i]*sizeof(double));
        
        v_degree_aux[i] = 0;
    }
    
    rewind(fin);
    fscanf(fin, "%d %d", &trash, &trash);
    
    for(i=0; i<L_T; i++){
        fscanf(fin, "%d %d %d", &n, &l, &m);
        //printf("%d %d %d\n", n, l, m);
        
        A_T1[n][v_degree_aux[n]]=l;
        A_T2[n][v_degree_aux[n]]=m;
        
        A_T1[l][v_degree_aux[l]]=n;
        A_T2[l][v_degree_aux[l]]=m;
        
        A_T1[m][v_degree_aux[m]]=l;
        A_T2[m][v_degree_aux[m]]=n;
        
        
        weights_T[n][v_degree_aux[n]] = 1.0;
        weights_T[l][v_degree_aux[l]] = 1.0;
        weights_T[m][v_degree_aux[m]] = 1.0;
        
        v_degree_aux[n] += 1;
        v_degree_aux[l] += 1;
        v_degree_aux[m] += 1;
    }
    
    free(v_degree_aux);
    fclose(fin);
}


void free_all(double *opinions, double *new_opinions){
    int i;

    free(v_degree);
    free(v_degree_T);
    
    for(i=0; i<N; i++){
        free(A[i]);
        free(A_neigh[i]);
        free(A_T1[i]);
        free(A_T2[i]);

        free(weights[i]);
        free(weights_T[i]);
    }
    free(A);
    free(A_neigh);
    free(A_T1);
    free(A_T2);

    free(weights);
    free(weights_T);
    free(opinions);
    free(new_opinions);
}


void initialize_opinions(double *opinions, double *new_opinions, double lmbda, double lmbda_T){
    int i;
    double alea, max;

    if(lmbda+lmbda_T < 1) max = 1.;
    else max = lmbda+lmbda_T;

    for(i=0; i<N; i++){
        alea=2*max*ran3();
        opinions[i] = new_opinions[i] = alea-max;
    }
}


void calculate_mean_var_lf(double *data, int N_comp, double *mean, double *dev){

    if(N_comp == 0){
        *mean = 0;
        *dev = 0;
        return;
    }

    int i, conv_check;
    double sum, sum_squared;

    sum = 0;
    sum_squared = 0;
    conv_check = 0;

    for (i=0; i<N_comp; i++){
        sum += data[i];
        sum_squared += data[i]*data[i];
        if(data[i] != 0) conv_check += 1;
    }

    if(conv_check == 0){ // if all values are zero
        *mean = 0;
        *dev = 0;
        return;
    }
    else if(conv_check == 1){
        *mean = sum/N_comp;
        *dev = 0;
        return;
    }

    *mean = sum/N_comp;
    *dev = ((sum_squared/N_comp)-(*mean)*(*mean))/N_comp;

    if(*dev<0) *dev = 0; // To avoid numerical errors
    else *dev = sqrt(*dev);
}


void compute_weights(double *network, double *norm_factors, double beta, double delta){
    // This function computes the weights for pairwise interactions, and store them in the "weights" matrix
    // Note that the content of "weights" is not properly normalized to 1 for efficiency. Rather, normalization 
    // is performed by the "norm_factors" array in function "compute_derivative".

    int i, j;

    for(i=0; i<N; i++){
        norm_factors[i] = 0;
    }

    for(i=0; i<N; i++)
    {
        for(j=0; j<v_degree[i]; j++)
        {
            if(i<A[i][j]){
                // The matrix "A_neigh" allows for computing the numerator in the weights function only once, 
                // saving the information of the link i-j for both agents i, j at the same time, thus reducing
                // significantly the simulation time
                weights[i][j] = weights[A[i][j]][A_neigh[i][j]] = pow(fabs(network[i]-network[A[i][j]]) + delta, -beta);

                norm_factors[i] += weights[i][j];
                norm_factors[A[i][j]] += weights[A[i][j]][A_neigh[i][j]];
            }
        }
    }
}


void compute_weights_T(double *network, double *norm_factors, double beta, double delta_T){
    // This function computes the weights for higher-order interactions, and store them in the "weights-T" matrix
    // Note that the content of "weights_T" is not properly normalized to 1 for efficiency. Rather, normalization 
    // is performed by the "norm_factors_T" array in function "compute_derivative".
    int i, j;
    
    for(i=0; i<N; i++){
        norm_factors[i] = 0;
    }
    
    for(i=0; i<N; i++){
        for(j=0; j<v_degree_T[i]; j++){
            weights_T[i][j] = pow(fabs(network[A_T1[i][j]]-network[i])+fabs(network[A_T2[i][j]] - network[i]) + delta_T, -beta);
            norm_factors[i] += weights_T[i][j];
        }
    }
}


void compute_derivative(double *network, double *slope, double lmbda, double lmbda_T, double beta, double delta, double delta_T){
    int i, j;
    double neighbors_opinions, neighbors_opinions_T;

    double *norm_factors, *norm_factors_T, *tanhs;

    tanhs = (double*)malloc(N*sizeof(double));
    norm_factors = (double*)malloc(N*sizeof(double));
    norm_factors_T = (double*)malloc(N*sizeof(double));

    if(lmbda != 0) // Conditions are in place for efficiency
        compute_weights(network, norm_factors, beta, delta);
    if(lmbda_T != 0)
        compute_weights_T(network, norm_factors_T, beta, delta_T);

    for(i=0; i<N; i++)
        tanhs[i] = tanh(network[i]);
    
    for(i=0; i<N; i++)
    {
        neighbors_opinions = 0;
        if(lmbda != 0)
            for(j=0; j<v_degree[i]; j++)
                neighbors_opinions += weights[i][j]*tanhs[A[i][j]]/norm_factors[i];
        
        neighbors_opinions_T = 0;
        if(lmbda_T != 0)
            for(j=0; j<v_degree_T[i]; j++)
                neighbors_opinions_T += weights_T[i][j]*tanh((network[A_T1[i][j]]+network[A_T2[i][j]])/2.)/norm_factors_T[i];

        slope[i] = -network[i] + lmbda*neighbors_opinions + lmbda_T*neighbors_opinions_T; //Equation for the opinion change of agent i
    }
    free(tanhs);
    free(norm_factors);
    free(norm_factors_T);
}


void actualize_network(double *network, double *new_network, double lmbda, double lmbda_T, double beta, double delta, double delta_T){
    // This function computes the opinion change using a fourth-order Runge-Kutta method

    int i, j;
    double *slope, *new_position;

    slope = (double*)malloc(N*sizeof(double)); // Opinion change
    new_position = (double*)malloc(N*sizeof(double)); // Auxiliary array for the intermediate opinions needed for the Runge-Kutta method

    /// Step 1
    compute_derivative(network, slope, lmbda, lmbda_T, beta, delta, delta_T);
    for(i=0; i<N; i++){
        new_network[i] = network[i] + slope[i]*dt/6.0;
        new_position[i] = network[i] + slope[i]*dt/2.0;
    }
    /// Step 2
    compute_derivative(new_position, slope, lmbda, lmbda_T, beta, delta, delta_T);
    for(i=0; i<N; i++){
        new_network[i] += slope[i]*dt/3.0;
        new_position[i] = network[i] + slope[i]*dt/2.0;
    }
    /// Step 3
    compute_derivative(new_position, slope, lmbda, lmbda_T, beta, delta, delta_T);
    for(i=0; i<N; i++){
        new_network[i] += slope[i]*dt/3.0;
        new_position[i] = network[i] + slope[i]*dt;
    }
    /// Step 4
    compute_derivative(new_position, slope, lmbda, lmbda_T, beta, delta, delta_T);
    for(i=0; i<N; i++){
        new_network[i] += slope[i]*dt/6.0;
    }
    free(new_position);
    free(slope);

    total_change = 0;
    for(i=0; i<N; i++){
        total_change += fabs(new_network[i]-network[i]); // For the early stopping condition
        network[i] = new_network[i];
    }
}


int belong(int *array, int N_comp, int number){
    int i;

    for(i=0; i<N_comp; i++){
        if(array[i] == number) return 1;
    }
    return 0;
}


void compute_opinions(double *data, int N_comp, double *mean, double *dev){

    int i, conv_check;
    double sum, sum_squared;

    sum = 0;
    sum_squared = 0;
    conv_check = 0;

    for (i=0; i<N_comp; i++){
        sum += data[i];
        sum_squared += data[i]*data[i];
        if(data[i] != 0) conv_check += 1;
    }

    if(conv_check == 0){ // If all values are zero, to avoid numerical errors
        *mean = 0;
        *dev = 0;
        return;
    }
    else if(conv_check == 1){
        *mean = sum/N_comp;
        *dev = 0;
        return;
    }

    *mean = sum/N_comp;
    *dev = (sum_squared*1.0/N_comp)-(*mean)*(*mean);

    if(*dev<0) *dev = 0; // To avoid numerical errors
    else *dev = sqrt(*dev);
}


double compute_exposure_HO(double *network, int N){
    int i, j;
    double aux, mean_expo;

    mean_expo = 0;

    for(i=0; i<N; i++){
        if(network[i]>0){
            aux = 0;
            for(j=0; j<v_degree_T[i]; j++)
                if(network[A_T1[i][j]]+network[A_T2[i][j]]<0) aux += 1;
            
            if(v_degree_T[i] != 0)  aux /= v_degree_T[i];
            mean_expo += aux;
        }

        if(network[i]<0){
            aux = 0;
            for(j=0; j<v_degree_T[i]; j++)
                if(network[A_T1[i][j]]+network[A_T2[i][j]]>0) aux += 1;
            
            if(v_degree_T[i] != 0) aux /= v_degree_T[i];
            mean_expo += aux;
        }
    }
    mean_expo /= N;
    return mean_expo;
}



int main(int argc, char **argv)
{
    ini_ran3();

    double lmbda, lmbda_T, gamma, beta;
    double delta, delta_T;

    gamma = 20.;

    if(argc==1){
        lmbda = 10.;
        lmbda_T = gamma-lmbda;
        beta = 0.5;
        printf("No arguments given. Falling bach to default values:\n");

    }
    else if(argc != 3){
        printf("Error, 2 parameters must be passed:\n");
        printf("\t lmbda, beta\n");
        printf("Only %d were given\n", argc-1);
        return 0;
    }
    else{    
        sscanf(argv[1], "%lf", &lmbda);
        sscanf(argv[2], "%lf", &beta);
        lmbda_T = gamma-lmbda;

        printf("Arguments read from command console:\n");
    }
    printf("lmbda: %lf, lmbda_T: %lf, beta: %lf\n", lmbda, lmbda_T, beta);
    printf("\n");

    if(lmbda != 0) delta = 2*lmbda*0.001;
    else delta = 1E-20; // Irrelevant, because the interaction term is always multiplied by lmbda=0
    if(lmbda_T != 0) delta_T = 2*lmbda_T*0.001;
    else delta_T = 1E-20; // Irrelevant, because the interaction term is always multiplied by lmbda_T=0


    // Reading network
    char network_type[10];
    int HO_degree; // Its value is set in "readNetwork"

    #ifdef RSC
    sprintf(network_type, "RSC");
    #endif // RSC
    #ifdef HGMRI
    sprintf(network_type, "HGMRI");
    #endif // HGMRI    
    #ifdef HG // For the case with inclusiveness
    sprintf(network_type, "HG");
    #endif // HG

    readNetwork(network_type, &HO_degree);

    if(N == 0){
        printf("Error, network not read correctly.\n");
        return 0;
    }

    
    // Initializing opinions
    double *opinions, *new_opinions;
    opinions = (double*)calloc(N, sizeof(double));
    new_opinions = (double*)calloc(N, sizeof(double));

    int i, j, l, n_steps;
    double mean_opinion, dev_opinion, initial_opinion, initial_dev_opinion;
    int n_stat, polarized;
    double exposure_T, trash;

    initialize_opinions(opinions, new_opinions, lmbda, lmbda_T);
    compute_opinions(opinions, N, &initial_opinion, &initial_dev_opinion);

    printf("Initial opinions: %lf +- %lf\n", initial_opinion, initial_dev_opinion);


    // Creating savefile
    FILE *f_results;
    char filename[100];

    sprintf(filename, "results/%s_%d_d=%d_opinions_l=%.2lf_b=%.2lf.txt", network_type, N, HO_degree, lmbda, beta);
    f_results = fopen(filename, "wt");
    create_header(f_results, network_type, HO_degree, initial_opinion, initial_dev_opinion, gamma, lmbda, lmbda_T, beta);
    fprintf(f_results, "# lmbda lmbda_T beta stat steps mean_opinion dev_opinion exposure\n");
    fclose(f_results);


    for(i=0; i<N_stats; i++)
    {
        initialize_opinions(opinions, new_opinions, lmbda, lmbda_T);
        compute_opinions(opinions, N, &mean_opinion, &dev_opinion);
        
        n_steps = 0;

        do
        {
            actualize_network(opinions, new_opinions, lmbda, lmbda_T, beta, delta, delta_T);
            compute_opinions(opinions, N, &mean_opinion, &dev_opinion);

            if(n_steps%1000 == 0 && n_steps != 0){
                printf("\t %lf +- %lf (stat: %d step %d)\n", mean_opinion, dev_opinion, n_stat, n_steps);
            }

            if(total_change < 1E-10) break; // early stopping. If total opinion change (summed for all agents) is smaller than 1E-10, we have read a stable state
            n_steps += 1;

        }while(n_steps < N_max_steps);

        if(dev_opinion>fabs(mean_opinion) && dev_opinion>0.1){ // That is, if the configuration is polarized
            polarized = 1;
            exposure_T = compute_exposure_HO(opinions, N);
        }
        else{
            exposure_T = 0;
        }
        
        f_results = fopen(filename, "at");
        fprintf(f_results, "%lf %lf %lf %d %d %lf %lf %lf\n", lmbda, lmbda_T, beta, n_stat, n_steps, mean_opinion, dev_opinion, exposure_T);
        printf("\t Run %d (l: %.3lf, l_T: %.3lf, b=%.3lf): %lf +- %lf (%d steps)\n", n_stat, lmbda, lmbda_T, beta, mean_opinion, dev_opinion, n_steps);
        fclose(f_results);

        n_stat += 1;
    }     

    free_all(opinions, new_opinions);
    return 0;
}
