#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>

//double inner_control_mesh(int i, int j, int b, double r, double delta_t, double m, std::vector< std::vector< std::vector<double> > >& W, std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector<double> >& V)

double* three_dim_index(double* matrix, int i, int j, int k, double m, int b);

double* two_dim_index(double* vector, int i, int j, double m, int b);

double inner_control_mesh(int i, int j, int b, double r, double delta_t, double m, double* W, double* X, double* V){

int m_int= (int)m;
double ControlMeshContinVal=0;

double sum=0, ContinVal=0, stock_expectation=0, true_stock_expectation=0, numerator=0, denominator=0, beta=0;
   
for(int k=0; k<b; k++){

	//sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]                         
	sum+= (*three_dim_index(W, (m_int-i), k, j, m, b)) * (*two_dim_index(V, (i-1), k, m, b));
}

ContinVal=(1/((double)b))*sum; 

sum=0;
   
for(int l=0; l<b; l++){
                   
        //sum+=(W[(m-i)][l][j])*exp(X[m-i][l][0]); //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]                         
	sum+=(*three_dim_index(W, (m_int-i), l, j, m, b)) * exp((*three_dim_index(X, (m_int-i), l, 0, m, b)));
}

stock_expectation=(1/((double)b))*sum; 

//true_stock_expectation=exp(X[m-i-1][j][0])*exp(r*delta_t);
true_stock_expectation=exp(*three_dim_index(X, (m_int-i-1), j, 0, m, b)) * exp(r*delta_t);

for(int p=0; p<b; p++){

//numerator += ( W[(m-i)][p][j] * exp( X[m-i][p][0] ) - stock_expectation ) * ( W[(m-i)][p][j] * V[i-1][p] - ContinVal );
numerator += ( (*three_dim_index(W, (m_int-i), p, j, m, b)) * exp( *three_dim_index(X, (m_int-i), p, 0, m, b) ) - stock_expectation ) * ( (*three_dim_index(W, (m_int-i), p, j, m, b)) * (*two_dim_index(V, (i-1), p, m, b)) - ContinVal );

//denominator += pow( ( W[(m-i)][p][j]  * exp( X[m-i][p][0] ) - stock_expectation ) , 2 );
denominator += pow( ( (*three_dim_index(W, (m_int-i), p, j, m, b))  * exp( *three_dim_index(X, (m_int-i), p, 0, m, b) ) - stock_expectation ) , 2 );

}

beta=numerator/denominator;

ControlMeshContinVal= ContinVal-beta*(stock_expectation-true_stock_expectation);

return ControlMeshContinVal;
}
