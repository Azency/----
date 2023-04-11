#ifndef GENERATE_DATE
#define GENERATE_DATE

#include "../Eigen/Dense"
#include "norm_inf_1.hpp"
#include <iostream>

using namespace Eigen;
using namespace std;

struct fun_res_gen{
    MatrixXf mat;
    double tau;

};


fun_res_gen generate_data(const int m, const int n, double alpha){
    MatrixXf mat = MatrixXf::Random(m,n) - MatrixXf::Constant(m,n,0.5);

    // cout<<mat<<endl;

    double tau = alpha*norm_inf_1(mat);

    fun_res_gen res = {mat, tau};

    return res;


};

#endif