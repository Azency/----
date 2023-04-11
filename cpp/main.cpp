#include <cstddef>
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include "src/generate_date.hpp"
#include "src/proj_inf1ball6.hpp"
#include <time.h>

int main(){

    clock_t t_begin, t_end;

    // t_begin = time(NULL);
    t_begin = clock();

    int m = 5000;
    int n = 1000;
    double alpha = 1e-3;

    fun_res_gen opt = generate_data(m, n, alpha);

    int max_iter = 100;
    double error = 1e-14;


    // MatrixXf a = *opt.ptr_mat;

    // cout<<"nihao"<<opt.tau<<endl;

    fun_return final = proj_inf1ball6(opt.mat, opt.tau, max_iter, error);

    // t_end = time(NULL);
    t_end = clock();

    //输出的矩阵是
    // MatrixXf X = *final.ptr_X;
    // cout<<"输出的X是"<<X.size()<<endl;

    //输出的向量是
    // VectorXf FX = *final.ptr_FX;
    // cout<<"输出的FX"<<FX<<endl;


    cout<<"消耗的总时间"<<(t_end - t_begin)/float(CLOCKS_PER_SEC)<<"s"<<endl;

    cout<<final.X.any()<<endl;

    return 0;


}