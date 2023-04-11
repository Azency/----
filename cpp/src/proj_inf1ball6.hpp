#ifndef PROJ_INF1BALL6
#define PROJ_INF1BALL6

#include "../Eigen/Dense"
#include "norm_inf_1.hpp"
#include "prox_inf5.hpp"
#include <cmath>


using namespace Eigen;

struct  fun_return{
    MatrixXf X;
    VectorXf FX;
};

void sort_vec(const VectorXf& vec, VectorXf& sorted_vec,  VectorXi& ind){
  ind=VectorXi::LinSpaced(vec.size(),0,vec.size()-1);//[0 1 2 3 ... N-1]
  auto rule=[vec](int i, int j)->bool{
    return vec(i)>vec(j);
  };//正则表达式,作为sort的谓词
  std::sort(ind.data(),ind.data()+ind.size(),rule);
  //data成员函数返回VectorXd的第一个元素的指针，类似于begin()
  sorted_vec.resize(vec.size());
  for(int i=0;i<vec.size();i++){
    sorted_vec(i)=vec(ind(i));
  }
};

fun_return proj_inf1ball6(MatrixXf B, double tau, int max_iter=100, double delt=10^(-8)){
// % project the matrix B onto the l-1-inf ball of radius tau, i.e.
// % X = argmin |X-B|_F
// % s.t. sum(max(abs(B), [] ,2)) <= tau
// % max_iter=100, delt=1e-8 by default
    MatrixXf X;
    VectorXf FX;

    double eps = 2.2204e-16;

    // if (nargin == 2){
    //     max_iter = 100;
    //     delt = 10^(-8);

    // }
    
    if(norm_inf_1(B) - tau < eps){
        X = B;
        FX = VectorXf::Zero(1);
        fun_return res = {X,FX};
        return res;
    }

    int m = B.rows();
    int n = B.cols();

    VectorXf x_inf = VectorXf::Zero(m);
    VectorXf kk = VectorXf::Zero(m);
    int r = n;
    double mu = 0;
    VectorXf z = VectorXf::Zero(m);
    FX = VectorXf::Zero(m);
    VectorXf G = VectorXf::Ones(m);
    MatrixXf sB = MatrixXf::Identity(m,n);
    MatrixXf aB = B.cwiseAbs();
    VectorXf saB = aB.rowwise().sum();
    MatrixXf B2 = MatrixXf::Zero(m, n);
    MatrixXf B3 = MatrixXf::Zero(m, n);
    MatrixXf aBt = aB.transpose();
    MatrixXf B2t = B2.transpose();


    VectorXi ind;
    VectorXf origin;
    VectorXf sorted;
    for(int ii=0;ii<m;ii++){
        origin = aBt.col(ii);
        sort_vec(origin,sorted, ind);
        B2t.col(ii) = sorted;
    }
    B2 = B2t.transpose();

    B3.col(0) = -saB;
    saB = saB - B2.col(0);

    for(int k=1;k<n;k++){
        B3.col(k) = -saB + (n+1-k+r)*B2.col(k-1);
        saB = saB - B2.col(k);
    }
    B3 = B3/r;



    B2t = B2.transpose();
    MatrixXf B3t = B3.transpose();

    int iter_;
    for(int iter=0;iter<max_iter;iter++){
        iter_ = iter;
        for(int ii=0;ii<m;ii++){
            if(G(ii)){
                fun_res temp = prox_inf5(B2t.col(ii), B3t.col(ii), z(ii)-mu/r, r, kk(ii));
                x_inf(ii) = temp.xn;
                kk(ii) = temp.kk;
            }
        }

        MatrixXf f = x_inf - MatrixXf:: Constant(x_inf.rows(),x_inf.cols(),eps);
        double m_prime = f.cwiseSign().sum();

        double temp = x_inf.sum() - tau;
        mu = mu + r/m_prime*temp;
        // % update z   
        z = x_inf - MatrixXf:: Constant(x_inf.rows(), x_inf.cols() ,1/m_prime*temp);

        FX(iter) = std:: abs(temp);

        if(iter>1){
            double eFX = 1e8*FX(iter-1) - 1e8*FX(iter);
            if(FX(iter) < delt || eFX<delt){
                break;
            }
        }
    }

    FX = FX.head(iter_);
    X = aB;
    for(int ii=0;ii<m;ii++){
        for(int j=0;j<X.cols();j++){
            if(X(ii,j)>x_inf(ii)){ X(ii,j) = x_inf(ii);}
        }
    }

    X = sB.cwiseProduct(X);

    fun_return res = {X, FX};

    return res;


}




#endif