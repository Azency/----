#ifndef PROX_INF5
#define PROX_INF5


#include "../Eigen/Dense"
#include "norm_inf_1.hpp"
using namespace Eigen;

struct fun_res{
    double xn;
    int kk;

};

fun_res prox_inf5(VectorXf b, VectorXf b_star, double y, double r, int k){
    int n = b.size();
    if(k==0){ 
        k = 1;
    }

    int kk;
    if(y<b_star(1))         { kk = 0;}
    else if (y>b_star(n))   { kk = n;}
    else if (y<b_star(k)) {
        for(int i=k; k>0; k--){
            if(y>b_star(k)){ 
                kk=i;
                break;
            }
        }
    }
    else{
        for(int i=k; i<=n; i++){
            if(y<b_star(i)){
                kk = i - 1;
                break;
            }
        }

    }

    double xn;
    if(kk>1){
        xn = b(kk-1) + r*(y-b_star(kk))/(n-kk+1+r);
    }
    else if (kk==1) {
        xn = r*(y-b_star(1))/(n+r);
    }
    else {
        xn = 0;
    }

    fun_res res = {xn,kk};
    return res;
}


#endif

    






