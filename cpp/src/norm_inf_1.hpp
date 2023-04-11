#ifndef NORM_INF_1
#define NORM_INF_1

#include "../Eigen/Dense"

double norm_inf_1(Eigen::MatrixXf mat){
    mat = mat.cwiseAbs();
    double res=0;
    for(int i=0;i<mat.rows();i++){
        res += mat.row(i).maxCoeff();
    }


    return res;
};


#endif
