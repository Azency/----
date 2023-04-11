#include <iostream>
#include "Eigen/Dense"
#include "Eigen/src/Core/Array.h"
#include "Eigen/src/Core/Matrix.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Matrix<int, 3, 3> Matrix3i;

int main()
{
    /*
    Matrix的初始化方法
    Eigen::Matrix<int, 3, 3>
    int 代表Matrix的数据类型，3，3 分别代表 rows， cols
    Matrix3i m1;
    m1(0,0) = 1
    m1(0,1) = 2
    m1(0,2) = 3
    ...

    或者用 m1 << 1,2,3 ...

    */

    Matrix3i m1;
    m1 << -1, 2, -3, 4, 5, 6, 7, 8, 9;
    cout << "m1 = \n" << m1 << endl;

    Matrix3i m2;
    m2 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    cout << "m2 = \n" << m2 << endl;

    MatrixXf mat = MatrixXf:: Random(3,3);

    Eigen::ArrayXXf arr = mat.array().abs();
    cout<<mat<<endl;

    Eigen::ArrayXXf app = arr.row(0);

    double p = 0;
    for(int i=0;i<mat.cols();i++){
        p += mat.cwiseAbs().row(i).maxCoeff();
    }

    // auto p = ;
    mat = mat.cwiseAbs();

    VectorXf a = mat.col(1);
    mat.col(0) = a + VectorXf::Constant(3,100);
    cout << "m1的i-无穷范数" <<mat<< endl;

    return 0;
}
