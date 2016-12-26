/*!
    \file  projective_main.cpp
    \brief main file with testing of implemented methods for estimation of 3D projective transform
    Implementation of Levenberg-Marquardt algorithm taken from https://gist.github.com/rbabich/3539146 (Ron Babich, May 2008).
    (variable-length arrays are removed).
    Theoretical foundation for implemented methods can be found in:
    Z. Zhang, "Estimating Projective Transformation Matrix (Collineation, Homography)", Microsoft Research Technical Report MSR-TR-2010-63, May 2010.
    available at  http://research.microsoft.com/apps/pubs/?id=131928
    Author     : S. V. Ponomarev
    Version    : 1.0
    Date       : 26.12.2016
*/
#include <iostream>
#include "methods.h"
#include "utils.h"

using std::cout;
using std::endl;

typedef int func_t(std::vector<cv::Vec4d>, std::vector<cv::Vec4d>, cv::Mat&);
typedef func_t* pfunc_t;

int main(int argc, char *argv[])
{
    // Generate 3d projective transformation matrix
    cv::Mat transform = cv::Mat::eye(cvSize(4, 4), CV_64F);
    transform.at<double>(0, 1) = transform.at<double>(3, 1) = transform.at<double>(3, 2) = 1.0;
    transform.at<double>(1, 0) = transform.at<double>(2, 3) = 2.0;
    transform.at<double>(0, 3) = transform.at<double>(2, 1) = 3.0;
    cout << "Original transformation matrix:" << endl;
    print_mat_precision(transform);
    cout << "------------------------------------------------" << endl;

    srand(time(NULL));
    // Generate pair of corresponding points w.r.t. generated 3d projective transform
    const int NUMBER_OF_POINTS = 100;
    std::vector<cv::Vec4d> source_points, target_points; // points in homogeneous coordinates
    for (unsigned int i = 0; i < NUMBER_OF_POINTS; i++) {
        double x = rand() % 100;
        double y = rand() % 100;
        double z = rand() % 100;
        cv::Vec4d v(x, y, z, 1.0);
        cv::Mat tmp = transform * cv::Mat(v);
        cv::Vec4d target(tmp.at<double>(0, 0), tmp.at<double>(1, 0), tmp.at<double>(2, 0), tmp.at<double>(3, 0));
        target /= target[3];
        source_points.push_back(v);
        target_points.push_back(target);
   }
    cv::Mat P;

    // Fill vector of functions
    std::vector<pfunc_t> vector_of_functions;
    vector_of_functions.push_back(&solveMethod1);
    vector_of_functions.push_back(&solveMethod2);
    vector_of_functions.push_back(&solveMethod3);
    vector_of_functions.push_back(&solveMethod4);
    vector_of_functions.push_back(&solveMethod5);
    vector_of_functions.push_back(&solveMethod6);

    // Calculate running time and root-mean-square error (RMS) for each method
    for (unsigned int i = 0; i < vector_of_functions.size(); i++) {
        const clock_t begin_time = clock();
        vector_of_functions[i](source_points, target_points, P);
        cout << "Method " << i + 1 << endl;
        cout << "Estimated projective transformation matrix: " << endl;
        print_mat_precision(P);
        cout << "Running time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
        cout.precision(8);
        cout << ", calculated RMS: " << calculate_RMS(source_points, target_points, P) << endl;
        cout << "------------------------------------------------" << endl;
    }
    return 0;
}
