/*!
    \file  utils.cpp
    \brief file with support functions
    Author     : S. V. Ponomarev
    Version    : 1.0
    Date       : 26.12.2016
*/
#include "utils.h"
#include <iomanip>

using std::cout;
using std::endl;

/*!
 * @fn double getNorm2(cv::Vec4d v1, cv::Vec4d v2)
 *
 * calculating Euclidian distance for two vectors (3d points in homogeneous coordinates)
 *
 * @param v1     - given first vector
 * @param v2     - given second vector
 *
 * @return double : Euclidian distance (norm2)
 *
*/
double getNorm2(cv::Vec4d v1, cv::Vec4d v2) {
    return sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1]) + (v1[2] - v2[2]) * (v1[2] - v2[2]));
}

/*!
 * @fn double calculate_RMS(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat P)
 *
 * calculating root-mean-square error for two sets of points using given projective transformation matrix P
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second set of points
 * @param P     - given 3d projective transformation matrix (4x4)
 *
 * @return double : root-mean-square error if input arguments are correct, -1 otherwise
 *
*/
double calculate_RMS(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if ((P.rows != P.cols) || (P.rows != 4)) {
        cout << "Error: projective matrix must be 4x4 size! Returning..." << endl;
        return -1;
    }
    double RMS = 0.;
    for (unsigned int i = 0; i < source_points.size(); i++) {
        cv::Mat tmp = P * (cv::Mat(source_points[i]));
        cv::Vec4d transformed(tmp.at<double>(0, 0), tmp.at<double>(1, 0), tmp.at<double>(2, 0), tmp.at<double>(3, 0));
        transformed /= transformed[3];
        RMS += getNorm2(transformed, target_points[i]);
    }
    RMS /= source_points.size();
    return RMS;
}

/*!
 * @fn void print_points(cv::Vec4d v1, cv::Vec4d v2)
 *
 * printing corresponding pair of points
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 *
 *
*/
void print_points(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points)
{
    for (unsigned int i = 0; i < source_points.size(); i++)
    {
        cout << "(" << source_points[i][0] << "," << source_points[i][1] << "," << source_points[i][2] << "," << source_points[i][3] << ")";
        cout << " - (" << target_points[i][0] << "," << target_points[i][1] << "," << target_points[i][2] << "," << target_points[i][3] << ")" << endl;
    }
}

/*!
 * @fn void print_mat_precision(cv::Mat mat, int prec)
 *
 * printing openCV mat object with specific precision
 *
 * @param mat     - given openCV matrix (double, type should be CV_64F)
 * @param prec    - given precision (number of decimal places, 4 by default)
 *
 *
*/
void print_mat_precision(cv::Mat mat, int prec)
{
    std::cout << std::fixed;
    for(int i = 0; i < mat.rows; i++)
    {
        cout << "[";
        for(int j = 0; j < mat.cols; j++)
        {
            cout << std::setprecision(prec) << mat.at<double>(i,j);
            if(j != mat.rows - 1)
                cout << ", ";
            else
                cout << "]" << endl;
        }
    }
}

/*!
 * @fn void print_vec_precision(cv::Vec4d v, int prec)
 *
 * printing openCV Vec4d object with specific precision
 *
 * @param mat     - given openCV vector
 * @param prec    - given precision (number of decimal places, 4 by default)
 *
 *
*/
void print_vec_precision(cv::Vec4d v, int prec) {
    cout << std::fixed;
    for (int i = 0; i < v.rows; i++)
        cout << std::setprecision(prec) << v[i] << " ";
    cout << endl;
}
