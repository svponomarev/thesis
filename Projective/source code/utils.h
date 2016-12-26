#pragma once
/*!
    \file  utils.h
    \brief file with support functions (declaration)
    Author     : S. V. Ponomarev
    Version    : 1.0
    Date       : 26.12.2016
*/
#include <opencv2/opencv.hpp>

double getNorm2(cv::Vec4d v1, cv::Vec4d v2);
double calculate_RMS(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat P);

void print_points(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points);
void print_mat_precision(cv::Mat mat, int prec = 4);
void print_vec_precision(cv::Vec4d v, int prec = 4);

