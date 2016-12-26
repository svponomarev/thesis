#pragma once
/*!
    \file  methods.h
    \brief file with all 6 implemented methods (declarations)
    Author     : S. V. Ponomarev
    Version    : 1.0
    Date       : 26.12.2016
*/
#include <opencv2/opencv.hpp>

int solveMethod1(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P);
int solveMethod2(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P);
int solveMethod3(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P);
int solveMethod4(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P);
int solveMethod5(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P);
int solveMethod6(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P);
