/*!
    \file  methods.cpp
    \brief file with all 6 implemented methods
    Implementation of Levenberg-Marquardt algorithm taken from https://gist.github.com/rbabich/3539146 (Ron Babich, May 2008).
    (variable-length arrays are removed).
    Theoretical foundation for implemented methods can be found in:
    Z. Zhang, "Estimating Projective Transformation Matrix (Collineation, Homography)", Microsoft Research Technical Report MSR-TR-2010-63, May 2010.
    available at  http://research.microsoft.com/apps/pubs/?id=131928
    Author     : S. V. Ponomarev
    Version    : 1.0
    Date       : 26.12.2016
*/
#include "methods.h"
#include "utils.h"
#include "levmarq.h"
#include <ctime>
using std::cout;
using std::endl;

/*!
 * @fn void method1_calculateA(std::vector<cv::Vec4d> points, cv::Mat &A)
 *
 * calculating collineation A which maps the standard reference points ([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,1,1,1]) to the given set of points
 *
 * @param points     - given set of points
 * @param A     - output calculated collineation
 *
*/
void method1_calculateA(std::vector<cv::Vec4d> points, cv::Mat &A) {
    cv::Mat M = cv::Mat::eye(cvSize(4, 4), CV_64F);
    for (unsigned int i = 0; i != points.size() - 1; i++)
        for (unsigned int j = 0; j != 4; j++)
            M.at<double>(j, i) = points[i][j];
    cv::Mat tmp = M.inv() * cv::Mat(points[4]);
    cv::Vec4d lambd(tmp.at<double>(0, 0), tmp.at<double>(1, 0), tmp.at<double>(2, 0), tmp.at<double>(3, 0)); // vector of scalar factors
    A = cv::Mat::eye(cvSize(4, 4), CV_64F);
     for (unsigned int i = 0; i != points.size() - 1; i++)
        for (unsigned int j = 0; j != 4; j++)
            A.at<double>(j, i) = points[i][j] * lambd[i];
}

/*!
 * @fn int solveMethod1(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P)
 *
 * calculating projective transform P using strictly minimum number of corresponding points (5 pairs). Simple, but unstable method.
 * if size of point sets more than 5, calculation is performed with 5 random chosen pairs.
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 * @param P     - output calculated projective transformation matrix
 *
 * @return int : 0 if projective matrix was successfully calculated, -1 otherwise
 *
*/
int solveMethod1(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if (source_points.size() < 5) {
        cout << "Error: number of points less than 5! Returning..." << endl;
        return -1;
    }
    cv::Mat A1, A2; // two matrices for collineation to reference points for two sets of points
    std::srand ( unsigned ( std::time(0) ) );
    if (source_points.size() > 5) { // use only 5 random pair of points from all given
        std::vector<cv::Vec4d> source_points5, target_points5;
        std::vector<int> indexes(source_points.size());
        for (unsigned int i = 0; i < source_points.size(); i++)
            indexes[i] = i;
        random_shuffle(indexes.begin(), indexes.end()); // randomize order of indexes
        for (unsigned int i = 0; i < 5; i++) {
            source_points5.push_back(source_points[indexes[i]]);
            target_points5.push_back(target_points[indexes[i]]);
        }
         method1_calculateA(source_points5, A1);
         method1_calculateA(target_points5, A2);
    }
    else {
        method1_calculateA(source_points, A1);
        method1_calculateA(target_points, A2);
    }
    P = A2 * A1.inv();
    P /= P.at<double>(3, 3); // normalize projective matrix so P[3][3] = 1
    return 0;
}

/*!
 * @fn int solveMethod2(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P)
 *
 * calculating projective transform P together with the scalar factors (non weighted least-squares solution).
 * The most precise method, but very slow for big size of point sets (1000 and more).
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 * @param P     - output calculated projective transformation matrix
 *
 * @return int : 0 if projective matrix was successfully calculated, -1 otherwise
 *
*/
int solveMethod2(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if (source_points.size() < 5) {
        cout << "Error: number of points less than 5! Returning..." << endl;
        return -1;
    }
    int rows = source_points.size() * 4;
    int cols = source_points.size() + 15;
    cv::Mat A = cv::Mat::zeros(cvSize(cols, rows), CV_64F);
    for (int i = 0; i < A.rows; i++) {
        double* Ai = A.ptr<double>(i);
        int pt_nmbr = i / 4;
        int shift = i % 4;
        for(int j = 0; j < A.cols; j++)
            if (j >= shift * 4 && j < shift * 4 + 4)
                Ai[j] = source_points[pt_nmbr][j - shift*4];
        if (pt_nmbr != 0)
            Ai[15 + pt_nmbr] = -target_points[pt_nmbr][shift];
    }
    cv::Mat b = cv::Mat::zeros(cvSize(1, rows), CV_64F);
    for (unsigned int i = 0; i != 4; i++) {
        double* bi = b.ptr<double>(i);
        bi[0] = target_points[0][i];
    }
    cv::Mat At = A.t();
    cv::Mat inv;
    cv::invert(At * A, inv);

    cv::Mat XM = cv::Mat::zeros(cvSize(1, cols), CV_64F);
    XM = inv * At * b; // resulting vector of 15 + n parameters
    P = cv::Mat::eye(cvSize(4, 4), CV_64F);
    for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
            P.at<double>(y, x) = XM.at<double>(y * 4 + x, 0);
    P /= P.at<double>(3, 3); // normalize projective matrix so P[3][3] = 1
    return 0;
}

/*!
 * @fn int solveMethod3(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P)
 *
 * calculating only projective transform P without the scalar factors using batch approach (calculation of eigenvectors).
 * non-weighted error function implementation
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 * @param P     - output calculated projective transformation matrix
 *
 * @return int : 0 if projective matrix was successfully calculated, -1 otherwise
 *
*/
int solveMethod3(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if (source_points.size() < 5) {
        cout << "Error: number of points less than 5! Returning..." << endl;
        return -1;
    }
    unsigned int n = source_points.size();
    cv::Mat B = cv::Mat::zeros(cvSize(16, 3), CV_64F);
    cv::Mat A = cv::Mat::zeros(cvSize(16, 16), CV_64F);
    for (unsigned int i = 0; i != n; i++) {
        for (int j = 0; j < 4; j++)
            B.at<double>(0, j) = source_points[i][j] * target_points[i][1];
        for (int j = 0; j < 4; j++)
            B.at<double>(0, j + 4) = source_points[i][j] * (-1) * target_points[i][0];
        for (int j = 0; j < 4; j++)
            B.at<double>(1, j) = source_points[i][j] * target_points[i][2];
        for (int j = 0; j < 4; j++)
            B.at<double>(1, j + 8) = source_points[i][j] * (-1) * target_points[i][0];
        for (int j = 0; j < 4; j++)
            B.at<double>(2, j) = source_points[i][j] * target_points[i][3];
        for (int j = 0; j < 4; j++)
            B.at<double>(2, j + 12) = source_points[i][j]  * (-1) * target_points[i][0];
        cv::Mat Bt = B.t();
        cv::Mat pr = Bt * B;
        A += pr;
    }
    cv::Mat eigenvalues, eigenvectors;
    cv::eigen(A, eigenvalues, eigenvectors);
    P = cv::Mat::eye(cvSize(4, 4), CV_64F);
    for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
            P.at<double>(y, x) = eigenvectors.at<double>(15, y * 4 + x); // the solution is the eigenvector corresponding to the smallest eigenvalue (last row)
    P /= P.at<double>(3, 3); // normalize projective matrix so P[3][3] = 1
    return 0;
}

/*!
 * @fn int solveMethod4(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P)
 *
 * calculating only projective transform P using iterative method based on the Kalman filtering technique.
 * Currently gives the worst RMS error w.r.t. other methods, result strongly depends on values in the covariation matrices (processNoiseCov, errorCovPost, measurementNoiseCov)
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 * @param P     - output calculated projective transformation matrix
 *
 * @return int : 0 if projective matrix was successfully calculated, -1 otherwise
 *
*/
int solveMethod4(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if (source_points.size() < 5) {
        cout << "Error: number of points less than 5! Returning..." << endl;
        return -1;
    }
    int n = source_points.size();
    cv::Mat P0;
    solveMethod1(source_points, target_points, P0); // Initial estimation of P based on 5 pair of points
    cv::Mat X = cv::Mat::zeros(cvSize(1, 15), CV_64F);
    for(int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if ((y == 3) && (x == 3))
                break;
            X.at<double>(y * 4 + x, 0) = P0.at<double>(y,x); // convert matrix to vector
        }
    }
    P = cv::Mat::eye(cvSize(4, 4), CV_64F);
    // Kalman filter initialization
    cv::KalmanFilter KF(15, 3, 0);
    setIdentity(KF.transitionMatrix); // Transition matrix = identity matrix
    setIdentity(KF.processNoiseCov, cv::Scalar::all(1e6)); // result is better with big values
    setIdentity(KF.errorCovPost, cv::Scalar::all(1e9)); // careful with overflow
    for(int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if ((y == 3) && (x == 3))
                break;
            KF.statePost.at<float>(y * 4 + x) = P0.at<double>(y,x); // KF.statePost = P0
        }
    }
    for (int pt = 5; pt < n; pt++) { // loop over the rest n - 5 points
        cv::Mat_<float> measurement(3,1); measurement.setTo(cv::Scalar(0));
        measurement.at<float>(0, 2) = target_points[pt][0] * source_points[pt][3]; // y_i = [0, 0, m'1*m4]
        // Measurement matrix = M_i
        for (int j = 0; j < 4; j++)
            KF.measurementMatrix.at<float>(0, j) = source_points[pt][j] * target_points[pt][1];
        for (int j = 0; j < 4; j++)
            KF.measurementMatrix.at<float>(0, j + 4) = source_points[pt][j] * (-1) * target_points[pt][0];
        for (int j = 0; j < 4; j++)
            KF.measurementMatrix.at<float>(1, j) = source_points[pt][j] * target_points[pt][2];
        for (int j = 0; j < 4; j++)
            KF.measurementMatrix.at<float>(1, j + 8) = source_points[pt][j] * (-1) * target_points[pt][0];
        for (int j = 0; j < 4; j++)
            KF.measurementMatrix.at<float>(2, j) = source_points[pt][j] * target_points[pt][3];
        for (int j = 0; j < 3; j++)
            KF.measurementMatrix.at<float>(2, j + 12) = source_points[pt][j]  * (-1) * target_points[pt][0];
        // measurementNoiseCov = W_i
        cv::Mat tmp = P0 * cv::Mat(source_points[pt]);
        cv::Vec4d v(tmp.at<double>(0, 0), tmp.at<double>(1, 0), tmp.at<double>(2, 0), tmp.at<double>(3, 0));
        cv::Mat Covariance = cv::Mat::eye(cvSize(4, 4), CV_32F); // Covariance matrix = Identity, if the uncertainty can be considered as identical and isotropic for all points
        cv::Mat Jacobian1 = cv::Mat::zeros(cvSize(4, 3), CV_32F);
        for (int y = 0; y < Jacobian1.rows; y++)
            for (int x = 0; x < Jacobian1.cols; x++)
                Jacobian1.at<float>(y, x) = target_points[pt][y + 1] * P0.at<double>(0,x) - target_points[pt][0] * P0.at<double>(y + 1,x);
        cv::Mat Jacobian2 = cv::Mat::zeros(cvSize(4, 3), CV_32F);
        Jacobian2.at<float>(0, 1) = Jacobian2.at<float>(1, 2) = Jacobian2.at<float>(2, 3) = v[0];
        for (int y = 0; y < 3; y++)
            Jacobian2.at<float>(y, 0) = -v[y + 1];
        KF.measurementNoiseCov = Jacobian1 * Covariance * Jacobian1.t() + Jacobian2 * Covariance * Jacobian2.t();

        KF.predict(); // first step - prediction
        cv::Mat estimated = KF.correct(measurement); // second step - correction
        for (int k = 0; k < estimated.rows; k++) {
            int row = k / 4;
            int col = k % 4;
            P.at<double>(row, col) = estimated.at<float>(k, 0); // convert vector to matrix
        }
    }
    return 0;
}

/*!
 * @fn int solveMethod5(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P)
 *
 * calculating projective transform P through normalization, point vectors and projective matrix should be all normalized, so norm = 1
 * the best solution is the eigenvector with the smallest eigenvalue.
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 * @param P     - output calculated projective transformation matrix
 *
 * @return int : 0 if projective matrix was successfully calculated, -1 otherwise
 *
*/
int solveMethod5(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if (source_points.size() < 5) {
        cout << "Error: number of points less than 5! Returning..." << endl;
        return -1;
    }
    unsigned int n = source_points.size();
     for (unsigned int i = 0; i != n; i++) {
        cv::normalize(source_points[i], source_points[i]); // points normalization
        cv::normalize(target_points[i], target_points[i]);
     }
    cv::Mat A = cv::Mat::zeros(cvSize(16, 16), CV_64F);
    for (unsigned int i = 0; i != n; i++) {
        cv::Mat M1 = cv::Mat::zeros(cvSize(4, 16), CV_64F);
        cv::Mat M2 = cv::Mat::zeros(cvSize(16, 4), CV_64F);
        cv::Mat Mi1 = cv::Mat::zeros(cvSize(1, 4), CV_64F);
        cv::Mat I = cv::Mat::eye(cvSize(4, 4), CV_64F);
        for (unsigned int j = 0; j < 4; j++) {
            Mi1.at<double>(j, 0) = target_points[i][j];
        }
        normalize(Mi1, Mi1); // matrix normalization
        cv::Mat MM = cv::Mat::zeros(cvSize(4, 4), CV_64F);
        for (int j = 0; j < M1.rows; j++) {
            int shift = j / 4;
            M1.at<double>(j, shift) = source_points[i][j % 4];
        }
        MM = I - Mi1 * Mi1.t();
        for (int y = 0; y < M2.rows; y++) {
            int shift = y % 4;
            for (int x = 0; x < M2.cols; x++) {
                if (x >= shift * 4 && x < shift * 4 + 4)
                    M2.at<double>(y, x) = source_points[i][x - shift*4];
            }
        }
        cv::Mat pr = M1 * MM * M2;
        A += pr; // A contains sum for all points
    }
    cv::Mat eigenvalues, eigenvectors;
    cv::eigen(A, eigenvalues, eigenvectors);
    P = cv::Mat::eye(cvSize(4, 4), CV_64F);
    for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
            P.at<double>(y, x) = eigenvectors.at<double>(15, y * 4 + x); // smallest eigenvalue = last row
    P /= P.at<double>(3, 3); // normalize projective matrix so P[3][3] = 1
    return 0;
}

/*!
 * @fn double method6_func(double *par, int x, void *fdata)
 *
 * calculation of objective function for estimation each individual ideal corresponding point (how good transformation describes correspondence for specific pair of points)
 *
 * @param par   - given vector of parameters (coefficients of projective transform)
 * @param x     - number of measurement (number of pair of corresponding points)
 * @param fdata - vector with all pairs of corresponding points
 *
 * @return double : value of function for given pair of corresponding points
 *
*/
double method6_func(double *par, int x, void *fdata)
{
    cv::Mat P = cv::Mat::zeros(cvSize(4,4), CV_64F);
    for (int i = 0; i < 16; i++) {
        int row = i / 4;
        int col = i % 4;
        P.at<double>(row, col) = par[i];
    }
    cv::Mat mi1x = cv::Mat(cvSize(1, 4), CV_64F); // point from first set in homogeneous coordinates (x - eXtended)
    cv::Mat mi2x = cv::Mat(cvSize(1, 4), CV_64F); // point from second set in homogeneous coordinates (x - eXtended)
    for (int i = 0; i < 4; i++) {
        mi1x.at<double>(i, 0) = ((double *)fdata)[x*8 + i];
        mi2x.at<double>(i, 0) = ((double *)fdata)[x*8 + 4 + i];
    }
    cv::Mat mci1x = mi1x + P.inv()*mi2x; // ideal point for point in first set (estimation) in homogeneous coordinates
    mci1x /= mci1x.at<double>(3, 0);
    cv::Mat mci2x = P * mci1x; // corresponding point in second set for ideal point for point in first set (estimation) in homogeneous coordinates
    mci2x /= mci2x.at<double>(3, 0);

    cv::Mat mci1, mci2; // ideal points in cartesian coordinates (x,y,z)
    mci1x(cv::Rect(0,0,1,3)).copyTo(mci1);
    mci2x(cv::Rect(0,0,1,3)).copyTo(mci2);

    cv::Mat P3x3, P4, P43; // specific parts of projective matrix P(4x4)
    P(cv::Rect(0,0,3,3)).copyTo(P3x3);
    P(cv::Rect(0,3,3,1)).copyTo(P43);
    P(cv::Rect(0,3,4,1)).copyTo(P4);

    cv::Mat division = P4 * mci1x;
    cv::Mat Jacobian = P3x3 - mci2 * P43;
    Jacobian = Jacobian / division;

    cv::Mat mi1, mi2; // points in cartesian coordinates (x,y,z)
    mi1x(cv::Rect(0,0,1,3)).copyTo(mi1);
    mi2x(cv::Rect(0,0,1,3)).copyTo(mi2);

    cv::Mat Covariance = cv::Mat::eye(cvSize(3, 3), CV_64F); // Covariance matrix = Identity, if the uncertainty can be considered as identical and isotropic for all points
    cv::Mat W = Covariance.inv()  + Jacobian.t() * Covariance.inv() * Jacobian;
    cv::Mat new_mci1 = W.inv() * (Covariance.inv() * mi1 + Jacobian.t() * Covariance.inv() * (mi2 - mci2 + Jacobian * mci1) ); // new estimation of ideal point for point in first set
    cv::Mat diff1 = W.inv() * (Jacobian.t()*Covariance.inv()*Jacobian*(mi1 - mci1) - Jacobian.t()*Covariance.inv()*(mi2 - mci2));
    cv::Mat diff2 = W.inv() * (Covariance.inv()*(mi1 - mci1) + Jacobian.t()*Covariance.inv()*(mi2 - mci2));
    cv::Mat equation = mi2 - mci2 - Jacobian * diff2;
    cv::Mat Fmat = diff1.t() * Covariance.inv() * diff1 + equation.t() * Covariance.inv() * equation;

    double F = Fmat.at<double>(0,0);
    return F;
}

/*!
 * @fn void method6_gradient(double *g, double *par, int x, void *fdata)
 *
 * numerical computation of function's gradient for method 6 and Levenberg-Marquardt algorithm
 *
 * @param g     - output gradient vector (size of vector = number of parameters)
 * @param par   - given vector of parameters (coefficients of projective transform)
 * @param x     - number of measurement (number of pair of corresponding points)
 * @param fdata - vector with all pairs of corresponding points
 *
*/
void method6_gradient(double *g, double *par, int x, void *fdata)
{
    double par0[16];
    for (int i = 0; i < 16; i++)
        par0[i] = par[i]; // save initial values of parameters
    double step = 1e-3;
    for (int i = 0; i < 16; i++) {
        par[i] = par0[i] + step;
        double f2 = method6_func(par, x, fdata);
        par[i] = par0[i] - step;
        double f1 = method6_func(par, x, fdata);
        par[i] = par0[i];
        g[i] = f2 - f1;
    }
}

/*!
 * @fn int solveMethod6(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P)
 *
 * calculating projective transform P through maximum likelihood estimation. nonlinear minimization of P is performed with the Levenberg-Marquardt algorithm.
 * slow iterative approach.
 *
 * @param source_points     - given first set of points
 * @param target_points     - given second pair of points
 * @param P     - output calculated projective transformation matrix
 *
 * @return int : 0 if projective matrix was successfully calculated, -1 otherwise
 *
*/
int solveMethod6(std::vector<cv::Vec4d> source_points, std::vector<cv::Vec4d> target_points, cv::Mat &P) {
    if (source_points.size() != target_points.size()) {
        cout << "Error: number of points doesn't match! Returning..." << endl;
        return -1;
    }
    if (source_points.size() < 5) {
        cout << "Error: number of points less than 5! Returning..." << endl;
        return -1;
    }
    unsigned int n = source_points.size();
    cv::Mat P0;
    solveMethod1(source_points, target_points, P0); // Initial estimation of P based on 5 pair of points
    // Filling parameters for levmarq algorithm
    double *f_data = new double[8 * source_points.size()]; // function data (coords of corresponding points)
    int k = 0;
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < 4; j++)
            f_data[k++] = source_points[i][j];
        for (unsigned int j = 0; j < 4; j++)
            f_data[k++] = target_points[i][j];
    }
    k = 0;
    double params[16]; // params - vector with coefficients of projective transformation matrix
    for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
            params[k++] = P0.at<double>(y, x);

    const int N_MEASUREMENTS = source_points.size();
    const int N_PARAMS = 16;

    double *t_data = new double[N_MEASUREMENTS]; // vector with measurements of function f.
    for (int i = 0; i < N_MEASUREMENTS; i++)
        t_data[i] = 0.; // The desired values is zeros, because f shows errors of correspondence between points using current P

    LMstat lmstat;
    levmarq_init(&lmstat);
    levmarq(N_PARAMS, params, N_MEASUREMENTS, t_data, NULL, &method6_func, &method6_gradient, f_data, &lmstat); // Levenberg-Marquardt algorithm
    delete [] f_data;
    delete [] t_data;

    cv::Mat P_res = cv::Mat::zeros(cvSize(4, 4), CV_64F);
    P = cv::Mat::eye(cvSize(4, 4), CV_64F);
    for (int i = 0; i < 16; i++) {
        int rows = i / 4;
        int cols = i % 4;
        P.at<double>(rows, cols) = params[i]; // convert vector to matrix
    }
    P /= P.at<double>(3, 3); // normalize projective matrix so P[3][3] = 1
    return 0;
}

