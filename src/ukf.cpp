#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .5; // picked arbitray

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;// picked arbitray

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // Not initialized until first process measurement
  is_initialized_ = false;
  
  // Set state dimension
  n_x_ = 5;
  
  // Set augmented dimension
  n_aug_ = 7;
  
  // Define spreading parameter
  lambda_ = 0;
  
  // Matrix to hold sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  // Vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  
  // Start time
  time_us_ = 0;
}

UKF::~UKF() {}

double normalize_radians(double radians) {
  if (radians > M_PI) {
    radians = radians - M_PI*2;
  } else if (radians < -M_PI) {
    radians = radians + M_PI*2;
  } else {
    return radians;
  }
  // Recurse because we could be more than one pi off.
  return normalize_radians(radians);
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  const auto& raw_measurement = meas_package.raw_measurements_;
  const double estimated_velocity = 4;
  const double estimated_yaw = .5;

  if (!is_initialized_) {
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = raw_measurement(0);
      double phi = raw_measurement(1);
      double rhodot = raw_measurement(2);
      
      x_ << rho * cos(phi), rho * sin(phi), estimated_velocity, rhodot * cos(phi), rhodot * sin(phi);
      
      //state covariance matrix
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, std_radphi_, 0,
            0, 0, 0, 0, std_radphi_;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << raw_measurement(0), raw_measurement(1), estimated_velocity, estimated_yaw, 0.0;
      
      //state covariance matrix. values can be tweaked and tuned.
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
    
  }
  
  // Calculate delta_t, store current time for future
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) { 
  lambda_ = 3 - n_aug_;
  VectorXd x_aug_ = VectorXd(7);
  MatrixXd P_aug_ = MatrixXd(7, 7);
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  
  //create augmented covariance matrix
  MatrixXd bottom_right = MatrixXd(2,2);
  bottom_right << std_a_*std_a_, 0,
        0, std_yawdd_*std_yawdd_;
  P_aug_.fill(0.0);
  P_aug_.bottomRightCorner(2, 2) = bottom_right;
  P_aug_.topLeftCorner(5, 5) = P_;
  
  //create square root matrix
  MatrixXd A_aug = P_aug_.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1) = x_aug_ + std::sqrt(lambda_+n_aug_)*A_aug.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - std::sqrt(lambda_+n_aug_)*A_aug.col(i);
  }

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //create vector for predicted state
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_pred += weights_(i) * Xsig_pred_.col(i);
  }
  
  //create covariance matrix for prediction
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    x_diff(3) = normalize_radians(x_diff(3));
    P_pred += weights_(i) * x_diff * x_diff.transpose();
  }

  x_ = x_pred;
  P_ = P_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  //set measurement dimension, lidar can measure px and py
  int n_z = 2;
  MatrixXd R = MatrixXd(2,2);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;

  // Everything below this line is basically
  // generic, except for the radian normalization handling.
  // ------------------------

  //predicted measurement sigma points
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //transform sigma points into measurement space
    VectorXd state_vec = Xsig_pred_.col(i);
    for (int j = 0; j < n_z; j++) {
      Zsig.col(i)(j) = state_vec(j);
    }
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  S += R;
  
  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  for (int i = 0; i < n_z; i++) {
    z(i) = meas_package.raw_measurements_(i);
  }

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = normalize_radians(x_diff(3));
    VectorXd z_diff = Zsig.col(i) - z_pred;
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  // residual
  VectorXd z_diff = z - z_pred;
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  MatrixXd R = MatrixXd(3,3);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);
  
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  double rho = 0;
  double phi = 0;
  double rho_d = 0;
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //transform sigma points into measurement space
    VectorXd state_vec = Xsig_pred_.col(i);
    double px = state_vec(0);
    double py = state_vec(1);
    double v = state_vec(2);
    double yaw = state_vec(3);
    
    rho = sqrt(px*px+py*py);
    phi = atan2(py,px);
    rho_d = (px*cos(yaw)*v+py*sin(yaw)*v) / rho;
    
    Zsig.col(i) << rho,
                   phi,
                   rho_d;
    
    //calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = normalize_radians(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  S += R;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);  
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = normalize_radians(x_diff(3));
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = normalize_radians(z_diff(1));
    Tc += weights_(i) * x_diff * z_diff.transpose();
    
  }
  
  // residual
  VectorXd z_diff = z - z_pred;
  z_diff(1) = normalize_radians(z_diff(1));
    
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
  
}
