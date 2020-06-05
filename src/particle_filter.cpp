/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random> // Need this for sampling from distributions
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // Set the number of particles
  
  // Declare the stand deviation for [standard deviation of x [m], 
  // standard deviation of y [m], standard deviation of yaw [rad]]
  std::default_random_engine gen;
  const double std_x = std[0];
  const double std_y = std[1];
  const double std_theta = std[2];

  // Creates a normal (Gaussian) distribution for x, y, and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  Particle particle;

  for(int i=0; i<num_particles; i++){
    // Initialize all parameter
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0f;

    particles.push_back(particle);
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  //const double std_x = std_pos[0];
  //const double std_y = std_pos[1];
  //const double std_theta = std_pos[2];

  for(int i=0; i<num_particles; i++){

    double theta_f = particles[i].theta + yaw_rate*delta_t;
    particles[i].x += velocity/yaw_rate*(sin(theta_f)-sin(particles[i].theta));
    particles[i].y += velocity/yaw_rate*(cos(particles[i].theta)-cos(theta_f));
    particles[i].theta = theta_f;

  }



}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     const Map &map_landmarks) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for(int i=0; i<predicted.size(); i++) {

    int dist_min_ind;
    double dist_j;
    double dist_min = dist(predicted[i].x, predicted[i].y, 
                           map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);
    

    for(int j=0; j<map_landmarks.landmark_list.size(); j++) {
      dist_j = dist(predicted[i].x, predicted[i].y, 
                    map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
      if(dist_j < dist_min) {
        dist_min = dist_j;
        dist_min_ind = j;
      }
    }

    predicted[i].id = dist_min_ind;
    
  }

}

vector<LandmarkObs> ParticleFilter::VehicleCoor2MapCoor(const vector<LandmarkObs> &observations) {
  /**
   * VehicleCoor2MapCoor convert the coordinate from VEHICLE to MAP
   * @param particles particles, in VEHICLE coordinate
   * @param observations observations, in VEHICLE coordinate
   */

  vector<LandmarkObs> observations_map;
  const double theta = -M_PI/2; // -90 degrees

  for(int i=0; i<observations.size(); i++){
    observations_map[i].x = particles[i].x + (cos(theta)*observations[i].x) - (sin(theta)*observations[i].y);
    observations_map[i].y = particles[i].y + (sin(theta)*observations[i].x) + (cos(theta)*observations[i].y);
  }

  return observations_map;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */


  vector<LandmarkObs> observations_map;
  observations_map = VehicleCoor2MapCoor(observations); //in MAP coordiate

  dataAssociation(observations_map, map_landmarks);

  //Define the observations standard deviation
  const double sigma_x = std_landmark[0];
  const double sigma_y = std_landmark[1];

  for(int i=0; i<observations_map.size(); i++){
        double x = observations_map[i].x;
        double y = observations_map[i].y;

        double mu_x = map_landmarks.landmark_list[0].x_f;
        double mu_y = map_landmarks.landmark_list[0].y_f;

        //Apply Multivariate-Gaussian probability density
        double exponent = exp(-(pow((x-mu_x),2)/(2*pow(sigma_x,2))+
                                pow((y-mu_y),2)/(2*pow(sigma_x,2))));
        
        particles[i].weight = 1/(2*M_PI*sigma_x*sigma_y) * exponent;
  }


}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::default_random_engine gen;
  std::uniform_int_distribution<int> distr(0,num_particles-1);
  
  vector<Particle> particles_3;
  int index = distr(gen); //need to redo
  double beta = 0.0f;

  double highest_weight = -1.0;
  for (int i = 0; i < num_particles; ++i) {
    if (particles[i].weight > highest_weight) {
      highest_weight = particles[i].weight;
    }
  }
  
  for (int i = 0; i < num_particles; ++i) {
    beta += distr(gen)/num_particles * 2 * highest_weight; //need to redo
    while(beta > particles[index].weight){
      beta -= particles[index].weight;
      index = (index+1) % num_particles;
    }
    particles_3.push_back(particles[index]);
  }

  particles = particles_3;
  

  

  /*
  for i in range(N):
      beta += random.random()*2*w_max
      while beta > w[index]:
          beta -= w[index]
          index = (index + 1) % N
      p3.append(p[index])
  p = p3
  */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}