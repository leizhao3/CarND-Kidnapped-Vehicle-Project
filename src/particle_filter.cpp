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
#include <iomanip>
#include <iterator>
#include <numeric>
#include <random> // Need this for sampling from distributions
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::setw;

void display_part_obsmap(Particle part, vector<LandmarkObs> &observations_map);

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * init Initializes particle filter by initializing particles to Gaussian
   *   distribution around first position and all the weights to 1.
   * @param x Initial x position [m] (simulated estimate from GPS)
   * @param y Initial y position [m]
   * @param theta Initial orientation [rad]
   * @param std[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   */

  num_particles = 100;  // Set the number of particles
  
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

  /*//Checking the initial parameters
  int width = 15;
  std::cout <<"During Initializing......................" << std::endl;
  std::cout << setw(width)<< "x" 
            << setw(width)<< "y" 
            << setw(width)<< "theta"
            << setw(width)<< "weight"
            << std::endl;*/

  for(int i=0; i<num_particles; i++){
    // Initialize all parameter
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0f;


    particles.push_back(particle);

    /*//Checking the initial parameters
    std::cout << setw(width)<< particles[i].x
              << setw(width)<< particles[i].y
              << setw(width)<< particles[i].theta
              << setw(width)<< particles[i].weight
              << std::endl;*/

  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * dataAssociation Finds which observations correspond to which landmarks 
   *   (likely by using a nearest-neighbors data association).
   * @param predicted Vector of predicted landmark observations
   * @param map_landmarks Map class containing map landmarks
   */



  /*//Check the prediction function
  int width = 15;
  std::cout <<"During Predicting......................" << std::endl;
  std::cout << setw(width)<< "x" 
            << setw(width)<< "y" 
            << setw(width)<< "theta" 
            << setw(width)<< "v" 
            << setw(width)<< "theta_dot" 
            << setw(width)<< "x_pred" 
            << setw(width)<< "y_pred" 
            << setw(width)<< "theta_pred" 
            << std::endl;*/
  
  std::default_random_engine gen;
  const double std_x = std_pos[0];
  const double std_y = std_pos[1];
  const double std_theta = std_pos[2];


  for(int i=0; i<num_particles; i++){
    /*//Check the prediction function
    std::cout << setw(width)<< particles[i].x
              << setw(width)<< particles[i].y
              << setw(width)<< particles[i].theta
              << setw(width)<< velocity
              << setw(width)<< yaw_rate;*/
    double theta_f;
    double x_f;
    double y_f;

    if(yaw_rate != 0) {
      theta_f = particles[i].theta + yaw_rate*delta_t;
      x_f = particles[i].x + velocity/yaw_rate*(sin(theta_f)-sin(particles[i].theta));
      y_f = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(theta_f));
    } else {
      theta_f = particles[i].theta + yaw_rate*delta_t;
      double distance = velocity*delta_t;
      x_f = particles[i].x + cos(distance);
      y_f = particles[i].y + sin(distance);
    }

    normal_distribution<double> dist_x(x_f, std_x);
    normal_distribution<double> dist_y(y_f, std_y);
    normal_distribution<double> dist_theta(theta_f, std_theta);

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

    /*//Check the prediction function
    std::cout << setw(width)<< particles[i].x
              << setw(width)<< particles[i].y
              << setw(width)<< particles[i].theta
              << std::endl;*/

  }



}

void ParticleFilter::dataAssociation(LandmarkObs &obs_map, 
                                     const Map &map_landmarks) {
  /**
   * dataAssociation Finds which observations correspond to which landmarks 
   *   (likely by using a nearest-neighbors data association).
   * @param obs_map predicted landmark observations in MAP coordiate
   * @param map_landmarks Map class containing map landmarks
   */

  
  int dist_min_ind = 0;
  double dist_j;
  double dist_min = dist(obs_map.x, obs_map.y, 
                          map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);
  

  for(int j=0; j<map_landmarks.landmark_list.size(); j++) {
    dist_j = dist(obs_map.x, obs_map.y, 
                  map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
    //std::cout << "dist_j is " << dist_j << std::endl;
    if(dist_j < dist_min) {
      dist_min = dist_j;
      dist_min_ind = j;
      /* Check the functionality of the data association
      std::cout << "Get smaller distance!!! " << std::endl;
      std::cout << "dist_min = " << dist_min << std::endl;
      std::cout << "dist_min_ind = " << dist_min_ind << std::endl;
      std::cout << "map_landmarks.landmark_list[dist_min_ind].id_i = " 
                << map_landmarks.landmark_list[dist_min_ind].id_i << std::endl;*/


    }
  }
  obs_map.id = map_landmarks.landmark_list[dist_min_ind].id_i;
  
  

}

vector<LandmarkObs> ParticleFilter::VehicleCoor2MapCoor(Particle part,
                                                        const vector<LandmarkObs> &observations) {
  /**
   * VehicleCoor2MapCoor convert the coordinate from VEHICLE to MAP
   * @param part single particle, in VEHICLE coordinate
   * @param observations observations, in VEHICLE coordinate
   */

    vector<LandmarkObs> observations_map = {};
    LandmarkObs observations_map_temp;

    for(int j=0; j<observations.size(); j++) {

      observations_map_temp.x = part.x + (cos(part.theta)*observations[j].x) - (sin(part.theta)*observations[j].y);
      observations_map_temp.y = part.y + (sin(part.theta)*observations[j].x) + (cos(part.theta)*observations[j].y);

      observations_map.push_back(observations_map_temp);
    }

  

  return observations_map;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * updateWeights Updates the weights for each particle based on the likelihood
   *   of the observed measurements. 
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2
   *   [Landmark measurement uncertainty [x [m], y [m]]]
   * @param observations Vector of landmark observations, in VEHICLE coordinate
   * @param map Map class containing map landmarks
   */

  //Trouble Shooting
  if(observations.size() <= 0) {
    std::cout << "no observations ==> no landmark is detected." << std::endl;
    return;
  }

  //Define the observations standard deviation
  const double sigma_x = std_landmark[0];
  const double sigma_y = std_landmark[1];

  /*//Map output dislplay
  std::cout << "Map is listed below"  << std::endl;
  std::cout << setw(15)<< "id_i"
            << setw(15)<< "x_f"
            << setw(15)<< "y_f"
            << std::endl;
  
  for(int j=0; j<map_landmarks.landmark_list.size(); j++) {
    std::cout << setw(15)<< map_landmarks.landmark_list[j].id_i
              << setw(15)<< map_landmarks.landmark_list[j].x_f
              << setw(15)<< map_landmarks.landmark_list[j].y_f
              << std::endl;
  }*/
  


  for(int i=0; i<num_particles; i++) {
    
    vector<LandmarkObs> observations_map;

    vector<int> associations = {}; 
    vector<double> sense_x = {};
    vector<double> sense_y = {};

    int width = 15;
    /*//Single particle & corresponding landmark display
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "i = " << i << std::endl;
    std::cout << setw(width)<< "j" 
              << setw(width)<< "Obs.x" 
              << setw(width)<< "Obs.y" 
              << setw(width)<< "Obs_map.x" 
              << setw(width)<< "Obs_map.y" 
              << setw(width)<< "Landmark ID" 
              << setw(width)<< "mu_x" 
              << setw(width)<< "mu_y" 
              << setw(width)<< "Posterior" 
              << std::endl;*/

    // Convert observations from VEHICLE to MAP coordiate. 
    observations_map = VehicleCoor2MapCoor(particles[i], observations); 

    double weight_temp = 1;
    double posterior = 1;
    const double normalizer = 1.0/(2.0*M_PI*sigma_x*sigma_y);

    //Look for landmark in the range
    Map map_landmarks_inrange = {};
    for(int k=0; k<map_landmarks.landmark_list.size(); k++) {
      double obs_distance = dist(particles[i].x, particles[i].y, 
                                  map_landmarks.landmark_list[k].x_f, 
                                  map_landmarks.landmark_list[k].y_f);
      if(obs_distance < sensor_range) {
        map_landmarks_inrange.landmark_list.push_back(map_landmarks.landmark_list[k]);
      }
    }

    for(int j=0; j<observations.size(); j++){

      dataAssociation(observations_map[j], map_landmarks_inrange);

      double x = observations_map[j].x;
      double y = observations_map[j].y;
      int landmark_index = -1;

      for(unsigned int m=0; m<map_landmarks_inrange.landmark_list.size(); m++) {
        if(observations_map[j].id == map_landmarks_inrange.landmark_list[m].id_i) {
          landmark_index = m;
        } 
      }

      double mu_x = map_landmarks_inrange.landmark_list[landmark_index].x_f;
      double mu_y = map_landmarks_inrange.landmark_list[landmark_index].y_f;

      //Apply Multivariate-Gaussian probability density
      double exponent = exp(-(pow((x-mu_x),2)/(2*pow(sigma_x,2))+
                              pow((y-mu_y),2)/(2*pow(sigma_x,2))));
    
      posterior = normalizer * exponent;
      weight_temp *= posterior;
      
      

      /*//Single particle & corresponding landmark display
      std::cout << setw(width)<< j
                << setw(width)<< observations[j].x << ","
                << setw(width)<< observations[j].y << ","
                << setw(width)<< observations_map[j].x << ","
                << setw(width)<< observations_map[j].y << ","
                << setw(width)<< map_landmarks.landmark_list[landmark_index].id_i
                << setw(width)<< mu_x
                << setw(width)<< mu_y;
      if(obs_range < sensor_range) {
        std::cout << setw(width)<< posterior;
      } else {
        std::cout << setw(width)<< "OoRange";
      }
      std::cout << std::endl;*/


      associations.push_back(observations_map[j].id);
      sense_x.push_back(observations_map[j].x);//x, in MAP coordiate
      sense_y.push_back(observations_map[j].y);//y, in MAP coordiate
    }
    particles[i].weight = weight_temp;

    SetAssociations(particles[i],associations, sense_x, sense_y);

    /*//Display single particle & observation in MAP coordiate
    std::cout << "particles[i].weight AFTER =" << particles[i].weight << std::endl;
    display_part_obsmap(particles[i], observations_map);
    */

  }


}

void ParticleFilter::resample() {
  /**
   * resample Resamples from the updated set of particles to form
   *   the new set of particles.
   */

  std::random_device rd;
  std::mt19937 gen(rd());

  vector<double> weights_Array = {};
  for(int i=0; i<num_particles; i++) {
    weights_Array.push_back(particles[i].weight);
  }

  std::discrete_distribution<> distr(weights_Array.begin(),weights_Array.end());

  vector<Particle> particles_3;

  
  //std::cout << "druing resampling............... Index Below" << std::endl;
  for(int i=0; i<num_particles; i++) {
    int index = distr(gen);
    //std::cout << index << " ";
    particles_3.push_back(particles[index]);
  }
  std::cout << std::endl;

  particles = particles_3;

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


void display_part_obsmap(Particle part, vector<LandmarkObs> &observations_map) {
  std::cout << "particles[i] = [" << part.x << "," << part.y << "," << part.theta << "]" << std::endl;
            
  std::cout << "observations_map.x =" << std::endl;
  std::cout << "[";
  for(int m=0; m<observations_map.size(); m++) {
    std::cout << observations_map[m].x << ",";
  }
  std::cout << "]" << std::endl;

  std::cout << "observations_map.y =" << std::endl;
  std::cout << "[";
  for(int m=0; m<observations_map.size(); m++) {
    std::cout << observations_map[m].y << ",";
  }
  std::cout << "]" << std::endl;

  std::cout << "landmark ID =" << std::endl;
  std::cout << "[";
  for(int m=0; m<observations_map.size(); m++) {
    std::cout << observations_map[m].id << ","; 
  }
  std::cout << "]" << std::endl;
}