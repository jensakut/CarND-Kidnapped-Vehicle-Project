/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 
 Exercise code included by Jens Kutschera
 July 2018
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <random>

#include "particle_filter.h"

#define PARTICLE_NUMBER 12;
#define EPS 0.001 


using namespace std;
	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = PARTICLE_NUMBER;
	default_random_engine gen;
	
	//create a normal distribution around gps position to draw particles from
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	// resize to number of particles. 
	particles.resize(num_particles);
	weights.resize(num_particles);
	double weight_init = 1/num_particles; 
	
	cout << num_particles << " samples" << endl;
	for (int i = 0; i < num_particles; ++i) {
		// Sample  and from these normal distrubtions like this: 
		// sample_x = dist_x(gen);
		// where "gen" is the random engine initialized earlier.

		particles[i].id=i;
		particles[i].x=dist_x(gen);
		particles[i].y=dist_y(gen);
		particles[i].theta=dist_theta(gen);
		particles[i].weight = weight_init;
		// Print your samples to the terminal.
		
		cout << "Sample " << particles[i].id << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
	}
	is_initialized = true;
	cout << "Initialized everything" << endl;
}

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
void ParticleFilter::prediction(double dt, double std_pos[], double v, double yd) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	cout << "Prediction: " << endl; 
	default_random_engine gen;
	normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);
	for (int i = 0; i < num_particles; ++i) {
		double theta=particles[i].theta; //better readability
		if (yd < EPS) {
			particles[i].x += v*dt*cos(theta)+dist_x(gen);
			particles[i].y += v*dt*sin(theta)+dist_y(gen);
			//theta unchanged when there isn't any yaw rate 
		}
		else {
			particles[i].x += v/yd*(sin(theta+yd*dt) - sin(theta))       + dist_x(gen);
			particles[i].y += v/yd*(cos(theta)       - cos(theta+yd*dt)) + dist_y(gen);
			particles[i].theta += yd*dt  +  dist_theta(gen);
		}
		// Print your samples to the terminal.
		cout << "Sample " << particles[i].id << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
	}
	
	
}
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	 
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < int(observations.size()); ++i) {
		double minDistance=numeric_limits<double>::max();
		
		for (int j = 0; j < int(predicted.size()); ++j) {
			// nearest neighbor implementation: 
			double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (distance < minDistance) {
				minDistance=distance;
				observations[i].id=predicted[j].id;
			}
		}
	}
	
}
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	cout << "update" << endl; 
	double sigma_x2 = std_landmark[0]*std_landmark[0];
	double sigma_y2 = std_landmark[1]*std_landmark[1];
	// only dependand on std_landmark
	double gaussNormalizer = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));

	for (int i = 0; i < num_particles; ++i) {
		//transform each observation into transformed_observations
		std::vector<LandmarkObs> transformed_observations;
		for (int j = 0; j < int(observations.size()); j++){
			double x = cos(particles[i].theta)*observations[j].x-sin(particles[i].theta)*observations[j].y+particles[i].x;
			double y = sin(particles[i].theta)*observations[j].x+cos(particles[i].theta)*observations[j].y+particles[i].y;
			transformed_observations.push_back({j,x,y});  
		}
		//transform each landmark in range into transformed_observations
		std::vector<LandmarkObs> predictions;

		for (int j = 0; j < int(map_landmarks.landmark_list.size()); ++j) {
			//variables for speed and readability
			float landmark_x = map_landmarks.landmark_list[j].x_f;
			float landmark_y = map_landmarks.landmark_list[j].y_f;
			int landmark_id  = map_landmarks.landmark_list[j].id_i;
			//only account for landmarks in range. Dist should be reasonably fast
			if (dist(landmark_x,landmark_y,particles[i].x,particles[i].y) <=sensor_range) {
				predictions.push_back({landmark_id, landmark_x, landmark_y});
			}
		}
		
		dataAssociation(predictions, transformed_observations);
		particles[i].weight=1.0;
		
		for (int j = 0; j<int(transformed_observations.size()); ++j){
			//the nearest neighbor landmark out of dataAssociation. 
			int associated_id = transformed_observations[j].id;
			// default the x and y values to something totally off in case the id cant be retrieved
			// this weight will be zero in that exception. 
			double predicted_x = numeric_limits<double>::max();
			double predicted_y = numeric_limits<double>::max();
			
			
			//grab x and y of this landmark. 
			for (int k = 0; k<int(predictions.size()); k++) {
				if (predictions[k].id == associated_id)  {
					predicted_x = predictions[k].x;
					predicted_y = predictions[k].y;
				}
			}
			if (predicted_x = numeric_limits<double>::max()) 
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl; 
			//the weight for this 
			//using the precomputed parts from above
			//exp = (x-mu_x)^2/sigma_x^2+(y-mu_y)^2/sigma_y^2 excluding -1/2 
			double exponent= (predicted_x - transformed_observations[j].y)*(predicted_x - transformed_observations[j].x)/sigma_x2+
							 (predicted_y - transformed_observations[j].y)*(predicted_y - transformed_observations[j].y)/sigma_y2;
			//weight = factor * e ^(-1*exponent/2) 
			//cout << "Norm " << gaussNormalizer << " exp " << exponent << endl; 
			particles[i].weight *= gaussNormalizer * exp(-exponent/2);
		}
		cout << "ID " << particles[i].id << " weight " << particles[i].weight << endl; 
	}			


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	/*
  std::default_random_engine generator;
  std::discrete_distribution<int> distribution {2,2,1,1,2,2,1,1,2,2};

  int p[10]={};

  for (int i=0; i<nrolls; ++i) {
    int number = distribution(generator);
    ++p[number];
  }

  std::cout << "a discrete_distribution:" << std::endl;
  for (int i=0; i<10; ++i)
    std::cout << i << ": " << std::string(p[i]*nstars/nrolls,'*') << std::endl;
*/
	// new preallocated vector with weights
	std::vector<double> weights(num_particles); 
	for (int i = 0; i<num_particles; ++i) {
		weights[i] = particles[i].weight;
	}
	//create random number generator
	default_random_engine generator;
	discrete_distribution<> distribution(weights.begin(), weights.end());
	//new set of particles
	std::vector<Particle> new_particles(num_particles);
	for (int i = 0; i<num_particles; ++i) {
		int number = distribution(generator);
		//fill new particles with drawn old particles
		new_particles[i]=particles[number];
		//visualize what has been drawn
		cout << "Sample " << new_particles[i].id << " " << new_particles[i].x << " " << new_particles[i].y << " " << new_particles[i].theta << endl;
		//set id fow the following cycle. 
		new_particles[i].id=i;
	}
	//set new particles
	particles=new_particles; 
	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	
	return particle; 
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
