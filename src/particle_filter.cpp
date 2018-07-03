/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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

#include "particle_filter.h"

#define PARTICLE_NUMBER 1000; //success requires 10 particles 
#define EPS 0.001 

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = PARTICLE_NUMBER;
	//static default_random_engine gen;
	
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

		//cout << "Sample " << particles[i].id << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
	}
	is_initialized = true;
	//cout << "Initialized everything" << endl;
}

void ParticleFilter::prediction(double dt, double std_pos[], double v, double yd) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//cout << "Prediction: " << endl; 
	//default_random_engine gen;
	normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);
	for (int i = 0; i < num_particles; ++i) {
		double theta=particles[i].theta; //better readability
		if (fabs(yd) < EPS) {
			particles[i].x += v*dt*cos(theta)+dist_x(gen);
			particles[i].y += v*dt*sin(theta)+dist_y(gen);
			//theta unchanged when there isn't any yaw rate 
			//cout << "0000000000000000000000000000000000000000000" << endl; 
		}
		else {
			double yddt=yd*dt;
			double v_yd=v/yd;
			particles[i].x += v_yd*(sin(theta+yddt) - sin(theta))       + dist_x(gen);
			particles[i].y += v_yd*(cos(theta)       - cos(theta+yddt)) + dist_y(gen);
			particles[i].theta += yddt  +  dist_theta(gen);
			//cout << "1111111111111111111111111111111111111111111111111111" << endl; 
		}
		// Print your samples to the terminal.
		//cout << "Sample " << particles[i].id << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
	}
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	//cout << observations.size()*predicted.size() << endl; 
	for (LandmarkObs& observation: observations) { 
		double minDist = numeric_limits<double>::max();
		for(LandmarkObs& prediction: predicted) {
			double distance = dist(prediction.x,prediction.y,observation.x,observation.y);
			if (distance < minDist) {
				minDist = distance; 
				observation.id=prediction.id;
			}				
			//cout << "x " << observation.x << " " << prediction.x << endl; 
			//cout << "y " << observation.x << " " << prediction.y << endl; 
		}
	}
}

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
	//cout << "updateWeights" << endl; 
	const double sigma_x2 = 2*std_landmark[0]*std_landmark[0];
	const double sigma_y2 = 2*std_landmark[1]*std_landmark[1];
	// only dependand on std_landmark
	const double gaussNormalizer = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
	
	
	// Cycle through particles [i]
	for (int i = 0; i<num_particles; i++) {
		//cout << "Particle Nr. " << i << endl; 
		//cycle through observations and transform them to map coordinates
		std::vector<LandmarkObs> transformed_observations;
		double px = particles[i].x;
		double py = particles[i].y;
		double pt = particles[i].theta;
		for (int j = 0; j < int(observations.size()); j++){
			double ox = observations[j].x;
			double oy = observations[j].y;
			double xMap = cos(pt)*ox-sin(pt)*oy + px;
			double yMap = sin(pt)*ox+cos(pt)*oy + py;
			transformed_observations.push_back({observations[j].id,xMap,yMap});  
		}
		
		//cycle through landmarks and exclude the ones out of range 
		// write them into a vector with the needed format. 
		std::vector<LandmarkObs> predictions;
		for (int k = 0; k < int(map_landmarks.landmark_list.size()); k++) {
			float lx = map_landmarks.landmark_list[k].x_f;
			float ly = map_landmarks.landmark_list[k].y_f;
			int   lid= map_landmarks.landmark_list[k].id_i;
			if (dist(px, py, lx, ly) <= sensor_range)
				predictions.push_back({lid, lx, ly});
		}
		
		dataAssociation(predictions, transformed_observations);
		/*for (LandmarkObs& observation: transformed_observations) { 
			cout << observation.id << " "; 
		}		
		cout << "after dataAssociation" << endl; */
		
		//set particle weight to 1
		particles[i].weight = 1.0; 
		
		for (LandmarkObs& observation: transformed_observations) { 
			//find the corresponding x and y value out of the prediction 
			double predx = numeric_limits<double>::max();
			double predy = numeric_limits<double>::max();
			for(LandmarkObs& prediction: predictions) {
				if (prediction.id == observation.id) {
					predx = prediction.x;
					predy = prediction.y;
					break;
				}
			}
			if (predx==numeric_limits<double>::max()) {
				cout << "noIDfound noIDfound noIDfound noIDfound noIDfound noIDfound noIDfound noIDfound noIDfound noIDfound" << endl; 
			}
			//the weight for this 
			//using the precomputed parts from above
			//exp = (x-mu_x)^2/sigma_x^2+(y-mu_y)^2/sigma_y^2 excluding -1/2 
			double exponent= (predx - observation.x)*(predx - observation.x)/sigma_x2+
							 (predy - observation.y)*(predy - observation.y)/sigma_y2;
			//weight = factor * e ^(-1*exponent/2) 
			//cout << "Norm " << gaussNormalizer << " exp " << exponent << endl; 
			particles[i].weight *= gaussNormalizer * exp(-exponent);		
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	
	// new preallocated vector with weights
	std::vector<double> weights(num_particles); 
	for (int i = 0; i<num_particles; ++i) {
		weights[i] = particles[i].weight;
	}
	discrete_distribution<> distribution(weights.begin(), weights.end());
	//new set of particles
	std::vector<Particle> new_particles(num_particles);
	for (int i = 0; i<num_particles; ++i) {
		int number = distribution(gen);
		//fill new particles with drawn old particles
		new_particles[i]=particles[number];
		//visualize what has been drawn
		//cout << "Sample " << new_particles[i].id << " " << new_particles[i].x << " " << new_particles[i].y << " " << new_particles[i].theta << endl;
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

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();
	
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