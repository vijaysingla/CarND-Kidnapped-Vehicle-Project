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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// This function initializes all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Random Gaussian noise is added to each particle.

	// Setting the standard deviation of x, y, theta
	double std_x,std_y,std_theta;

	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	// defining no. of particles
	num_particles = 50;

	// creating a random generator
	std::default_random_engine gen;

	// Creating normal distributions for x, y, theta
	std::normal_distribution<double> dist_x(x,std_x);
	std::normal_distribution<double> dist_y(y,std_y);
	std::normal_distribution<double> dist_theta(theta,std_theta);

	double weight = 1.0;

	//Initializing particles and weights using GPS data and adding random noise
	Particle particle;
	for (int i =0;i< num_particles;i++)
	{
	particle.id = i;
	particle.x = dist_x(gen);
	particle.y = dist_y(gen);
	particle.theta = dist_theta(gen);
	particle.weight = 1;
	particles.push_back(particle);
	weights.push_back(weight);
	}

	// Setting the initialization boolean to 1 after initializing
	is_initialized = 1;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// This function is being used to predict the next state for particles after delta_t using control parameters
	//such as velocity and yaw rate

	// Setting the standard deviation of x, y, theta
	double std_x,std_y,std_theta;
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];

	//  Creating a random generator
	std::default_random_engine gen;

	//Creating normal distributions for x, y, theta with zero mean
	std::normal_distribution<double> dist_x(0,std_x);
	std::normal_distribution<double> dist_y(0,std_y);
	std::normal_distribution<double> dist_theta(0,std_theta);

	// defining the variables to be used in prediction calculations
	double v_by_yawrate,yawdt,vdt;

	if (fabs(yaw_rate) > 0.001)
	{
		v_by_yawrate = velocity/yaw_rate;

	}
	else
	{
		v_by_yawrate =0;

	}
	vdt = velocity*delta_t;
	yawdt = yaw_rate*delta_t;

	// prediction step is being done for each particle
	for (int i =0;i< num_particles;i++)
	{
		// prediction step when yaw rate is not equal to zero
		if (fabs(yaw_rate) > 0.001)
		{

			particles[i].x = particles[i].x + v_by_yawrate*(sin(particles[i].theta + yawdt)-sin(particles[i].theta));
			particles[i].y = particles[i].y + v_by_yawrate*(cos(particles[i].theta)-cos(particles[i].theta + yawdt));
			particles[i].theta = particles[i].theta + yawdt;
		}
		// prediction step when yaw rate is equal to zero
		else
		{
			particles[i].x = particles[i].x + vdt*cos(particles[i].theta);
			particles[i].y = particles[i].y + vdt*sin(particles[i].theta);
		}



	// Adding Gaussian noise to the particles
	particles[i].x = particles[i].x + dist_x(gen);
	particles[i].y = particles[i].y+dist_y(gen);
	particles[i].theta = particles[i].theta + dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// THis function can be defined to  find the predicted measurement that is closest
	//to each observed measurement and assign the observed measurement to this particular landmark.
	//   THis function is not defined as this is not being used in the code.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	// This function is used to update the weights of each particle using a mult-variate Gaussian distribution.
	// Information about this distribution can be found here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

	// NOTE: The observations are given in the VEHICLE'S coordinate system. THe particles are located
	//   according to the MAP'S coordinate system. Transform between the two systems is carried out in this function
	//   This transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// Defining the standard deviation for landmarks measurements
	double sig_x,sig_y;
	sig_x = std_landmark[0];
	sig_y = std_landmark[1];

	// Defining some variable independent of particle to be used for multi-variate gaussian distribution
	double sig2_x,sig2_y,gauss_norm;
	sig2_x = 1/(2*sig_x*sig_x);
	sig2_y = 1/(2*sig_y*sig_y);
	gauss_norm = 1/(2*M_PI*sig_x*sig_y);

	// Defining the size of observations vector
	int obs_size = observations.size();



	// Weight Update Step
	for (int i =0; i < num_particles; i++)
	{
		// Clearing out particles associations, sense_x and sense_y from previous step

		particles[i].associations.clear();
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		// Resetting all the weights to 1
		particles[i].weight =1;
		weights[i] = 1;

		// Defining some variables to be used for homogenous transformation
		double c_theta,s_theta;
		c_theta = cos(particles[i].theta);
		s_theta = sin(particles[i].theta);

		// Defining a vector for transformed observations
		vector<LandmarkObs> obs_trans ;

		// Homogeneous transformation is perfomed for every observation corresponding to each particle
		for (int j =0 ;j< obs_size;j++)
		{

			LandmarkObs trans;
			trans.x = particles[i].x + (c_theta*observations[j].x) -(s_theta*observations[j].y);
			trans.y = particles[i].y + (c_theta*observations[j].y) +(s_theta*observations[j].x);
			obs_trans.push_back(trans);

		}

		//Finding the closest landmark association to every transformed observation and updating the weight
		//for every particle
		for (int j =0 ;j< obs_size;j++)
		{
			double prev_dst,new_dst;
			prev_dst = sensor_range;

			// Closest landmark observation for each transformed observation is min_dist_index
			int min_dist_index;
			min_dist_index =0;
			// Finding the closest landmark association to each transformed observation
			for (int k =0;k<map_landmarks.landmark_list.size();k++)

			{

				new_dst = dist(obs_trans[j].x,obs_trans[j].y,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
				if (new_dst < prev_dst)
				{
					prev_dst = new_dst;
					min_dist_index = k;

				}

			}

			// Setting the particles associations, sense_x,sense_y
			particles[i].associations.push_back(map_landmarks.landmark_list[min_dist_index].id_i);
			particles[i].sense_x.push_back(obs_trans[j].x);
			particles[i].sense_y.push_back(obs_trans[j].y);

			// Inputs for  multi variate gaussian distribution
			double mu_x,mu_y,x,y;
			// Landmark map locations
			mu_x= map_landmarks.landmark_list[min_dist_index].x_f;
			mu_y= map_landmarks.landmark_list[min_dist_index].y_f;
			// Transformed landmark observations
			x = obs_trans[j].x;
			y = obs_trans[j].y;

			// Multi Variate Gaussian distribution output
			double multi_variate_out;
			double exponent = ((x-mu_x)*(x-mu_x)*sig2_x) + ((y-mu_y)*(y-mu_y)*sig2_y);
			multi_variate_out = gauss_norm*exp(-exponent);

			// weights_update
 			particles[i].weight = particles[i].weight*multi_variate_out;
			weights[i] = particles[i].weight;
		}

	}


}

void ParticleFilter::resample() {
	// This function is used to  resample particles with replacement with probability proportional to their weight.
	// std::discrete_distribution is used to achieve that
	//  http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Creating a random generator
	std::default_random_engine gen;

	// Creating a discrete distribution for all the particle weights
	std:: discrete_distribution<> resample(weights.begin(),weights.end());

	// Creating a vector of resampled particles
	vector<Particle> particles_resample;

	// For N particles, drawing a particle depending on particle weights
	for (int i =0;i<num_particles;i++)
	{
		particles_resample.push_back(particles[resample(gen)]);

	}
	// Setting particles vector equal to vector of resampled particles
	particles = particles_resample;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	// This function can be defined and used to set particles associations,sense_x,sense_y
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
	// This functions is not being used as particle associations,sense_x,sense_y are defined in Updateweights funcitons

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
