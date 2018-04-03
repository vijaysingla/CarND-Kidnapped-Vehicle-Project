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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	double std_x,std_y,std_theta;
	// Setting the standard definiton of x, y, theta
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	// defining no. of particles
	num_particles = 20;
	// creating a random generator
	std::default_random_engine gen;

	// Creating normal distributions for x, y, theta
	std::normal_distribution<double> dist_x(x,std_x);
	std::normal_distribution<double> dist_y(y,std_y);
	std::normal_distribution<double> dist_theta(theta,std_theta);

	double weight = 1.0;

	//Initializing particles and weights
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
	//std:: cout<< "Initialization Done"<<endl;
	is_initialized = 1;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Setting the standard definiton of x, y, theta
	double std_x,std_y,std_theta;
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];
	std::default_random_engine gen;

	//Creating normal distributions for x, y, theta
	std::normal_distribution<double> dist_x(0,std_x);
	std::normal_distribution<double> dist_y(0,std_y);
	std::normal_distribution<double> dist_theta(0,std_theta);

	double v_by_yawrate,yawdt;
	if (fabs(yaw_rate) > 0.001)
	{
		v_by_yawrate = velocity/yaw_rate;

	}
	else
	{
		v_by_yawrate =0;

	}
	double vdt = velocity*delta_t;
	yawdt = yaw_rate*delta_t;
	for (int i =0;i< num_particles;i++)
	{
		double c_theta,s_theta,c_theta_dt,s_theta_dt;
		c_theta = cos(particles[i].theta);
		s_theta = sin(particles[i].theta);
		c_theta_dt = cos(particles[i].theta + yawdt);
		s_theta_dt = sin(particles[i].theta + yawdt);

		if (fabs(yaw_rate) > 0.001)
		{

			particles[i].x = particles[i].x + v_by_yawrate*(s_theta_dt-s_theta);
			particles[i].y = particles[i].y + v_by_yawrate*(c_theta-c_theta_dt);
			particles[i].theta = particles[i].theta + yawdt;
		}
		else
		{
			particles[i].x = particles[i].x + vdt*c_theta;
			particles[i].y = particles[i].y + vdt*s_theta;
		}




	particles[i].x = particles[i].x + dist_x(gen);
	particles[i].y = particles[i].y+dist_y(gen);
	particles[i].theta = particles[i].theta + dist_theta(gen);
	}
	//std:: cout<< "Prediction Done"<<endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

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

	double sig_x,sig_y;
	sig_x = std_landmark[0];
	sig_y = std_landmark[1];

	double sig2_x,sig2_y;
	sig2_x = 1/(2*sig_x*sig_x);
	sig2_y = 1/(2*sig_y*sig_y);

	double gauss_norm;
	gauss_norm = 1/(2*M_PI*sig_x*sig_y);
	int obs_size = observations.size();

	//Homogenous tranformation

	// Observation coordinates
	for (int i =0; i < num_particles; i++)
	{
		particles[i].associations.clear();
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		particles[i].weight =1;
		weights[i] = 1;

		double c_theta,s_theta;
		c_theta = cos(particles[i].theta);
		s_theta = sin(particles[i].theta);


		vector<LandmarkObs> obs_trans ;

		for (int j =0 ;j< obs_size;j++)
		{
			// Homogeneous transformation
			LandmarkObs trans;
			trans.x = particles[i].x + (c_theta*observations[j].x) -(s_theta*observations[j].y);
			trans.y = particles[i].y + (c_theta*observations[j].y) +(s_theta*observations[j].x);
			obs_trans.push_back(trans);

		}
//		std::cout<<"transformation done"<<endl;

		for (int j =0 ;j< obs_size;j++)
		{
			double prev_dst,new_dst;

			prev_dst = sensor_range;
			int min_dist_index;
			min_dist_index =0;
//			cout<<"list_size ="<<map_landmarks.landmark_list.size()<<endl;
			for (int k =0;k<map_landmarks.landmark_list.size();k++)

			{

				new_dst = dist(obs_trans[j].x,obs_trans[j].y,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
				if (new_dst < prev_dst)
				{
					prev_dst = new_dst;
					min_dist_index = k; //map_landmarks.landmark_list[k].id_i;

				}

			}
//			std::cout<< "i =" << i <<endl;
			//std::cout<< "k =" << min_dist_index <<endl;
			particles[i].associations.push_back(map_landmarks.landmark_list[min_dist_index].id_i);
			particles[i].sense_x.push_back(obs_trans[j].x);
			particles[i].sense_y.push_back(obs_trans[j].y);

			double mu_x,mu_y,x,y;

			mu_x= map_landmarks.landmark_list[min_dist_index].x_f;
			mu_y= map_landmarks.landmark_list[min_dist_index].y_f;
			x = obs_trans[j].x;
			y = obs_trans[j].y;

//			std::cout<< "x="<<x<<"  "<<"y="<<y<<"  "<<"mu_x="<<mu_x<<"  "<<"mu_y="<<mu_y<<"  "<<"sigx="<<sig_x<<" "<<"sigy="<<sig_y <<endl;
			double multi_variate_out;

			double exponent = ((x-mu_x)*(x-mu_x)*sig2_x) + ((y-mu_y)*(y-mu_y)*sig2_y);
			multi_variate_out = gauss_norm*exp(-exponent);
//			std::cout<<"variate ="<<multi_variate_out<<endl;
 			particles[i].weight = particles[i].weight*multi_variate_out;
//			std::cout << "weight= "<<particles[i].weight <<endl;
			weights[i] = particles[i].weight;
		}

	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::default_random_engine gen;
	std:: discrete_distribution<> resample(weights.begin(),weights.end());

	vector<Particle> particles_resample;
	for (int i =0;i<num_particles;i++)
	{
		particles_resample.push_back(particles[resample(gen)]);

	}
	particles = particles_resample;

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
