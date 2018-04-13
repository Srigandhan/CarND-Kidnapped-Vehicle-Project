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
	

	num_particles = 200;
	default_random_engine gen;
	
	// These lines creates a normal (Gaussian) distribution for x,y and theta.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	//Initialize each particle based on the gps position and using the guassian distribution with Standard Deviation given.
	// Also initialize the Weights vector.
	
	for (int i=0;i<num_particles;i++) {
		Particle p;
		p.id = i+1;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(1.0);
	}	
	
	is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//cout<<"Prediction Started";
	default_random_engine gen;
	
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	//Predict the postion using the previous time stamp values.
	
    for(Particle& p: particles){
        if( fabs(yaw_rate) < 0.0001){
            //if velocity is too close to zero
            p.x += velocity * delta_t * cos(p.theta);
            p.y += velocity * delta_t * sin(p.theta);
        } else{
            p.x += velocity / yaw_rate * ( sin( p.theta + yaw_rate*delta_t ) - sin(p.theta) );
            p.y += velocity / yaw_rate * ( cos( p.theta ) - cos( p.theta + yaw_rate*delta_t ) );
            p.theta += yaw_rate * delta_t;
        }

		//Adding Random Gaussian Noise
		p.x += dist_x(gen);
		p.y += dist_y(gen);
		p.theta += dist_theta(gen);		
	}	


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	
	
	for(int i=0;i<observations.size();i++)
	{
		double min = std::numeric_limits<double>::max();
		int id;
		for(int j=0;j<predicted.size();j++)
		{
			double distance = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
			
			if(distance<min)
			{
				min=distance;
				//Assign the Index to the Obs ID to fetch the Landmarks associated easily afterwards when needed.
				observations[i].id=j;
			}
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
	
	for (int i=0;i<num_particles;i++) {
		
		double px = particles[i].x;
		double py = particles[i].y;
		double ptheta = particles[i].theta;
	
		vector<LandmarkObs> obsmap;
		
		//Convert Observations from Vehicle Coordinates to Map Cooradinates.
		for(int j=0; j<observations.size(); j++)
		{
			LandmarkObs l;
			l.x = px + (cos(ptheta)*observations[j].x) - (sin(ptheta)*observations[j].y);
			l.y = py + (sin(ptheta)*observations[j].x) + (cos(ptheta)*observations[j].y);
			obsmap.push_back(l);
		}

		
		vector<LandmarkObs> LMs;
		
		
		//Find and store MAP Landmarks that are inside the sensor range
		
		for(int j=0;j<map_landmarks.landmark_list.size();j++)
		{
			double lx = map_landmarks.landmark_list[j].x_f;
			double ly = map_landmarks.landmark_list[j].y_f;
			double lid = map_landmarks.landmark_list[j].id_i;

			if((fabs(px-lx)<=sensor_range) && (fabs(py-ly)<=sensor_range))
				{
				LandmarkObs s;
				s.id = lid;
				s.x = lx;
				s.y = ly;
				LMs.push_back(s);
				}

		}

		dataAssociation(LMs,obsmap);
		
		particles[i].weight=1.0;
	
		//Calculate the weights of each observations and multiply the probabilities.
		for(int j=0;j<obsmap.size();j++)
		{
			
			double constant = 1 / (2*M_PI*std_landmark[0]*std_landmark[1]);
			double obsx = obsmap[j].x;
			double obsy = obsmap[j].y;

			double mux = LMs[obsmap[j].id].x;
			double muy = LMs[obsmap[j].id].y;

			
			double e =(pow((obsx-mux),2)/(2*pow(std_landmark[0],2))) + (pow((obsy-muy),2)/(2*pow(std_landmark[1],2)));


			
			double weight = constant*exp(-e);
			
			particles[i].weight*=weight;

		}
		
		//Store the weight of the current particle in the weights vector to use it in resampling.
		weights[i] = particles[i].weight;
	
	}

	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//Using discrete_distribution to resample.
	
	default_random_engine gen;
	vector<Particle> resampled;
	discrete_distribution<int> dist_w(weights.begin(), weights.end());
	for (auto &particle : particles){
		resampled.push_back(particles[dist_w(gen)]);
	}
	particles = resampled;
	
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
