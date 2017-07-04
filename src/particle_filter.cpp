  
//  particle_filter.cpp
//  Initial version (about 50%) by Tiffany Huang -  2016-12-16
//  Functions added (about 50%) by Matthew Younkins - 2017-07-01
//  All poorly thought out algorithms by Matthew Younkins
//	...but it was Tiffany that used the std namespace

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

using namespace std;					//How could you, Tiffany?

std::default_random_engine gen;
std::normal_distribution<double> distn(0, 1); 

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// 1) Sets number of particles,
	// 2) initializes all particles to first position (based on GPS x,y,theta, uncertainty)
	// 3) Sets all initial weights to 1. 
	// 4) Adds random Gaussian noise to each particle.
	// More info is in particle_filter.h !  
	
	
	// How many particles should we use?
	// I did my development with 100, per suggestion, but lower numbers are possible.
	// Some Results:
	//
	//   N  X-err 	Y-err 	Yaw-err
	//10000 							***Error:  "You ran out of time"
	// 1000 0.108	0.098	0.004
	//	100 0.114	0.107	0.004
	// 	 20	0.136	0.128	0.005
	//	 10 0.174	0.151	0.006
	//	 10 0.181   0.149   0.006		<- Note change due to stochastic process
	//    9 							The whole kit-and-kaboodle fails... no blue dot movement
	//	  8								***Y-error is larger than maximum
	//	  4 							***Y-error is larger than maximum
	num_particles 	= 	16;				//Set the number of particles
  
  
	std::default_random_engine gen;
	std::normal_distribution<double> N_x(x, std[0]);   		//Function to add random noise in x
	std::normal_distribution<double> N_y(y, std[1]);		//Function to add random noise in y
	std::normal_distribution<double> N_theta(theta, std[2]);//Function to add random noise in theta
	
	Particle particle;
	
	for (int i = 0; i < num_particles; i++)					// Similar to method in https://youtu.be/-3HI3Iw3Z9g
	{			
		particle.id = i;									// Just an index
		particle.x = N_x(gen);								// Adds random x noise to particle
		particle.y = N_y(gen);								// Adds random y noise to particle
		particle.theta = N_theta(gen);						// Adds random theta noise
		particle.weight = 1.0;								//initialize weight to 1
		particles.push_back(particle);						
		weights.push_back(1);
	}
	is_initialized = true;
}			


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
// This function defines the measurements for each particle and adds noise. 
// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
// http://www.cplusplus.com/reference/random/default_random_engine/
// For reference, please see Term2, Lesson14, Section 8			
	
	std::default_random_engine gen;								 
	double 	new_x;
	double 	new_y;
	double 	new_theta;															  

	for (int i = 0; i < particles.size(); i++)
	{							
		if (fabs(yaw_rate) > 0.0001)
		{
			double velocityDivYawRate = velocity/yaw_rate;
			new_theta = particles[i].theta + yaw_rate * delta_t;
			//note that the functions below don't need to know delta_t because it's implicit in new_theta
			new_x = particles[i].x + velocityDivYawRate*(sin(new_theta) - sin(particles[i].theta));
			new_y = particles[i].y + velocityDivYawRate*(cos(particles[i].theta) - cos(new_theta));	   
		}
		else
		{
			new_theta = particles[i].theta;	
			new_x = particles[i].x + velocity * cos(particles[i].theta) * delta_t;
			new_y = particles[i].y + velocity * sin(particles[i].theta) * delta_t;		
		}

 		std::normal_distribution<double> N_x(new_x, std_pos[0]);	//need to add jitter to create new points
		std::normal_distribution<double> N_y(new_y, std_pos[1]);
		std::normal_distribution<double> N_theta(new_theta, std_pos[2]);	
		
		particles[i].x 		= N_x(gen);
		particles[i].y 		= N_y(gen);
		particles[i].theta 	= N_theta(gen);
	}
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations)
{
// This function finds the predicted measurement that is closest to each observed measurement.
// Then it assigns the observed measurement to that landmark.
// For reference, use Lesson 14.13
	for (int o = 0; o < observations.size(); o++)	//o = index to iterate through observations
	{
		for (int p = 0; p < predicted.size(); p++)	//p = index to iterate through predictions
		{
			double distanceToPointNew = dist(observations[o].x, observations[o].y, predicted[p].x, predicted[p].y);
			double distanceToPoint;
			if(p == 0)
			{
				distanceToPoint = distanceToPointNew;
				observations[o].id = p;
			}
			else if(distanceToPointNew < distanceToPoint)
			{
				distanceToPoint = distanceToPointNew;
                observations[o].id = p; // Index of matching landmark in the predicted_map
			}
		}
	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],std::vector<LandmarkObs> observations, Map map_landmarks) {
// This function updates the weighting applied to each particle with a mult-variate Gaussian distribution.
// https://en.wikipedia.org/wiki/Multivariate_normal_distribution
// Observations are initially given in the VEHICLE'S coordinate system.  
// Particles are located in the MAP's coordinate system.
// A transform (both rotation and translation) is used to convert between the two systems.
// This is very similar to what was accomplished in Lesson14.
// Extra info: https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
// pow(extra,2):   Equation 3.33 in http://planning.cs.uiuc.edu/node99.html

	for (int p = 0; p < particles.size(); p++)  						//iterate through every particle
	{
		
		std::vector<LandmarkObs> whereWeObservedThings;					//Where we found things
		for (int o = 0; o < observations.size(); o++)					//Transform each observation to global coordinates
		{
			LandmarkObs transformedCoordinatesOfObservation;
			//perform the space transformation from vehicle to map
			transformedCoordinatesOfObservation.x = particles[p].x+(observations[o].x*cos(particles[p].theta)-observations[o].y*sin(particles[p].theta));
			transformedCoordinatesOfObservation.y = particles[p].y+(observations[o].x*sin(particles[p].theta)+observations[o].y*cos(particles[p].theta));
			whereWeObservedThings.push_back(transformedCoordinatesOfObservation);
		}

		std::vector<LandmarkObs> whereThingsOughtToBe;
		
		for (int f = 0; f < map_landmarks.landmark_list.size(); f++)	//iterate through each landmark
		{
			if (dist(particles[p].x, particles[p].y, map_landmarks.landmark_list[f].x_f, map_landmarks.landmark_list[f].y_f) <= sensor_range)
			{
				LandmarkObs thisPoint;
				thisPoint.x = map_landmarks.landmark_list[f].x_f;		//just need to translate from landmark list to LandmarkObs parameters
				thisPoint.y = map_landmarks.landmark_list[f].y_f;
				thisPoint.id = map_landmarks.landmark_list[f].id_i;
				whereThingsOughtToBe.push_back(thisPoint);
			}
		}
		dataAssociation(whereThingsOughtToBe, whereWeObservedThings);

		particles[p].weight = 1.0;
		
		double preCalcStdLandmarkX = 1.0/std_landmark[0]*std_landmark[0];
		double preCalcStdLandmarkY = 1.0/std_landmark[1]*std_landmark[1];
		double preCalcMultiplier   = 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]);
		
		for (int o = 0; o < whereWeObservedThings.size(); o++)			//iterate through where we've made observations
		{
			double deltaX 	= whereWeObservedThings[o].x - whereThingsOughtToBe[ whereWeObservedThings[o].id ].x;
			double deltaY 	= whereWeObservedThings[o].y - whereThingsOughtToBe[ whereWeObservedThings[o].id ].y;
			double tempWeight = preCalcMultiplier * exp(-0.5 * (deltaX*deltaX*preCalcStdLandmarkX + deltaY*deltaY*preCalcStdLandmarkY));
			particles[p].weight *= tempWeight;
		}
		weights[p] 	= particles[p].weight;
	}
}


void ParticleFilter::resample()
{
	// Resamples particles (choose with replacement), with probability proportional to their weight. 
	// Some info here: http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std::default_random_engine gen;												
	std::discrete_distribution<int> distribution(weights.begin(), weights.end());
	std::vector<Particle> resample_particles;
	
	for (int i = 0; i < num_particles; i++) resample_particles.push_back(particles[distribution(gen)]);
	particles  = resample_particles;
}



//no change to the code below -- Matthew Younkins

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
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
	std::vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	std::vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	std::vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
