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

// declare a random engine to be used across multiple and various method calls
static default_random_engine gen;
//Debugging
//static int counter2 = 0;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 10;

	normal_distribution<double> N_x_init(0.0, std[0]);
	normal_distribution<double> N_y_init(0.0, std[1]);
	normal_distribution<double> N_theta_init(0.0, std[2]);
	Debugging_bool = false;

	for (int i = 0; i < num_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = x + N_x_init(gen);
		p.y = y + N_y_init(gen);
		p.theta = theta + N_theta_init(gen);
		p.weight = 1.0;
		particles.push_back(p);
	}

	is_initialized = true;
	cout << "particles vector size: " << particles.size() << endl;
	std::cout << particles[2].id << endl << particles[2].x << endl << particles[2].y
		<< particles[2].theta << endl << particles[2].weight << endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//describing a normal distribution with std deviation of the control measurements
	normal_distribution<double> x_control(0, std_pos[0]);
	normal_distribution<double> y_control(0, std_pos[1]);
	normal_distribution<double> theta_control(0, std_pos[2]);

	if (yaw_rate == 0)
	{
		for (int i = 0; i < num_particles; i++)
		{
			//these equations are to be used when the yaw_rate = 0
			particles[i].x += velocity*delta_t*cos(particles[i].theta);
			particles[i].y += velocity*delta_t*sin(particles[i].theta);

			//Adding noise to the positions and directions of the particles
			particles[i].x += x_control(gen);
			particles[i].y += y_control(gen);
			particles[i].theta += theta_control(gen);

		}
	}

	else
	{
		double factor = velocity / yaw_rate;

		for (int i = 0; i < num_particles; i++)
		{
			//updating the position of the particles according to the given control measurements
			double theta_f = particles[i].theta + yaw_rate*delta_t;
			particles[i].x += factor*(sin(theta_f) - sin(particles[i].theta));
			particles[i].y += factor*(cos(particles[i].theta) - cos(theta_f));
			particles[i].theta = theta_f;

			//Adding control noise to the position and diretion of the predicted particles
			particles[i].x += x_control(gen);
			particles[i].y += y_control(gen);
			particles[i].theta += theta_control(gen);
		}
	}


	//std::cout << "prediction works";


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//current_dist is given a random large number value
	double Large_Number = 10000;
	double prev_dist = Large_Number;
	double distance;

	//cout << "before dataA" << endl << "lm, obs size: " << endl;
	//cout << predicted.size() << " , " << observations.size() << endl;

	for (int i = 0; i < observations.size(); i++)
	{
		prev_dist = Large_Number;
		LandmarkObs o = observations[i];
		int map_id = -1;
		for (int j = 0; j < predicted.size(); j++)
		{
			LandmarkObs p = predicted[j];
			distance = dist(o.x, o.y, p.x, p.y);
			if (distance < prev_dist) {
				prev_dist = distance;
				map_id = predicted[j].id;
			}
		}
		observations[i].id = map_id;

		//Each observation is now associated with the nearest particle, with the observation
		//id given the value of the landmark id.
	}
	//std::cout << "data association works";

	//cout << "before dataA" << endl << "lm, obs size: " << endl;
	//cout << predicted.size() << " , " << observations.size() << endl;

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

	//calculate normalization term


	//Debugging
	int obs_len = observations.size();
	std::vector<LandmarkObs> observation_vector;
	//Debugging: cout << "obs len: " << obs_len;

	//created a new vector observation_vector to try and avoid any garbage values
	for (int i = 0; i < obs_len; i++)
	{
		observation_vector.push_back(observations[i]);
	}




	//Debugging cout << "observation_vector size: " << observation_vector.size() << endl;

	//Debugging
	/*
	for (int i = 0; i < observations.size(); i++)
	{
	cout<< i << " : (" << observations[i].x << "," << observations[i].y << ")" << endl;
	}
	*/




	//std::cout << "Start of UW";
	double gauss_norm = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));

	//clearing all the previously added weights
	weights.clear();


	for (int i = 0; i < num_particles; i++)
	{
		//Debugging: cout << "for the " << i << "particle" << endl;

		//edit
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		//Create a vector of all the landmarks within the sensor range
		vector<LandmarkObs> nearby_landmarks;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
		{
			double lm_x = map_landmarks.landmark_list[j].x_f;
			double lm_y = map_landmarks.landmark_list[j].y_f;
			int lm_id = map_landmarks.landmark_list[j].id_i;
			//edit
			if (dist(p_x, p_y, lm_x, lm_y)<sensor_range)
			{
				nearby_landmarks.push_back(LandmarkObs{ lm_id, lm_x, lm_y });

			}
		}

		//Create a vector of all the transformed observations for the particle
		vector<LandmarkObs> transformed_obs;
		transformed_obs.clear();
		for (int k = 0; k < observations.size(); k++)
		{

			double x_transformed = observations[k].x * cos(particles[i].theta)
				- observations[k].y * sin(particles[i].theta) + particles[i].x;

			double y_transformed = observations[k].x * sin(particles[i].theta)
				+ observations[k].y * cos(particles[i].theta) + particles[i].y;

			transformed_obs.push_back(LandmarkObs{ observations[k].id,
				x_transformed, y_transformed });
		}



		//Debugging
		/*
		cout << "	before dataA" << endl << "	nearby_lm, obs size: " << endl;
		cout << "	" << nearby_landmarks.size() << " , " << transformed_obs.size() << endl;
		cout << "	after dataA" << endl << "	lm, obs size: " << endl;
		cout << "	" << nearby_landmarks.size() << " , " << transformed_obs.size() << endl;
		*/

		//Debugging
		//counter2 += 1;

		dataAssociation(nearby_landmarks, transformed_obs);

		//Debugging
		/*
		cout << "nearby Landmarks size and vector along with id: " << nearby_landmarks.size() << endl;
		for (int i = 0; i < nearby_landmarks.size(); i++)
		{
		cout << i << " : (" << nearby_landmarks[i].x << "," << nearby_landmarks[i].y << "), id: "
		<< nearby_landmarks[i].id << endl;
		}
		cout << "	transformed_obs size and vector: " << transformed_obs.size() << endl;

		for (int l = 0; l < transformed_obs.size(); l++)
		{
		cout << "	(" << transformed_obs[l].x << " , " << transformed_obs[l].y
		<< ") associated lm id: " << transformed_obs[l].id << endl;
		}

		*/


		//This function will give each transformed observation an id value which will be the
		//id value of the nearest landmark 

		//Our next step is to find the Multi-Variable Gaussian of each transformed observation
		//with respect to its associated landmark.

		//We multiply this value of every transformed observation of every landmark to give us
		//the weight of the particle.

		particles[i].weight = 1.0;

		//Debugging
		/*
		cout << "	transformed obs size: " << transformed_obs.size() << endl;
		cout << "	transformed_obs[19]: " << transformed_obs[19].x << endl;
		cout << "	transformed_obs[32]: " << transformed_obs[32].x << endl;
		cout << "		for transformed observations which have been associated with landmarks:" << endl;
		*/

		if (Debugging_bool == true)
		{
			cout << i << "th particle pos: (" << particles[i].x << "," << particles[i].y << ")" << endl;
			cout << "obs, trans_obs vector size: " << observations.size() << "," << transformed_obs.size() << endl;
		}


		for (int j = 0; j < transformed_obs.size(); j++)
		{
			int obs_id = transformed_obs[j].id;
			for (int k = 0; k < nearby_landmarks.size(); k++)
			{
				if (nearby_landmarks[k].id == obs_id)
				{
					//calculate exponent
					double mu_x = nearby_landmarks[k].x;
					double mu_y = nearby_landmarks[k].y;
					double obs_x = transformed_obs[j].x;
					double obs_y = transformed_obs[j].y;

					double exponent = ((pow((obs_x - mu_x), 2)) / (2 * pow(std_landmark[0], 2))) + ((pow((obs_y - mu_y), 2)) / (2 * pow(std_landmark[1], 2)));
					double multivariate = gauss_norm*exp(-exponent);


					//Debugging
					if (Debugging_bool == true)
					{
						cout << "obs: (" << observations[j].x << "," << observations[j].y <<
							") transformed:(" << transformed_obs[j].x << "," <<
							transformed_obs[j].y << ") lm_id:" << obs_id << " lm:(" << nearby_landmarks[k].x
							<< "," << nearby_landmarks[k].y;
						cout << ") gaussian value:" << multivariate << endl;
					}

					

					particles[i].weight *= multivariate;

				}
			}
		}
		//Debugging: cout << i << "th, particle weight: " << particles[i].weight << endl;


		if (Debugging_bool ==true) cout << "	Total " << i << " particle weight: " << particles[i].weight << endl;


		//Now we have the weight of every particle
		weights.push_back(particles[i].weight);

	}

	//normalizing the weights
	double sum_of_weights = 0.0;
	//cout << "weights vector is: ";
	for (int l = 0; l < weights.size(); l++)
	{
		//cout  << weights[l] << " ";
	}

	//Debugging
	/*
	sum_of_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
	cout << "sum of weights: " << sum_of_weights << endl;

	for (int i = 0; i < weights.size(); i++)
	{
	weights[i] = weights[i] / sum_of_weights;

	}
	*/

	//std::cout << "update weights works";


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	/*std::discrete_distribution<int> d(weights.begin(), weights.end()); // Define a discrete distribution
	std::vector<Particle> new_particles; // Resampled particles holder
	std::default_random_engine gen;
	for (int i = 0; i< num_particles; i++) {
	auto index = d(gen);
	new_particles.push_back(std::move(particles[index]));
	}
	//std::cout << "resample works";
	*/

	/*Debugging:
	cout << "before resampling, weights vector: " << endl;
	for (int j = 0; j < num_particles; j++)
	{
		cout << "	" << weights[j] << endl;
	}
	*/
	if (Debugging_bool == true)
	{
		cout << "before resampling, x,y:" << endl;
		for (int j = 0; j < num_particles; j++)
		{
			cout << "	" << particles[j].x << "," << particles[j].y << endl;
		}
	}




	// Vector for new particles
	vector<Particle> new_particles(num_particles);
	
	
	
	// Use discrete distribution to return particles by weight
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	for (int i = 0; i < num_particles; ++i) {
		int number = distribution(gen);
		new_particles[i] = particles[number];
	
	}
	//weights.clear();
	particles = new_particles;

	/*Debugging*/
	if (Debugging_bool == true)
	{
		cout << "after resampling, x,y:" << endl;
		for (int j = 0; j < num_particles; j++)
		{
			cout << "	" << particles[j].x << "," << particles[j].y << endl;
		}
	}

	


	//default_random_engine gen;
	//std::discrete_distribution<int> d(weights.begin(), weights.end());
	//std::vector<Particle> new_particles;
	//
	//for (unsigned i = 0; i < num_particles; i++)
	//{
	//	int ind = d(gen);
	//	new_particles.push_back(std::move(particles[ind]));
	//}
	//particles = std::move(new_particles);

	//default_random_engine gen;

	//// Discrete distribution with all the weight values
	//discrete_distribution<> weights_prob(weights.begin(), weights.end());
	//// initialise new particle vector
	//vector<Particle> newParticles;
	//// resample particles
	//for (int i = 0; i < num_particles; ++i) {
	//	newParticles.push_back(std::move(particles[weights_prob(gen)]));
	//}
	//
	//particles = std::move(newParticles);

	//resampled_particles.clear();
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::discrete_distribution<> d(weights.begin(), weights.end());
	//for (int n = 0; n<num_particles; ++n) {
	//	resampled_particles.push_back(particles[d(gen)]);
	//}

	// clear all particles and replace with resampled ones
	//particles.erase(particles.begin(),particles.end());
	//particles = resampled_particles;


}

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

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}