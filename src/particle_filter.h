/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"
#include <random>

struct Particle {

	int id;
    // State
	double x;
	double y;
	double theta;
    // Importance weight
	double weight;
    // Id list of the landmarks the particle is associated to
    std::vector<int> associations;
    // Positions of the associated landmarks in map coordinates
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};



class ParticleFilter {
	
	// Number of particles to draw
	int num_particles; 
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;

public:
	
	// Set of current particles
    std::vector<Particle> particles;

	// Constructor
	// @param M Number of particles
	ParticleFilter() : num_particles(0), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

    enum STATE_STD_E{
        STD_X = 0,
        STD_Y,
        STD_YAW,
        STD_SIZE
    };

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
     * @param std_landmark[] Array of dimension 2 [standard deviation of x [m],
     *   standard deviation of y [m]
     * @param observations Vector of landmark observations (Measurements)
	 * @param map Map class containing map landmarks
	 */
    void updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations,
            Map &map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);
	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}

private:
    /**
     * @brief FindParticleAssociations
     * @param p particle to find the associations of
     * @param predicted predicted measurements in map coordinates
     * @param map_landmarks container of the landmarks in map coordinates
     */
    void FindParticleAssociations(Particle &p,const std::vector<LandmarkObs> predicted,
                                  Map& map_landmarks);
    static constexpr int NO_ASSOCIATION_FOUND_ID = -1;
    /**
     * @brief ComputeParticleWeight computes the weight of a given particle and its measurement
     * associations using the predicted measurements input in map coordinates
     * @param p the particle
     * @param predicted predicted measurements in map coordinates
     */
    void ComputeParticleWeight(Particle &p, const std::vector<LandmarkObs> predicted,double std_landmark[]);
    /**
     * @brief GaussianProbability : std_x and std_y of the particle filter must have been
     * init before calling this function
     * @param x measurement
     * @param y measurement
     * @param ux average on x
     * @param uy average on y
     * @return the probability
     */
    double GaussianProbability(const double x, const double y, const double ux, const double uy, double std_landmark[]);
    void NormalizeWeights();
};



#endif /* PARTICLE_FILTER_H_ */
