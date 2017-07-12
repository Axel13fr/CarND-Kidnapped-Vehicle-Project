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
    num_particles = 1000;
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    std::default_random_engine gen;
    std::normal_distribution<> gauss_x(x,std[STD_X]);
    std::normal_distribution<> gauss_y(y,std[STD_Y]);
    std::normal_distribution<> gauss_yaw(theta,std[STD_YAW]);

    for(int i = 0; i < num_particles; i++){
        Particle p{i,x + gauss_x(gen),y + gauss_y(gen),theta + gauss_yaw(gen),1.0,{},{},{}};
        particles.push_back(p);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    std::default_random_engine gen;
    std::normal_distribution<> gauss_x(0,std_pos[STD_X]);
    std::normal_distribution<> gauss_y(0,std_pos[STD_Y]);
    std::normal_distribution<> gauss_yaw(0,std_pos[STD_YAW]);

    for(auto& p : particles){
        // First update particles based on elapsed time
        const auto yaw = p.theta;

        //avoid division by zero
        if (fabs(yaw_rate) > 0.001) {
            p.x = p.x + velocity/yaw_rate * ( sin (yaw + yaw_rate*delta_t) - sin(yaw));
            p.y = p.y + velocity/yaw_rate * ( cos(yaw) - cos(yaw+yaw_rate*delta_t) );
            p.theta += yaw_rate*delta_t;
        }
        else {
            p.x += velocity*delta_t*cos(yaw);
            p.y += velocity*delta_t*sin(yaw);
            //p.theta += 0;
        }

        // Add process noise
        p.x += gauss_x(gen);
        p.y += gauss_y(gen);
        p.theta += gauss_yaw(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    for(auto obs : observations){

    }

}

void ParticleFilter::FindParticleAssociations(Particle& p,const std::vector<LandmarkObs> predicted,Map& map_landmarks)
{
    // Twice the sensor range
    double dist = 200;
    for(const auto& lmark : predicted){
        LandmarkObs closest;
        for(const auto& map_lmark : map_landmarks.landmark_list){
            // Save closest map land mark
            auto new_dist = min(dist,lmark.distToMapLandMark(map_lmark));
            if(new_dist < dist){
                closest = {map_lmark.id_i,static_cast<double>(map_lmark.x_f),
                                          static_cast<double>(map_lmark.y_f)};
                dist = new_dist;
            }
        }
        p.associations.push_back(closest.id);
        p.sense_x.push_back(closest.x);
        p.sense_y.push_back(closest.y);
    }
}

double ParticleFilter::GaussianProbability(const double x,const double y,const double ux,const double uy)
{
    const double factor = 1/(2.0*M_PI*std_x*std_y);
    const auto x_d2 = (x - ux)*(x - ux);
    const auto y_d2 = (y - uy)*(y - uy);

    return factor*exp(-( x_d2/(2*std_x*std_x) + y_d2/(2*std_y*std_y)) );
}

void ParticleFilter::ComputeParticleWeight(Particle& p,const std::vector<LandmarkObs> predicted)
{
    p.weight = 1;
    // for each associated map points
    for(size_t i = 0; i < predicted.size() ; i++){
        p.weight *= GaussianProbability(p.sense_x[i], // map associated landmark
                                        p.sense_y[i],
                                        predicted[i].x, // measurement
                                        predicted[i].y);
    }
}

void ParticleFilter::NormalizeWeights()
{
    double sum = 0;
    for(const auto& p : particles){
        sum += p.weight;
    }
    weights.clear();
    for(auto& p : particles){
        p.weight /= sum;
        // Copy of the weights for reuse in resampling: optimization
        weights.push_back(p.weight);
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map& map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system.
    std_x = std_landmark[0];
    std_y = std_landmark[1];

    for(auto p : particles){

        std::vector<LandmarkObs> predicted;
        // 1. Transform the observation from veh coordinates to map coordinates
        for(const auto& obs_veh : observations){
            // observation to map coordinates from particle (vehicle) coordinates
            // O is the center of the map coord
            // Obs is the observation position
            // Part is the particle position
            LandmarkObs obs_map;
            // -->           -->            -->
            // OObs in map = OPart in map + PartObs in map
            //           -->                ---->
            //	       = OPart in map + Rot*PartObs in part
            // with Rot = [ cos(theta) -sin(theta) ]
            //	          [ sin(theta)  cos(theta) ]
            obs_map.x = p.x + cos(p.theta)*obs_veh.x - sin(p.theta)*obs_veh.y;
            obs_map.y = p.y + sin(p.theta)*obs_veh.x + cos(p.theta)*obs_veh.y;
            predicted.push_back(obs_map);
        }

        // 2. Find closest predicted measurement
        FindParticleAssociations(p,predicted,map_landmarks);

        //3. Compute weight
        ComputeParticleWeight(p,predicted);

    }

    //4. Normalize weights to finish the update
    //NormalizeWeights(); //-->  Not needed yhen using discrete_distribution !

}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::default_random_engine gen;
    // will draw an index in the range [0, number of particles], each with a probability
    // corresponding to the weights
    std::discrete_distribution<size_t> rand(weights.begin(),weights.end());

    // Draw particles randomly and save them
    std::vector<Particle> resampled;
    for(size_t i = 0 ; i < particles.size() ; i++){
        resampled.push_back(particles[rand(gen)]);
    }
    // Replace old particles with new resampled selection
    particles = resampled;
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
