/**
 * @file      nbody.cpp
 *
 * @author    Name Surname \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xlogin00@fit.vutbr.cz
 *
 * @brief     PCG Assignment 2
 *
 * @version   2023
 *
 * @date      04 October   2023, 09:00 (created) \n
 */

#include <cfloat>
#include <cmath>

#include "nbody.h"

/* Constants */
constexpr float G                  = 6.67384e-11f;
constexpr float COLLISION_DISTANCE = 0.01f;

/*********************************************************************************************************************/
/*                TODO: Fullfill Partile's and Velocitie's constructors, destructors and methods                     */
/*                                    for data copies between host and device                                        */
/*********************************************************************************************************************/

/**
 * @brief Constructor
 * @param N - Number of particles
 */
Particles::Particles(const unsigned N)
{

}

/// @brief Destructor
Particles::~Particles()
{

}

/**
 * @brief Copy particles from host to device
 */
void Particles::copyToDevice()
{

}

/**
 * @brief Copy particles from device to host
 */
void Particles::copyToHost()
{

}

/**
 * @brief Constructor
 * @param N - Number of particles
 */
Velocities::Velocities(const unsigned N)
{

}

/// @brief Destructor
Velocities::~Velocities()
{

}

/**
 * @brief Copy velocities from host to device
 */
void Velocities::copyToDevice()
{

}

/**
 * @brief Copy velocities from device to host
 */
void Velocities::copyToHost()
{

}

/*********************************************************************************************************************/

/**
 * Calculate gravitation velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void calculateGravitationVelocity(Particles& p, Velocities& tmpVel, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Calculate gravitation velocity, see reference CPU version,                             */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/


}// end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate collision velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void calculateCollisionVelocity(Particles& p, Velocities& tmpVel, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Calculate collision velocity, see reference CPU version,                               */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/


}// end of calculate_collision_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Update particles
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void updateParticles(Particles& p, Velocities& tmpVel, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Update particles position and velocity, see reference CPU version,                     */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/


}// end of update_particle
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate particles center of mass
 * @param p    - particles
 * @param com  - pointer to a center of mass
 * @param lock - pointer to a user-implemented lock
 * @param N    - Number of particles
 */
void centerOfMass(Particles& p, float4& com, int* lock, const unsigned N)
{

}// end of centerOfMass
//----------------------------------------------------------------------------------------------------------------------

/**
 * CPU implementation of the Center of Mass calculation
 * @param particles - All particles in the system
 * @param N         - Number of particles
 */
float4 centerOfMassRef(MemDesc& memDesc)
{
  float4 com{};

  for (std::size_t i{}; i < memDesc.getDataSize(); i++)
  {
    const float3 pos = {memDesc.getPosX(i), memDesc.getPosY(i), memDesc.getPosZ(i)};
    const float  w   = memDesc.getWeight(i);

    // Calculate the vector on the line connecting current body and most recent position of center-of-mass
    // Calculate weight ratio only if at least one particle isn't massless
    const float4 d = {pos.x - com.x,
                      pos.y - com.y,
                      pos.z - com.z,
                      ((memDesc.getWeight(i) + com.w) > 0.0f)
                        ? ( memDesc.getWeight(i) / (memDesc.getWeight(i) + com.w))
                        : 0.0f};

    // Update position and weight of the center-of-mass according to the weight ration and vector
    com.x += d.x * d.w;
    com.y += d.y * d.w;
    com.z += d.z * d.w;
    com.w += w;
  }

  return com;
}// enf of centerOfMassRef
//----------------------------------------------------------------------------------------------------------------------
