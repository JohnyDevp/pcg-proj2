/**
 * @file      nbody.h
 *
 * @author    Jan Holan \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xholan11@fit.vutbr.cz
 *
 * @brief     PCG Assignment 2
 *
 * @version   2023
 *
 * @date      9 December   2024 \n
 */

#ifndef NBODY_H
#define NBODY_H

#include "h5Helper.h"
#include "Vec.h"

/**
 * @brief Particles data structure
 */
struct Particles
{
  /// @brief Default constructor not allowed
  Particles() = delete;

  /**
   * @brief Constructor
   * @param N - Number of particles
   */
  Particles(const unsigned N);

  /// @brief Copy constructor not allowed
  Particles(const Particles &) = delete;

  /// @brief Move constructor not allowed
  Particles(Particles &&) = delete;

  /// @brief Destructor
  ~Particles();

  /// @brief Copy assignment operator not allowed
  Particles &operator=(const Particles &) = delete;

  /// @brief Move assignment operator not allowed
  Particles &operator=(Particles &&) = delete;

  /**
   * @brief Copy particles from host to device
   */
  void copyToDevice();

  /**
   * @brief Copy particles from device to host
   */
  void copyToHost();

  /********************************************************************************************************************/
  /* TODO: Particles data structure optimized for use on GPU. Use float3 and float4 structures defined in file Vec.h  */
  /********************************************************************************************************************/
  float4 *pos;
  float3 *vel;
  int N;
};

/**
 * @brief Velocities data structure (to be used as buffer for partial results)
 */
struct Velocities
{
  /// @brief Default constructor not allowed
  Velocities() = delete;

  /**
   * @brief Constructor
   * @param N - Number of particles
   */
  Velocities(const unsigned N);

  /// @brief Copy constructor not allowed
  Velocities(const Velocities &) = delete;

  /// @brief Move constructor not allowed
  Velocities(Velocities &&) = delete;

  /// @brief Destructor
  ~Velocities();

  /// @brief Copy assignment operator not allowed
  Velocities &operator=(const Velocities &) = delete;

  /// @brief Move assignment operator not allowed
  Velocities &operator=(Velocities &&) = delete;

  /**
   * @brief Copy velocities from host to device
   */
  void copyToDevice();

  /**
   * @brief Copy velocities from device to host
   */
  void copyToHost();

  /********************************************************************************************************************/
  /* TODO: Velocities data structure optimized for use on GPU. Use float3 and float4 structures defined in file Vec.h */
  /********************************************************************************************************************/
  float3 *vel;
  int N;
};

/**
 * Calculate gravitation velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void calculateGravitationVelocity(Particles &p,
                                  Velocities &tmpVel,
                                  const unsigned N,
                                  float dt);

/**
 * Calculate collision velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void calculateCollisionVelocity(Particles &p,
                                Velocities &tmpVel,
                                const unsigned N,
                                float dt);

/**
 * Update particles
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void updateParticles(Particles &p,
                     Velocities &tmpVel,
                     const unsigned N,
                     float dt);

/**
 * Calculate particles center of mass
 * @param p    - particles
 * @param com  - pointer to a center of mass
 * @param lock - pointer to a user-implemented lock
 * @param N    - Number of particles
 */
void centerOfMass(Particles &p,
                  float4 &com,
                  int *lock,
                  const unsigned N);

/**
 * CPU implementation of the Center of Mass calculation
 * @param memDesc - Memory descriptor of particle data on CPU side
 */
float4 centerOfMassRef(MemDesc &memDesc);

#endif /* NBODY_H */
