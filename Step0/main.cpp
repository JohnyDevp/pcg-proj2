/**
 * @file      main.cpp
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

#include <cmath>
#include <cstdio>
#include <chrono>
#include <string>

#include "nbody.h"
#include "h5Helper.h"

/**
 * Main rotine
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv)
{
  if (argc != 7)
  {
    std::printf("Usage: %s <N> <dt> <steps> <write intesity> <input> <output>\n", argv[0]);
    std::exit(1);
  }

  // Number of particles
  const unsigned N = static_cast<unsigned>(std::stoul(argv[1]));
  // Length of time step
  const float dt = std::stof(argv[2]);
  // Number of steps
  const unsigned steps = static_cast<unsigned>(std::stoul(argv[3]));
  // Write frequency
  const unsigned writeFreq = static_cast<unsigned>(std::stoul(argv[4]));

  // Log benchmark setup
  std::printf("       NBODY GPU simulation\n"
              "N:                       %u\n"
              "dt:                      %f\n"
              "steps:                   %u\n",
              N, dt, steps);

  const std::size_t recordsCount = (writeFreq > 0) ? (steps + writeFreq - 1) / writeFreq : 0;

  Particles  particles{N};
  Velocities tmpVelocities{N};

  /********************************************************************************************************************/
  /*                                     TODO: Fill memory descriptor parameters                                      */
  /********************************************************************************************************************/

  /*
   * Caution! Create only after CPU side allocation
   * parameters:
   *                            Stride of two            Offset of the first
   *       Data pointer       consecutive elements        element in FLOATS,
   *                          in FLOATS, not bytes            not bytes
   */
  // MemDesc md(&particles.pos->x, 4, 0,
  //            &particles.pos->y, 4, 0,
  //            &particles.pos->z, 4, 0,
  //            &particles.pos->w, 4, 0,
  //            &particles.vel->x, 3, 0,
  //            &particles.vel->y, 3, 0,
  //            &particles.vel->z, 3, 0,
  //            N,
  //            recordsCount);

  MemDesc md(particles.posX,           1,                         0,
             particles.posY,           1,                         0,
             particles.posZ,           1,                         0,
             particles.velX,           1,                         0,
             particles.velY,           1,                         0,
             particles.velZ,           1,                         0,
             particles.weight,         1,                         0,
             N,
             recordsCount);

  // Initialisation of helper class and loading of input data
  H5Helper h5Helper(argv[5], argv[6], md);

  try
  {
    h5Helper.init();
    h5Helper.readParticleData();
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  }

  /********************************************************************************************************************/
  /*                                     TODO: Memory transfer CPU -> GPU                                             */
  /********************************************************************************************************************/
  // particles.copyToDevice();
  // tmpVelocities.copyToDevice();

  // Start measurement
  const auto start = std::chrono::steady_clock::now();

  // #pragma acc kernels
  for (unsigned s = 0u; s < steps; ++s)
  {
/******************************************************************************************************************/
/*                                        TODO: GPU computation                                                   */
/******************************************************************************************************************/

      calculateGravitationVelocity(particles, tmpVelocities, N, dt);

      calculateCollisionVelocity(particles, tmpVelocities, N, dt);

      updateParticles(particles, tmpVelocities, N, dt);
  }

  // End measurement
  const auto end = std::chrono::steady_clock::now();

  // Approximate simulation wall time
  const float elapsedTime = std::chrono::duration<float>(end - start).count();
  std::printf("Time: %f s\n", elapsedTime);

  /********************************************************************************************************************/
  /*                                     TODO: Memory transfer GPU -> CPU                                             */
  /********************************************************************************************************************/
  // particles.copyToHost();
  // tmpVelocities.copyToHost();

  // Compute reference center of mass on CPU
  const float4 refCenterOfMass = centerOfMassRef(md);

  std::printf("Reference center of mass: %f, %f, %f, %f\n",
              refCenterOfMass.x,
              refCenterOfMass.y,
              refCenterOfMass.z,
              refCenterOfMass.w);

  std::printf("Center of mass on GPU: %f, %f, %f, %f\n", 0.f, 0.f, 0.f, 0.f);

  // Writing final values to the file
  h5Helper.writeComFinal(refCenterOfMass);
  h5Helper.writeParticleDataFinal();
} // end of main
//----------------------------------------------------------------------------------------------------------------------
