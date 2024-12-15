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

  Particles particles[2]{Particles{N}, Particles{N}};

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
  MemDesc md(&particles[0].pos->x, 4, 0,
             &particles[0].pos->y, 4, 0,
             &particles[0].pos->z, 4, 0,
             &particles[0].vel->x, 3, 0,
             &particles[0].vel->y, 3, 0,
             &particles[0].vel->z, 3, 0,
             &particles[0].pos->w, 4, 0,
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
  /*                   TODO: Allocate memory for center of mass buffer. Remember to clear it.                         */
  /********************************************************************************************************************/

  unsigned blocks = 32;
  float4 *comBuffer = (float4 *)malloc(blocks * sizeof(float4));
  for (unsigned i = 0u; i < blocks; ++i)
  {
    comBuffer[i] = {0.0f, 0.0f, 0.0f, 0.0f};
  }
#pragma acc enter data copyin(comBuffer[0 : blocks])

  /********************************************************************************************************************/
  /*                                      TODO: Set openacc stream ids                                                */
  /********************************************************************************************************************/

  /********************************************************************************************************************/
  /*                                     TODO: Memory transfer CPU -> GPU                                             */
  /********************************************************************************************************************/
  for (unsigned i = 0u; i < N; ++i)
  {
    particles[1].pos[i] = particles[0].pos[i];
    particles[1].vel[i] = particles[0].vel[i];
  }

  particles[0].copyToDevice();
  particles[1].copyToDevice();

  // Lambda for checking if we should write current step to the file
  auto shouldWrite = [writeFreq](unsigned s) -> bool
  {
    return writeFreq > 0u && (s % writeFreq == 0u);
  };

  // Lamda for getting record number
  auto getRecordNum = [writeFreq](unsigned s) -> unsigned
  {
    return s / writeFreq;
  };

  // Start measurement
  const auto start = std::chrono::steady_clock::now();

  /********************************************************************************************************************/
  /*            TODO: Edit the loop to work asynchronously and overlap computation with data transfers.               */
  /*                  Use shouldWrite lambda to determine if data should be outputted to file.                        */
  /*                           if (shouldWrite(s, writeFreq)) { ... }                                                 */
  /*                        Use getRecordNum lambda to get the record number.                                         */
  /********************************************************************************************************************/
  for (unsigned s = 0u; s < steps; ++s)
  {
    const unsigned srcIdx = s % 2;       // source particles index
    const unsigned dstIdx = (s + 1) % 2; // destination particles index

    /******************************************************************************************************************/
    /*                                        TODO: GPU computation                                                   */
    /******************************************************************************************************************/
    calculateVelocity(particles[srcIdx], particles[dstIdx], N, dt);

    if (shouldWrite(s))
    {
      auto recordNum = getRecordNum(s);

      particles[0].copyToHost();
      particles[1].copyToHost();

      h5Helper.writeParticleData(recordNum);

      centerOfMass(particles[dstIdx], comBuffer, N);
      
#pragma acc update self(comBuffer[0])

      float4 comFinal = {comBuffer[0].x, comBuffer[0].y, comBuffer[0].z, comBuffer[0].w};
      h5Helper.writeCom(comFinal, recordNum);
    }
  }

  const unsigned resIdx = steps % 2; // result particles index

  /********************************************************************************************************************/
  /*                          TODO: Invocation of center of mass kernel, do not forget to add                         */
  /*                              additional synchronization and set appropriate stream                               */
  /********************************************************************************************************************/

  // final reduction
  centerOfMass(particles[resIdx], comBuffer, N);

  // End measurement
  const auto end = std::chrono::steady_clock::now();

  // Approximate simulation wall time
  const float elapsedTime = std::chrono::duration<float>(end - start).count();
  std::printf("Time: %f s\n", elapsedTime);

/********************************************************************************************************************/
/*                                     TODO: Memory transfer GPU -> CPU                                             */
/********************************************************************************************************************/
#pragma acc update self(comBuffer[0])
  float4 comFinal = {comBuffer[0].x, comBuffer[0].y, comBuffer[0].z, comBuffer[0].w};

  particles[0].copyToHost();
  particles[1].copyToHost();
  // if the result end up in the second array, copy it back to the first one
  // because the md is set to the first array
  if (resIdx == 1)
  {
    // Copy result from device to host
    for (unsigned i = 0u; i < N; ++i)
    {
      particles[0].pos[i] = particles[1].pos[i];
      particles[0].vel[i] = particles[1].vel[i];
    }
  }

  // Compute reference center of mass on CPU
  const float4 refCenterOfMass = centerOfMassRef(md);

  std::printf("Reference center of mass: %f, %f, %f, %f\n",
              refCenterOfMass.x,
              refCenterOfMass.y,
              refCenterOfMass.z,
              refCenterOfMass.w);

  std::printf("Center of mass on GPU: %f, %f, %f, %f\n",
              comFinal.x,
              comFinal.y,
              comFinal.z,
              comFinal.w);

  // Writing final values to the file
  h5Helper.writeComFinal(comFinal);
  h5Helper.writeParticleDataFinal();

/********************************************************************************************************************/
/*                                TODO: Free center of mass buffer memory                                           */
/********************************************************************************************************************/
#pragma acc exit data delete (comBuffer[0 : blocks])
  free(comBuffer);
} // end of main
//----------------------------------------------------------------------------------------------------------------------
