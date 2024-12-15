/**
 * @file      nbody.cpp
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

#include <cfloat>
#include <cmath>

#include "nbody.h"
#include "Vec.h"

/* Constants */
constexpr float G = 6.67384e-11f;
constexpr float COLLISION_DISTANCE = 0.01f;

/*********************************************************************************************************************/
/*                TODO: Fullfill Partile's and Velocitie's constructors, destructors and methods                     */
/*                                    for data copies between host and device                                        */
/*********************************************************************************************************************/

/**
 * @brief Constructor
 * @param N - Number of particles
 */
Particles::Particles(const unsigned N) : N(N)
{
  pos = new float4[N];
  vel = new float3[N];

#pragma acc enter data copyin(this[0 : 1])
#pragma acc enter data create(pos[0 : N], vel[0 : N])
}

/// @brief Destructor
Particles::~Particles()
{

#pragma acc exit data delete (this[0 : 1])
#pragma acc exit data delete (pos[0 : N], vel[0 : N])

  delete[] pos;
  delete[] vel;
}

/**
 * @brief Copy particles from host to device
 */
void Particles::copyToDevice()
{
#pragma acc update device(pos[0 : N], vel[0 : N])
}

/**
 * @brief Copy particles from device to host
 */
void Particles::copyToHost()
{
#pragma acc update host(pos[0 : N], vel[0 : N])
}

/*********************************************************************************************************************/

/**
 * Calculate velocity
 * @param pIn  - particles input
 * @param pOut - particles output
 * @param N    - Number of particles
 * @param dt   - Size of the time step
 */
void calculateVelocity(Particles &pIn, Particles &pOut, const unsigned N, float dt)
{
/*******************************************************************************************************************/
/*                    TODO: Calculate gravitation velocity, see reference CPU version,                             */
/*                            you can use overloaded operators defined in Vec.h                                    */
/*******************************************************************************************************************/
#pragma acc parallel loop present(pIn, pOut) async(1)
  for (unsigned i = 0u; i < N; ++i)
  {
    float3 newVel_first = {0, 0, 0};
    float3 newVel_second = {0, 0, 0};

    const float4 currentPos = pIn.pos[i];
    const float3 currentVel = pIn.vel[i];

    for (unsigned j = 0u; j < N; ++j)
    {
      const float4 otherPos = pIn.pos[j];
      const float3 otherVel = pIn.vel[j];

      const float4 d = otherPos - currentPos;
      const float r2 = d.x * d.x + d.y * d.y + d.z * d.z;
      const float r = std::sqrt(r2);

      const float f = G * currentPos.w * otherPos.w / r2 + std::numeric_limits<float>::min();

      const float r_w_min = r + std::numeric_limits<float>::min();
      // calculate gravitational force
      newVel_first.x += (r_w_min > COLLISION_DISTANCE) ? d.x / r_w_min * f : 0.f;
      newVel_first.y += (r_w_min > COLLISION_DISTANCE) ? d.y / r_w_min * f : 0.f;
      newVel_first.z += (r_w_min > COLLISION_DISTANCE) ? d.z / r_w_min * f : 0.f;

      // calculate collision force
      newVel_second.x += (r > 0.f && r < COLLISION_DISTANCE)
                             ? ((((currentPos.w - otherPos.w) * currentVel.x + 2.f * otherPos.w * otherVel.x) / (currentPos.w + otherPos.w)) - currentVel.x)
                             : 0.f;
      newVel_second.y += (r > 0.f && r < COLLISION_DISTANCE)
                             ? ((((currentPos.w - otherPos.w) * currentVel.y + 2.f * otherPos.w * otherVel.y) / (currentPos.w + otherPos.w)) - currentVel.y)
                             : 0.f;
      newVel_second.z += (r > 0.f && r < COLLISION_DISTANCE)
                             ? ((((currentPos.w - otherPos.w) * currentVel.z + 2.f * otherPos.w * otherVel.z) / (currentPos.w + otherPos.w)) - currentVel.z)
                             : 0.f;
    }

    newVel_first.x *= dt / currentPos.w;
    newVel_first.y *= dt / currentPos.w;
    newVel_first.z *= dt / currentPos.w;

    const float3 finvel = newVel_first + newVel_second;

    pOut.vel[i] = pIn.vel[i] + finvel;

    pOut.pos[i].x = pIn.pos[i].x + pOut.vel[i].x * dt;
    pOut.pos[i].y = pIn.pos[i].y + pOut.vel[i].y * dt;
    pOut.pos[i].z = pIn.pos[i].z + pOut.vel[i].z * dt;
  }

} // end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate particles center of mass
 * @param p         - particles
 * @param comBuffer - pointer to a center of mass buffer
 * @param N         - Number of particles
 */
void centerOfMass(Particles &p, float4 *comBuffer, const unsigned N)
{
  /********************************************************************************************************************/
  /*                 TODO: Calculate partiles center of mass inside center of mass buffer                             */
  /********************************************************************************************************************/

  const unsigned blocks = 32;
  const unsigned computeMassStream = 3;

// Parallel loop with separate reductions for each component
#pragma acc parallel loop present(comBuffer) copyin(p) async(computeMassStream)
  for (unsigned i = 0; i < blocks; i++)
  {
    // clear the buffer on the GPU
    comBuffer[i] = {0.0f, 0.0f, 0.0f, 0.0f};

    const int elementsPerBlock = (N / blocks) + 1;
    const int startIdx = elementsPerBlock * i;

#pragma acc loop seq
    for (unsigned j = startIdx; j - startIdx < elementsPerBlock; j++)
    {
      if (j >= N)
        continue;

      const float4 b = p.pos[j]; // Access the float4 position
      float dW = (comBuffer[i].w + b.w) > 0.f ? (b.w / (comBuffer[i].w + b.w)) : 0.f;

      comBuffer[i].x += (b.x - comBuffer[i].x) * dW;
      comBuffer[i].y += (b.y - comBuffer[i].y) * dW;
      comBuffer[i].z += (b.z - comBuffer[i].z) * dW;
      comBuffer[i].w += b.w;
    }
  }
  // final reduction

#pragma acc parallel loop seq async(computeMassStream)
  for (unsigned i = 1; i < blocks; i++)
  {
    const float4 b = comBuffer[i]; // Access the float4 position
    float dW = (comBuffer[0].w + b.w) > 0.f ? (b.w / (comBuffer[0].w + b.w)) : 0.f;

    comBuffer[0].x += (b.x - comBuffer[0].x) * dW;
    comBuffer[0].y += (b.y - comBuffer[0].y) * dW;
    comBuffer[0].z += (b.z - comBuffer[0].z) * dW;
    comBuffer[0].w += b.w;
  }
} // end of centerOfMass
//----------------------------------------------------------------------------------------------------------------------

/**
 * CPU implementation of the Center of Mass calculation
 * @param particles - All particles in the system
 * @param N         - Number of particles
 */
float4 centerOfMassRef(MemDesc &memDesc)
{
  float4 com{};

  for (std::size_t i{}; i < memDesc.getDataSize(); i++)
  {
    const float3 pos = {memDesc.getPosX(i), memDesc.getPosY(i), memDesc.getPosZ(i)};
    const float w = memDesc.getWeight(i);

    // Calculate the vector on the line connecting current body and most recent position of center-of-mass
    // Calculate weight ratio only if at least one particle isn't massless
    const float4 d = {pos.x - com.x,
                      pos.y - com.y,
                      pos.z - com.z,
                      ((memDesc.getWeight(i) + com.w) > 0.0f)
                          ? (memDesc.getWeight(i) / (memDesc.getWeight(i) + com.w))
                          : 0.0f};

    // Update position and weight of the center-of-mass according to the weight ration and vector
    com.x += d.x * d.w;
    com.y += d.y * d.w;
    com.z += d.z * d.w;
    com.w += w;
  }

  return com;
} // enf of centerOfMassRef
//----------------------------------------------------------------------------------------------------------------------
