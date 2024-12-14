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
#pragma acc parallel loop present(pIn, pOut)
  for (unsigned i = 0u; i < N; ++i)
  {
    float3 newVel = {0, 0, 0};

    const float4 currentPos = pIn.pos[i];
    const float3 currentVel = pIn.vel[i];

    for (unsigned j = 0u; j < N; ++j)
    {
      const float4 otherPos = pIn.pos[j];
      const float3 otherVel = pIn.vel[j];

      const float4 d = otherPos - currentPos;
      float r = d.abs();

      // calculate gravitational force
      const float fr = (G * dt * otherPos.w) / (r * r * r + std::numeric_limits<float>::min());
      newVel.x += (r > COLLISION_DISTANCE) ? (d.x * fr) : 0.f;
      newVel.y += (r > COLLISION_DISTANCE) ? (d.y * fr) : 0.f;
      newVel.z += (r > COLLISION_DISTANCE) ? (d.z * fr) : 0.f;

      // calculate collision force
      newVel.x += (r > 0.f && r < COLLISION_DISTANCE)
                      ? ((((currentPos.w - otherPos.w) * currentVel.x + 2.f * otherPos.w * otherVel.x) / (currentPos.w + otherPos.w)) - currentVel.x)
                      : 0.f;
      newVel.y += (r > 0.f && r < COLLISION_DISTANCE)
                      ? ((((currentPos.w - otherPos.w) * currentVel.y + 2.f * otherPos.w * otherVel.y) / (currentPos.w + otherPos.w)) - currentVel.y)
                      : 0.f;
      newVel.z += (r > 0.f && r < COLLISION_DISTANCE)
                      ? ((((currentPos.w - otherPos.w) * currentVel.z + 2.f * otherPos.w * otherVel.z) / (currentPos.w + otherPos.w)) - currentVel.z)
                      : 0.f;
    }

    pOut.vel[i] = pIn.vel[i] + newVel;

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
    // Separate reduction variables for each component
    float comX = 0.0f;
    float comY = 0.0f;
    float comZ = 0.0f;
    float comW = 0.0f;

    // Parallel loop with separate reductions for each component
    // #pragma acc parallel loop firstprivate(comX,comY,comZ,comW) present(p) reduction(+:comX, comY, comZ, comW) present(p)
    #pragma acc kernels
    for (unsigned i = 0; i < N; i++) {
        const float4 b = p.pos[i]; // Access the float4 position
        float dW = (comW + b.w) > 0.f ? (b.w / (comW + b.w)) : 0.f;

        comX += (b.x - comX) * dW;
        comY += (b.y - comY) * dW;
        comZ += (b.z - comZ) * dW;
        comW += b.w;
    }

    // Atomic updates to combine the local reductions into the global center of mass buffer
    #pragma acc atomic update
    comBuffer->x += comX;

    #pragma acc atomic update
    comBuffer->y += comY;

    #pragma acc atomic update
    comBuffer->z += comZ;

    #pragma acc atomic update
    comBuffer->w += comW;

// #pragma acc data present(p) copy(comBuffer)
//   {
//     // Declare local reduction variables for x, y, z, and w
//     float x_sum = 0.0f, y_sum = 0.0f, z_sum = 0.0f, w_sum = 0.0f;

//     // float4 localCom = {0.0f, 0.0f, 0.0f, 0.0f};

//     // Parallel reduction for all particles
// #pragma acc parallel loop reduction(+ : x_sum, y_sum, z_sum, w_sum)
//     for (std::size_t i = 0; i < N; i++)
//     {
//       const float4 pos = {p.pos[i].x, p.pos[i].y, p.pos[i].z, p.pos[i].w};
//       const float4 d = {pos.x - x_sum,
//                         pos.y - y_sum,
//                         pos.z - z_sum,
//                         ((pos.w + w_sum) > 0.0f)
//                             ? (pos.w / (pos.w + w_sum))
//                             : 0.0f};
//       x_sum += d.x * d.w;
//       y_sum += d.y * d.w;
//       z_sum += d.z * d.w;
//       w_sum += pos.w;
//     }

// // Write back the reduced results to comBuffer (single-threaded section)
// #pragma acc serial present(comBuffer)
//     {
//       comBuffer->x += x_sum;
//       comBuffer->y += y_sum;
//       comBuffer->z += z_sum;
//       comBuffer->w += w_sum;
//     }
//   }

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
