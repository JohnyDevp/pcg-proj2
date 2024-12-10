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


void calculateGravitationVelocity(Particles &p, Particles &pOut, const unsigned N, float dt)
{
/*******************************************************************************************************************/
/*                    TODO: Calculate gravitation velocity, see reference CPU version,                             */
/*                            you can use overloaded operators defined in Vec.h                                    */
/*******************************************************************************************************************/
// #pragma acc parallel loop present(p, tmpVel)
  for (unsigned i = 0u; i < N; ++i)
  {
    float3 newVel = {0, 0, 0};

    const float4 currentPos = p.pos[i];

    for (unsigned j = 0u; j < N; ++j)
    {
      const float4 otherPos = p.pos[j];
      const float4 d = otherPos - currentPos;

      const float r = d.abs() + std::numeric_limits<float>::min();

      const float fr = (G * dt * otherPos.w) / (r * r * r + std::numeric_limits<float>::min());
      newVel.x += (r > COLLISION_DISTANCE) ? (d.x * fr) : 0.f;
      newVel.y += (r > COLLISION_DISTANCE) ? (d.y * fr) : 0.f;
      newVel.z += (r > COLLISION_DISTANCE) ? (d.z * fr) : 0.f;
    }

    pOut.vel[i] = newVel;
  }
} // end of calculate_gravitation_velocity
void calculateCollisionVelocity(Particles &p, Particles &pOut, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Calculate collision velocity, see reference CPU version,                               */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/

// #pragma acc parallel loop present(p, tmpVel)
  for (unsigned i = 0u; i < N; ++i)
  {
    float3 newVel = {0, 0, 0};

    float4 currentPos = p.pos[i];
    float3 currentVel = p.vel[i];

    for (unsigned j = 0u; j < N; ++j)
    {
      const float4 otherPos = p.pos[j];
      const float3 otherVel = p.vel[j];

      const float4 d = otherPos - currentPos;
      const float r = d.abs();

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
    pOut.vel[i] += newVel;
  }

} // end of calculate_collision_velocity
void updateParticles(Particles &p, Particles &pOut, const unsigned N, float dt)
{
/*******************************************************************************************************************/
/*                    TODO: Update particles position and velocity, see reference CPU version,                     */
/*                            you can use overloaded operators defined in Vec.h                                    */
/*******************************************************************************************************************/
// #pragma acc parallel loop present(p, tmpVel)
  for (unsigned i = 0u; i < N; ++i)
  {
    p.vel[i] += pOut.vel[i];

    p.pos[i].x += p.vel[i].x * dt;
    p.pos[i].y += p.vel[i].y * dt;
    p.pos[i].z += p.vel[i].z * dt;
  } // end of update_particle
  //----------------------------------------------------------------------------------------------------------------------
}


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
calculateGravitationVelocity(pIn, pOut, N, dt);
calculateCollisionVelocity(pIn, pOut, N, dt);
updateParticles(pIn, pOut, N, dt);

  // for (unsigned i = 0u; i < N; ++i)
  // {
  //   float3 newVel = {0, 0, 0};

  //   const float4 currentPos = pIn.pos[i];
  //   float3 currentVel = pIn.vel[i];

  //   for (unsigned j = 0u; j < N; ++j)
  //   {
  //     const float4 otherPos = pIn.pos[j];
  //     const float4 d = otherPos - currentPos;

  //     const float r = d.abs() + std::numeric_limits<float>::min();

  //     const float fr = (G * dt * otherPos.w) / (r * r * r + std::numeric_limits<float>::min());
  //     newVel.x += (r > COLLISION_DISTANCE) ? (d.x * fr) : 0.f;
  //     newVel.y += (r > COLLISION_DISTANCE) ? (d.y * fr) : 0.f;
  //     newVel.z += (r > COLLISION_DISTANCE) ? (d.z * fr) : 0.f;
      
  //     // ==============================
  //     const float3 otherVel = pIn.vel[j];
  //     const float weight_diff = (2.f * otherPos.w) / (currentPos.w - otherPos.w);

  //     newVel.x += (r > 0.f && r < COLLISION_DISTANCE)
  //                     ? (otherVel.x * weight_diff)
  //                     : 0.f;
  //     newVel.y += (r > 0.f && r < COLLISION_DISTANCE)
  //                     ? (otherVel.y * weight_diff)
  //                     : 0.f;
  //     newVel.z += (r > 0.f && r < COLLISION_DISTANCE)
  //                     ? (otherVel.z * weight_diff)
  //                     : 0.f;
  //   }

  //   pOut.vel[i] = pIn.vel[i] + newVel;
  //   pOut.pos[i].x = pIn.pos[i].x + pIn.vel[i].x * dt;
  //   pOut.pos[i].y = pIn.pos[i].y + pIn.vel[i].y * dt;
  //   pOut.pos[i].z = pIn.pos[i].z + pIn.vel[i].z * dt;
  // }

  // =================================================================================================================

} // end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------



/**
 * Calculate particles center of mass
 * @param p    - particles
 * @param com  - pointer to a center of mass
 * @param lock - pointer to a user-implemented lock
 * @param N    - Number of particles
 */
void centerOfMass(Particles &p, float4 *comBuffer, const unsigned N)
{

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
