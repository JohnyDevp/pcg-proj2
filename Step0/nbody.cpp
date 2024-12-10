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

  posX = new float[N];
  posY = new float[N];
  posZ = new float[N];
  weight = new float[N];
  velX = new float[N];
  velY = new float[N];
  velZ = new float[N];
  ;
  vel = new float3[N];

#pragma acc enter data copyin(this[0 : 1])
#pragma acc enter data create(pos[0 : N], vel[0 : N])
#pragma acc enter data create(posX[0 : N], posY[0 : N], posZ[0 : N], weight[0 : N], velX[0 : N], velY[0 : N], velZ[0 : N])
}

/// @brief Destructor
Particles::~Particles()
{

#pragma acc exit data delete (this[0 : 1])
#pragma acc exit data delete (posX[0 : N], posY[0 : N], posZ[0 : N], weight[0 : N], velX[0 : N], velY[0 : N], velZ[0 : N])
#pragma acc exit data delete (pos[0 : N], vel[0 : N])

  delete[] posX;
  delete[] posY;
  delete[] posZ;
  delete[] weight;
  delete[] velX;
  delete[] velY;
  delete[] velZ;

  delete[] pos;
  delete[] vel;
}

/**
 * @brief Copy particles from host to device
 */
void Particles::copyToDevice()
{
#pragma acc update device(pos[0 : N], vel[0 : N])
#pragma acc update device(posX[0 : N], posY[0 : N], posZ[0 : N], weight[0 : N], velX[0 : N], velY[0 : N], velZ[0 : N])
}

/**
 * @brief Copy particles from device to host
 */
void Particles::copyToHost()
{
#pragma acc update host(pos[0 : N], vel[0 : N])
#pragma acc update host(posX[0 : N], posY[0 : N], posZ[0 : N], weight[0 : N], velX[0 : N], velY[0 : N], velZ[0 : N])
}

/**
 * @brief Constructor
 * @param N - Number of particles
 */
Velocities::Velocities(const unsigned N) : N(N)
{
  vel = new float3[N];

  x = new float[N];
  y = new float[N];
  z = new float[N];

#pragma acc enter data copyin(this[0 : 1])
#pragma acc enter data create(vel[0 : N])
#pragma acc enter data create(x[0 : N], y[0 : N], z[0 : N])
}

/// @brief Destructor
Velocities::~Velocities()
{
#pragma acc exit data delete (this[0 : 1])

#pragma acc exit data delete (vel[0 : N])
#pragma acc exit data delete (x[0 : N], y[0 : N], z[0 : N])

  delete[] x;
  delete[] y;
  delete[] z;

  delete[] vel;
}

/**
 * @brief Copy velocities from host to device
 */
void Velocities::copyToDevice()
{
#pragma acc update device(vel[0 : N])
#pragma acc update device(x[0 : N], y[0 : N], z[0 : N])
}

/**
 * @brief Copy velocities from device to host
 */
void Velocities::copyToHost()
{
#pragma acc update host(vel[0 : N])
#pragma acc update host(x[0 : N], y[0 : N], z[0 : N])
}

/*********************************************************************************************************************/

/**
 * Calculate gravitation velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */

void calculateGravitationVelocity__(Particles &p, Velocities &tmpVel, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Calculate gravitation velocity, see reference CPU version,                             */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/

  float4 *const pPos = p.pos;
  float3 *const pVel = p.vel;
  float3 *const tVel = tmpVel.vel;

  for (unsigned i = 0u; i < N; ++i)
  {
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    const float posX = pPos[i].x;
    const float posY = pPos[i].y;
    const float posZ = pPos[i].z;
    const float weight = pPos[i].w;

    for (unsigned j = 0u; j < N; ++j)
    {
      const float otherPosX = pPos[j].x;
      const float otherPosY = pPos[j].y;
      const float otherPosZ = pPos[j].z;
      const float otherWeight = pPos[j].w;

      const float dx = otherPosX - posX;
      const float dy = otherPosY - posY;
      const float dz = otherPosZ - posZ;

      const float r2 = dx * dx + dy * dy + dz * dz;
      const float r = std::sqrt(r2) + std::numeric_limits<float>::min();

      const float f = G * weight * otherWeight / r2 + std::numeric_limits<float>::min();

      newVelX += (r > COLLISION_DISTANCE) ? dx / r * f : 0.f;
      newVelY += (r > COLLISION_DISTANCE) ? dy / r * f : 0.f;
      newVelZ += (r > COLLISION_DISTANCE) ? dz / r * f : 0.f;
    }

    newVelX *= dt / weight;
    newVelY *= dt / weight;
    newVelZ *= dt / weight;

    tVel[i].x = newVelX;
    tVel[i].y = newVelY;
    tVel[i].z = newVelZ;
  }
} // end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Calculate collision velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void calculateCollisionVelocity__(Particles &p, Velocities &tmpVel, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Calculate collision velocity, see reference CPU version,                               */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/
  float4 *const pPos = p.pos;
  float3 *const pVel = p.vel;
  float3 *const tVel = tmpVel.vel;

  // # pragma omp parallel for firstprivate(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ, N, dt)
  for (unsigned i = 0u; i < N; ++i)
  {
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    const float posX = pPos[i].x;
    const float posY = pPos[i].y;
    const float posZ = pPos[i].z;
    const float weight = pPos[i].w;

    const float velX = pVel[i].x;
    const float velY = pVel[i].y;
    const float velZ = pVel[i].z;

    // #   pragma omp simd aligned(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ: dataAlignment)
    for (unsigned j = 0u; j < N; ++j)
    {
      const float otherPosX = pPos[j].x;
      const float otherPosY = pPos[j].y;
      const float otherPosZ = pPos[j].z;
      const float otherWeight = pPos[j].w;

      const float otherVelX = pVel[j].x;
      const float otherVelY = pVel[j].y;
      const float otherVelZ = pVel[j].z;

      const float dx = otherPosX - posX;
      const float dy = otherPosY - posY;
      const float dz = otherPosZ - posZ;

      const float r2 = dx * dx + dy * dy + dz * dz;
      const float r = std::sqrt(r2);

      newVelX += (r > 0.f && r < COLLISION_DISTANCE)
                     ? (((weight * velX - otherWeight * velX + 2.f * otherWeight * otherVelX) / (weight + otherWeight)) - velX)
                     : 0.f;
      newVelY += (r > 0.f && r < COLLISION_DISTANCE)
                     ? (((weight * velY - otherWeight * velY + 2.f * otherWeight * otherVelY) / (weight + otherWeight)) - velY)
                     : 0.f;
      newVelZ += (r > 0.f && r < COLLISION_DISTANCE)
                     ? (((weight * velZ - otherWeight * velZ + 2.f * otherWeight * otherVelZ) / (weight + otherWeight)) - velZ)
                     : 0.f;
    }

    tVel[i].x += newVelX;
    tVel[i].y += newVelY;
    tVel[i].z += newVelZ;
  }

} // end of calculate_collision_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * Update particles
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
void updateParticles__(Particles &p, Velocities &tmpVel, const unsigned N, float dt)
{
  /*******************************************************************************************************************/
  /*                    TODO: Update particles position and velocity, see reference CPU version,                     */
  /*                            you can use overloaded operators defined in Vec.h                                    */
  /*******************************************************************************************************************/

  float4 *const pPos = p.pos;
  float3 *const pVel = p.vel;
  float3 *const tVel = tmpVel.vel;

  // # pragma omp parallel for simd \
          firstprivate(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ, N, dt) \
          aligned(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ: dataAlignment)
  for (unsigned i = 0u; i < N; ++i)
  {
    float posX = pPos[i].x;
    float posY = pPos[i].y;
    float posZ = pPos[i].z;

    float velX = pVel[i].x;
    float velY = pVel[i].y;
    float velZ = pVel[i].z;

    const float newVelX = tVel[i].x;
    const float newVelY = tVel[i].y;
    const float newVelZ = tVel[i].z;

    velX += newVelX;
    velY += newVelY;
    velZ += newVelZ;

    posX += velX * dt;
    posY += velY * dt;
    posZ += velZ * dt;

    pPos[i].x = posX;
    pPos[i].y = posY;
    pPos[i].z = posZ;

    pVel[i].x = velX;
    pVel[i].y = velY;
    pVel[i].z = velZ;
  } // end of update_particle
  //----------------------------------------------------------------------------------------------------------------------
}

/**
 * Calculate particles center of mass
 * @param p    - particles
 * @param com  - pointer to a center of mass
 * @param lock - pointer to a user-implemented lock
 * @param N    - Number of particles
 */
void centerOfMass(Particles &p, float4 &com, int *lock, const unsigned N)
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

// ===================================================================================================
// ===================================================================================================
// ===================================================================================================

void calculateGravitationVelocity(Particles &p, Velocities &tmpVel, const unsigned N, float dt)
{
  float *const pPosX = p.posX;
  float *const pPosY = p.posY;
  float *const pPosZ = p.posZ;
  float *const pVelX = p.velX;
  float *const pVelY = p.velY;
  float *const pVelZ = p.velZ;
  float *const pWeight = p.weight;

  float *const tmpVelX = tmpVel.x;
  float *const tmpVelY = tmpVel.y;
  float *const tmpVelZ = tmpVel.z;

  // # pragma omp parallel for firstprivate(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ, N, dt)
  for (unsigned i = 0u; i < N; ++i)
  {
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    const float posX = pPosX[i];
    const float posY = pPosY[i];
    const float posZ = pPosZ[i];
    const float weight = pWeight[i];

    // #   pragma omp simd aligned(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ: dataAlignment)
    for (unsigned j = 0u; j < N; ++j)
    {
      const float otherPosX = pPosX[j];
      const float otherPosY = pPosY[j];
      const float otherPosZ = pPosZ[j];
      const float otherWeight = pWeight[j];

      const float dx = otherPosX - posX;
      const float dy = otherPosY - posY;
      const float dz = otherPosZ - posZ;

      const float r2 = dx * dx + dy * dy + dz * dz;
      const float r = std::sqrt(r2) + std::numeric_limits<float>::min();

      const float f = G * weight * otherWeight / r2 + std::numeric_limits<float>::min();

      newVelX += (r > COLLISION_DISTANCE) ? dx / r * f : 0.f;
      newVelY += (r > COLLISION_DISTANCE) ? dy / r * f : 0.f;
      newVelZ += (r > COLLISION_DISTANCE) ? dz / r * f : 0.f;
    }

    newVelX *= dt / weight;
    newVelY *= dt / weight;
    newVelZ *= dt / weight;

    tmpVelX[i] = newVelX;
    tmpVelY[i] = newVelY;
    tmpVelZ[i] = newVelZ;
  }
} // end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

void calculateCollisionVelocity(Particles &p, Velocities &tmpVel, const unsigned N, float dt)
{
  float *const pPosX = p.posX;
  float *const pPosY = p.posY;
  float *const pPosZ = p.posZ;
  float *const pVelX = p.velX;
  float *const pVelY = p.velY;
  float *const pVelZ = p.velZ;
  float *const pWeight = p.weight;

  float *const tmpVelX = tmpVel.x;
  float *const tmpVelY = tmpVel.y;
  float *const tmpVelZ = tmpVel.z;

  // # pragma omp parallel for firstprivate(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ, N, dt)
  for (unsigned i = 0u; i < N; ++i)
  {
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    const float posX = pPosX[i];
    const float posY = pPosY[i];
    const float posZ = pPosZ[i];
    const float velX = pVelX[i];
    const float velY = pVelY[i];
    const float velZ = pVelZ[i];
    const float weight = pWeight[i];

    // #   pragma omp simd aligned(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ: dataAlignment)
    for (unsigned j = 0u; j < N; ++j)
    {
      const float otherPosX = pPosX[j];
      const float otherPosY = pPosY[j];
      const float otherPosZ = pPosZ[j];
      const float otherVelX = pVelX[j];
      const float otherVelY = pVelY[j];
      const float otherVelZ = pVelZ[j];
      const float otherWeight = pWeight[j];

      const float dx = otherPosX - posX;
      const float dy = otherPosY - posY;
      const float dz = otherPosZ - posZ;

      const float r2 = dx * dx + dy * dy + dz * dz;
      const float r = std::sqrt(r2);

      newVelX += (r > 0.f && r < COLLISION_DISTANCE)
                     ? (((weight * velX - otherWeight * velX + 2.f * otherWeight * otherVelX) / (weight + otherWeight)) - velX)
                     : 0.f;
      newVelY += (r > 0.f && r < COLLISION_DISTANCE)
                     ? (((weight * velY - otherWeight * velY + 2.f * otherWeight * otherVelY) / (weight + otherWeight)) - velY)
                     : 0.f;
      newVelZ += (r > 0.f && r < COLLISION_DISTANCE)
                     ? (((weight * velZ - otherWeight * velZ + 2.f * otherWeight * otherVelZ) / (weight + otherWeight)) - velZ)
                     : 0.f;
    }

    tmpVelX[i] += newVelX;
    tmpVelY[i] += newVelY;
    tmpVelZ[i] += newVelZ;
  }
} // end of calculate_collision_velocity
//----------------------------------------------------------------------------------------------------------------------

void updateParticles(Particles &p, Velocities &tmpVel, const unsigned N, float dt)
{
  float *const pPosX = p.posX;
  float *const pPosY = p.posY;
  float *const pPosZ = p.posZ;
  float *const pVelX = p.velX;
  float *const pVelY = p.velY;
  float *const pVelZ = p.velZ;
  float *const pWeight = p.weight;

  float *const tmpVelX = tmpVel.x;
  float *const tmpVelY = tmpVel.y;
  float *const tmpVelZ = tmpVel.z;

  // # pragma omp parallel for simd \
//          firstprivate(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ, N, dt) \
//          aligned(pPosX, pPosY, pPosZ, pVelX, pVelY, pVelZ, pWeight, tmpVelX, tmpVelY, tmpVelZ: dataAlignment)
  for (unsigned i = 0u; i < N; ++i)
  {
    float posX = pPosX[i];
    float posY = pPosY[i];
    float posZ = pPosZ[i];

    float velX = pVelX[i];
    float velY = pVelY[i];
    float velZ = pVelZ[i];

    const float newVelX = tmpVelX[i];
    const float newVelY = tmpVelY[i];
    const float newVelZ = tmpVelZ[i];

    velX += newVelX;
    velY += newVelY;
    velZ += newVelZ;

    posX += velX * dt;
    posY += velY * dt;
    posZ += velZ * dt;

    pPosX[i] = posX;
    pPosY[i] = posY;
    pPosZ[i] = posZ;

    pVelX[i] = velX;
    pVelY[i] = velY;
    pVelZ[i] = velZ;
  }
} // end of update_particle
//----------------------------------------------------------------------------------------------------------------------
