/**
 * @file      Vec.h
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

#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <tuple>

// Forward declarations
template<typename T> struct Vec3;
template<typename T> struct Vec4;

// Type aliases
using float3 = Vec3<float>;
using float4 = Vec4<float>;

/**
 * @brief 3D vector
 * @tparam T Type of vector components
 */
template<typename T>
struct Vec3
{
  Vec3() = default;
  constexpr Vec3(T value);
  constexpr Vec3(T x, T y, T z);
  Vec3(const Vec3& other) = default;
  Vec3(Vec3&& other) = default;
  ~Vec3() = default;

  Vec3& operator=(const Vec3& other) = default;
  Vec3& operator=(Vec3&& other) = default;

  constexpr Vec3& operator+=(const Vec3& other);
  constexpr Vec3& operator-=(const Vec3& other);
  constexpr Vec3& operator*=(const Vec3& other);
  constexpr Vec3& operator/=(const Vec3& other);

  constexpr Vec3& operator+=(T value);
  constexpr Vec3& operator-=(T value);
  constexpr Vec3& operator*=(T value);
  constexpr Vec3& operator/=(T value);

  constexpr T abs() const;

  T x;
  T y;
  T z;
};

/**
 * @brief 4D vector
 * @tparam T Type of vector components
 */
template<typename T>
struct Vec4
{
  Vec4() = default;
  constexpr Vec4(T value);
  constexpr Vec4(T x, T y, T z, T w);
  explicit constexpr Vec4(const Vec3<T>& other, T w);
  Vec4(const Vec4& other) = default;
  Vec4(Vec4&& other) = default;

  Vec4& operator=(const Vec4& other) = default;
  Vec4& operator=(Vec4&& other) = default;

  constexpr Vec4& operator+=(const Vec4& other);
  constexpr Vec4& operator-=(const Vec4& other);
  constexpr Vec4& operator*=(const Vec4& other);
  constexpr Vec4& operator/=(const Vec4& other);

  constexpr Vec4& operator+=(T value);
  constexpr Vec4& operator-=(T value);
  constexpr Vec4& operator*=(T value);
  constexpr Vec4& operator/=(T value);

  template<typename... Us>
  constexpr std::tuple<Us...> decompose() const;
  constexpr T abs() const;

  T x;
  T y;
  T z;
  T w;
};

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>::Vec3(T value)
: x(value), y(value), z(value)
{}

template<typename T>
constexpr Vec3<T>::Vec3(T x, T y, T z)
: x(x), y(y), z(z)
{}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator+=(const Vec3& other)
{
  x += other.x;
  y += other.y;
  z += other.z;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator-=(const Vec3& other)
{
  x -= other.x;
  y -= other.y;
  z -= other.z;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator*=(const Vec3& other)
{
  x *= other.x;
  y *= other.y;
  z *= other.z;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator/=(const Vec3& other)
{
  x /= other.x;
  y /= other.y;
  z /= other.z;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator+=(T value)
{
  x += value;
  y += value;
  z += value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator-=(T value)
{
  x -= value;
  y -= value;
  z -= value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator*=(T value)
{
  x *= value;
  y *= value;
  z *= value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T>& Vec3<T>::operator/=(T value)
{
  x /= value;
  y /= value;
  z /= value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator+(const Vec3<T>& lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) += rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator-(const Vec3<T>& lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) -= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator*(const Vec3<T>& lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) *= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator/(const Vec3<T>& lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) /= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator+(const Vec3<T>& lhs, T value)
{
  return Vec3(lhs) += value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator-(const Vec3<T>& lhs, T value)
{
  return Vec3(lhs) -= value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator*(const Vec3<T>& lhs, T value)
{
  return Vec3(lhs) *= value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator/(const Vec3<T>& lhs, T value)
{
  return Vec3(lhs) /= value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator+(T lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) += rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator-(T lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) -= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator*(T lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) *= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec3<T> operator/(T lhs, const Vec3<T>& rhs)
{
  return Vec3(lhs) /= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr T Vec3<T>::abs() const
{
  return std::sqrt(x * x + y * y + z * z);
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>::Vec4(T value)
: x(value), y(value), z(value), w(value)
{}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>::Vec4(T x, T y, T z, T w)
: x(x), y(y), z(z), w(w)
{}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>::Vec4(const Vec3<T>& other, T w)
: x(other.x), y(other.y), z(other.z), w(w)
{}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator+=(const Vec4& other)
{
  x += other.x;
  y += other.y;
  z += other.z;
  w += other.w;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator-=(const Vec4& other)
{
  x -= other.x;
  y -= other.y;
  z -= other.z;
  w -= other.w;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator*=(const Vec4& other)
{
  x *= other.x;
  y *= other.y;
  z *= other.z;
  w *= other.w;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator/=(const Vec4& other)
{
  x /= other.x;
  y /= other.y;
  z /= other.z;
  w /= other.w;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator+=(T value)
{
  x += value;
  y += value;
  z += value;
  w += value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator-=(T value)
{
  x -= value;
  y -= value;
  z -= value;
  w -= value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator*=(T value)
{
  x *= value;
  y *= value;
  z *= value;
  w *= value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T>& Vec4<T>::operator/=(T value)
{
  x /= value;
  y /= value;
  z /= value;
  w /= value;
  return *this;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator+(const Vec4<T>& lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) += rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator-(const Vec4<T>& lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) -= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator*(const Vec4<T>& lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) *= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator/(const Vec4<T>& lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) /= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator+(const Vec4<T>& lhs, T value)
{
  return Vec4(lhs) += value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator-(const Vec4<T>& lhs, T value)
{
  return Vec4(lhs) -= value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator*(const Vec4<T>& lhs, T value)
{
  return Vec4(lhs) *= value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator/(const Vec4<T>& lhs, T value)
{
  return Vec4(lhs) /= value;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator+(T lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) += rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator-(T lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) -= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator*(T lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) *= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr Vec4<T> operator/(T lhs, const Vec4<T>& rhs)
{
  return Vec4(lhs) /= rhs;
}

#pragma acc routine seq
template<typename T>
constexpr T Vec4<T>::abs() const
{
  return std::sqrt(x * x + y * y + z * z + w * w);
}

#endif /* VEC_H */
