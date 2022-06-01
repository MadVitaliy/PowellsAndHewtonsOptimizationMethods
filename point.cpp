//
// Created by Vitaliy on 30.05.2022.
//

#include "point.h"
#include <cmath>
#include <vector>
#include <cassert>
#include <iomanip>
#include <iostream>

point::point(value_type c0, value_type c1) noexcept: coords{c0, c1}
  {}

[[nodiscard]] constexpr const typename point::value_type& point::operator[](std::size_t i) const noexcept
  {
  return coords[i];
  }

[[nodiscard]] typename point::value_type& point::operator[](std::size_t i) noexcept
  {
  return coords[i];
  }

std::ostream& operator<<(std::ostream& os, const point& p)
  {
  os << p[0] << " " << p[1];
  return os;
  }

[[nodiscard]] point& operator*=(point& p, point::value_type value) noexcept
  {
  p[0] *= value;
  p[1] *= value;
  return p;
  }

[[nodiscard]] point operator*(const point& p, point::value_type value) noexcept
  {
  point result(p);
  result *= value;
  return result;
  }

[[nodiscard]] point& operator/=(point& p, point::value_type value) noexcept
  {
  p[0] /= value;
  p[1] /= value;
  return p;
  }

[[nodiscard]] point operator/(const point& p, point::value_type value) noexcept
  {
  point result(p);
  result /= value;
  return result;
  }


[[nodiscard]] point& operator+=(point& p1, const point& p2) noexcept
  {
  p1[0] += p2[0];
  p1[1] += p2[1];
  return p1;
  }

[[nodiscard]] point operator+(const point& p1, const point& p2) noexcept
  {
  point result(p1);
  result += p2;
  return result;
  }


[[nodiscard]] point& operator-=(point& p1, const point& p2) noexcept
  {
  p1[0] -= p2[0];
  p1[1] -= p2[1];
  return p1;
  }

[[nodiscard]] point operator-(const point& p1, const point& p2) noexcept
  {
  point result(p1);
  result -= p2;
  return result;
  }

[[nodiscard]] point::value_type operator*(const point& p1, const point& p2) noexcept
  {
  return p1[0] * p2[0] + p1[1] * p2[1];
  }

[[nodiscard]] sqmatrix ColumnOnRow(const point& p1, const point& p2) noexcept
  {
  return {p1[0] * p2[0], p1[0] * p2[1], p1[1] * p2[0], p1[1] * p2[1]};
  }

sqmatrix::sqmatrix(value_type c0, value_type c1, value_type c2, value_type c3) noexcept: coords{c0, c1, c2, c3}
  {}

constexpr const typename sqmatrix::value_type& sqmatrix::operator[](std::size_t i) const noexcept
  {
  return coords[i];
  }

typename sqmatrix::value_type& sqmatrix::operator[](std::size_t i) noexcept
  {
  return coords[i];
  }

constexpr const typename sqmatrix::value_type& sqmatrix::get(std::size_t i, std::size_t j) const noexcept
  {
  return coords[i * dimension + j];
  }

typename sqmatrix::value_type& sqmatrix::get(std::size_t i, std::size_t j) noexcept
  {
  return coords[i * dimension + j];
  }

std::ostream& operator<<(std::ostream& os, const sqmatrix& m)
  {
  os << std::left << std::setw(15) << m[0] << std::setw(15) << m[1] << std::endl;
  os << std::left << std::setw(15) << m[2] << std::setw(15) << m[3] << std::endl;
  return os;
  }

//2x2*2x1 = 2x1
[[nodiscard]] point operator*(const sqmatrix& m, const point& p) noexcept
  {
  return {m[0] * p[0] + m[1] * p[1], m[2] * p[0] + m[3] * p[1]};
  }

//1x2*2x2 = 1x2
[[nodiscard]] point operator*(const point& p, const sqmatrix& m) noexcept
  {
  return {p[0] * m[0] + p[1] * m[2], p[0] * m[1] + p[1] * m[3]};
  }

[[nodiscard]] sqmatrix& operator*=(sqmatrix& m, sqmatrix::value_type value) noexcept
  {
  m[0] *= value;
  m[1] *= value;
  m[2] *= value;
  m[3] *= value;
  return m;
  }

[[nodiscard]] sqmatrix operator*(const sqmatrix& m, sqmatrix::value_type value) noexcept
  {
  sqmatrix result(m);
  result *= value;
  return result;
  }

[[nodiscard]] sqmatrix& operator/=(sqmatrix& m, sqmatrix::value_type value) noexcept
  {
  m[0] /= value;
  m[1] /= value;
  m[2] /= value;
  m[3] /= value;
  return m;
  }

[[nodiscard]] sqmatrix operator/(const sqmatrix& m, sqmatrix::value_type value) noexcept
  {
  sqmatrix result(m);
  result /= value;
  return result;
  }

[[nodiscard]] sqmatrix& operator+=(sqmatrix& m1, const sqmatrix& m2) noexcept
  {
  m1[0] += m2[0];
  m1[1] += m2[1];
  m1[2] += m2[2];
  m1[3] += m2[3];
  return m1;
  }

[[nodiscard]] sqmatrix operator+(const sqmatrix& m1, const sqmatrix& m2) noexcept
  {
  sqmatrix result(m1);
  result += m2;
  return result;
  }

[[nodiscard]] sqmatrix& operator-=(sqmatrix& m1, const sqmatrix& m2) noexcept
  {
  m1[0] -= m2[0];
  m1[1] -= m2[1];
  m1[2] -= m2[2];
  m1[3] -= m2[3];
  return m1;
  }

[[nodiscard]] sqmatrix operator-(const sqmatrix& m1, const sqmatrix& m2) noexcept
  {
  sqmatrix result(m1);
  result -= m2;
  return result;
  }

[[nodiscard]] sqmatrix Invert(const sqmatrix& m)
  {
  return sqmatrix(m[3], -m[1], -m[2], m[0]) / (m[0] * m[3] - m[1] * m[2]);
  }

[[nodiscard]] bool IsEqual(const point& p1, const point& p2, double epsilon) noexcept
  {
  const auto diff = p1 - p2;
  return std::abs(diff[0]) < epsilon && std::abs(diff[1]) < epsilon;
  }

[[nodiscard]] bool IsEqual(const sqmatrix& m1, const sqmatrix& m2, double epsilon) noexcept
  {
  const auto diff = m1 - m2;
  return std::abs(diff[0]) < epsilon && std::abs(diff[1]) < epsilon && std::abs(diff[2]) < epsilon &&
         std::abs(diff[3]) < epsilon;
  }

[[nodiscard]] double LengthSqrt(const point& p1) noexcept
  {
  return p1[0] * p1[0] + p1[1] * p1[1];
  }

[[nodiscard]] double LengthSqrt(const point& p1, const point& p2) noexcept
  {
  return LengthSqrt(p1 - p2);
  }

[[nodiscard]] double Length(const point& p1) noexcept
  {
  return std::sqrt(LengthSqrt(p1));
  }

[[nodiscard]] double Length(const point& p1, const point& p2) noexcept
  {
  return std::sqrt(LengthSqrt(p1 - p2));
  }

[[nodiscard]] point Normalize(const point& p) noexcept
  {
  const auto length = Length(p);
  return {p[0] / length, p[1] / length};
  }

[[nodiscard]] point DivideInRation(const point& p1, const point& p2, double ration) noexcept
  {
  return {(p1[0] + ration * p2[0]) / (1 + ration), (p1[1] + ration * p2[1]) / (1 + ration)};
  }

[[nodiscard]] bool IsInBBox(const bbox& bb, const point& p) noexcept
  {
  return false;
  }

[[nodiscard]] edge IntersectBBox(const bbox& bb, const point& p, const point& dir)
  {
  double dp[4];
  dp[0] = (bb.p1[0] - p[0]) / dir[0];//dpx1
  dp[1] = (bb.p1[1] - p[1]) / dir[1];//dpy1
  dp[2] = (bb.p2[0] - p[0]) / dir[0];//dpx2
  dp[3] = (bb.p2[1] - p[1]) / dir[1];//dpy2

  std::vector<double> first_intersection_scalers;
  first_intersection_scalers.reserve(2);
  std::vector<double> second_intersection_scalers;
  second_intersection_scalers.reserve(2);

  for (const auto scaler: dp)
    scaler < 0 ? first_intersection_scalers.push_back(scaler) : second_intersection_scalers.push_back(scaler);

  assert(first_intersection_scalers.size() == 2 && second_intersection_scalers.size() == 2);


  const auto less_abs = [](double d1, double d2)
    { return std::abs(d1) < std::abs(d2); };

  const auto first_point_scaler = std::min(first_intersection_scalers[0], first_intersection_scalers[1], less_abs);
  const auto second_point_scaler = std::min(second_intersection_scalers[0], second_intersection_scalers[1], less_abs);

  return {p + dir * first_point_scaler, p + dir * second_point_scaler};
  }
