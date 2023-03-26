//
// Created by Vitaliy on 28.05.2022.
//

#ifndef POWELLSMETHOD_POINT_H
#define POWELLSMETHOD_POINT_H

#include "PowellsMethod.Math/PowellsMethod.Math.API.h"

#include <cstddef>
#include <limits>
#include <ostream>

struct POWELLSMETHODMATH_API point
  {
  constexpr static std::size_t dimension = 2;
  using value_type = double;

  value_type coords[dimension] = {0.0, 0.0};

  constexpr point() = default;

  constexpr point(const point& p) = default;

  point(value_type c0, value_type c1) noexcept;

  [[nodiscard]] constexpr const value_type& operator[](std::size_t i) const noexcept;

  [[nodiscard]] value_type& operator[](std::size_t i) noexcept;
  };

struct POWELLSMETHODMATH_API sqmatrix
  {
  constexpr static std::size_t dimension = 2;
  using value_type = double;

  value_type coords[dimension * dimension] = {0.0, 0.0};

  constexpr sqmatrix() = default;

  constexpr sqmatrix(const sqmatrix& p) = default;

  sqmatrix(value_type c0, value_type c1, value_type c2, value_type c3) noexcept;

  [[nodiscard]] constexpr const value_type& operator[](std::size_t i) const noexcept;

  [[nodiscard]] value_type& operator[](std::size_t i) noexcept;

  [[nodiscard]] constexpr const value_type& get(std::size_t i, std::size_t j) const noexcept;

  [[nodiscard]] value_type& get(std::size_t i, std::size_t j) noexcept;

  };

struct POWELLSMETHODMATH_API edge
  {
  point p1;
  point p2;
  };

struct POWELLSMETHODMATH_API bbox
  {
  point p1;
  point p2;
  };


std::ostream& POWELLSMETHODMATH_API operator<<(std::ostream& os, const point& p);

[[nodiscard]] point& POWELLSMETHODMATH_API operator*=(point& p, point::value_type value) noexcept;

[[nodiscard]] point POWELLSMETHODMATH_API operator*(const point& p, point::value_type value) noexcept;

[[nodiscard]] point& POWELLSMETHODMATH_API operator/=(point& p, point::value_type value) noexcept;

[[nodiscard]] point POWELLSMETHODMATH_API operator/(const point& p, point::value_type value) noexcept;

[[nodiscard]] point& POWELLSMETHODMATH_API operator+=(point& p1, const point& p2) noexcept;

[[nodiscard]] point POWELLSMETHODMATH_API operator+(const point& p1, const point& p2) noexcept;

[[nodiscard]] point& POWELLSMETHODMATH_API operator-=(point& p1, const point& p2) noexcept;

[[nodiscard]] point POWELLSMETHODMATH_API operator-(const point& p1, const point& p2) noexcept;

//scalar
[[nodiscard]] point::value_type POWELLSMETHODMATH_API operator*(const point& p1, const point& p2) noexcept;

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API ColumnOnRow(const point& p1, const point& p2) noexcept;


//
std::ostream& POWELLSMETHODMATH_API operator<<(std::ostream& os, const sqmatrix& m);

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API operator*(const sqmatrix& m, sqmatrix::value_type value) noexcept;

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API operator*(sqmatrix::value_type value, const sqmatrix& m) noexcept;

//2x2*2x1 = 2x1
[[nodiscard]] point POWELLSMETHODMATH_API operator*(const sqmatrix& m, const point& p) noexcept;

//1x2*2x2 = 1x2
[[nodiscard]] point POWELLSMETHODMATH_API operator*(const point& p, const sqmatrix& m) noexcept;

[[nodiscard]] sqmatrix& POWELLSMETHODMATH_API operator/=(sqmatrix& m, sqmatrix::value_type value) noexcept;

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API operator/(const sqmatrix& m, sqmatrix::value_type value) noexcept;

[[nodiscard]] sqmatrix& POWELLSMETHODMATH_API operator+=(sqmatrix& m1, const sqmatrix& m2) noexcept;

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API operator+(const sqmatrix& m1, const sqmatrix& m2) noexcept;

[[nodiscard]] sqmatrix& POWELLSMETHODMATH_API operator-=(sqmatrix& m1, const sqmatrix& m2) noexcept;

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API operator-(const sqmatrix& m1, const sqmatrix& m2) noexcept;

[[nodiscard]] sqmatrix POWELLSMETHODMATH_API Invert(const sqmatrix& m);

//
[[nodiscard]] bool
POWELLSMETHODMATH_API IsEqual(const point& p1, const point& p2, double epsilon = std::numeric_limits<double>::epsilon()) noexcept;

[[nodiscard]] bool
POWELLSMETHODMATH_API IsEqual(const sqmatrix& m1, const sqmatrix& m2, double epsilon = std::numeric_limits<double>::epsilon()) noexcept;

[[nodiscard]] double  POWELLSMETHODMATH_APILengthSqrt(const point& p1) noexcept;

[[nodiscard]] double POWELLSMETHODMATH_API LengthSqrt(const point& p1, const point& p2) noexcept;

[[nodiscard]] double POWELLSMETHODMATH_API Length(const point& p1) noexcept;

[[nodiscard]] double POWELLSMETHODMATH_API Length(const point& p1, const point& p2) noexcept;

[[nodiscard]] point POWELLSMETHODMATH_API Normalize(const point& p) noexcept;

[[nodiscard]] point POWELLSMETHODMATH_API DivideInRation(const point& p1, const point& p2, double ration) noexcept;

[[nodiscard]] bool POWELLSMETHODMATH_API IsInBBox(const bbox& bb, const point& p) noexcept;

[[nodiscard]] edge POWELLSMETHODMATH_API IntersectBBox(const bbox& bb, const point& p, const point& dir);

#endif //POWELLSMETHOD_POINT_H