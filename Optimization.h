//
// Created by Vitaliy on 01.06.2022.
//

#ifndef POWELLSMETHOD_OPTIMIZATION_H
#define POWELLSMETHOD_OPTIMIZATION_H

#include <functional>

struct point;
struct sqmatrix;
struct bbox;

namespace Optimization
  {
  point Gradient(const std::function<double(const point&)>& i_fun, const point& i_p, double i_d = 0.001);

  sqmatrix Hessian(const std::function<double(const point&)>& i_fun, const point& i_p, double i_d = 0.001);

  double Derivative(const std::function<double(const point&)>& i_fun, const point& i_p, const point& i_direction,
                    double i_d = 0.000001);

  point OneDimensionalOptGolden(const std::function<double(const point&)>& i_fun, const point& i_p1, const point& i_p2,
                                double i_epsilon = 0.0001);

  point PowellsOptimization(const std::function<double(const point&)>& i_fun, const point& i_p0, const bbox& i_bb,
                            double i_epsilon = 0.0001, bool i_print = false);

  point
  NewtonsOptimization(const std::function<double(const point&)>& i_fun, const point& i_p0, double i_epsilon = 0.0001,
                      bool i_print = false);
  }


#endif //POWELLSMETHOD_OPTIMIZATION_H
