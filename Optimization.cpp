//
// Created by Vitaliy on 01.06.2022.
//

#include "Optimization.h"

#include "point.h"
#include <iostream>

namespace Optimization
  {


  point Gradient(const std::function<double(const point&)>& i_fun, const point& i_p, double i_d)
    {
    const auto fpx1 = i_fun({i_p[0] + i_d, i_p[1]});
    const auto fmx1 = i_fun({i_p[0] - i_d, i_p[1]});

    const auto fpx2 = i_fun({i_p[0], i_p[1] + i_d});
    const auto fmx2 = i_fun({i_p[0], i_p[1] - i_d});

    return {(fpx1 - fmx1) / (2 * i_d), (fpx2 - fmx2) / (2 * i_d)};
    }


  sqmatrix Hessian(const std::function<double(const point&)>& i_fun, const point& i_p, double i_d)
    {

    const auto x = i_p[0], y = i_p[1];
    const auto fxy = i_fun({x, y});
    const auto fp2dx = i_fun({x + 2 * i_d, y});
    const auto fm2dx = i_fun({x - 2 * i_d, y});
    const auto fp2dy = i_fun({x, y + 2 * i_d});
    const auto fm2dy = i_fun({x, y - 2 * i_d});
    const auto fpdxpdy = i_fun({x + i_d, y + i_d});
    const auto fmdxpdy = i_fun({x - i_d, y + i_d});
    const auto fpdxmdy = i_fun({x + i_d, y - i_d});
    const auto fmdxmdy = i_fun({x - i_d, y - i_d});

    sqmatrix result;
    result[0] = ((fp2dx - fxy) / (2 * i_d) - (fxy - fm2dx) / (2 * i_d)) / (2 * i_d); //d^2f/dx^2
    result[3] = ((fp2dy - fxy) / (2 * i_d) - (fxy - fm2dy) / (2 * i_d)) / (2 * i_d); //d^2f/dy^2
    result[1] = result[2] =
            ((fpdxpdy - fmdxpdy) / (2 * i_d) - (fpdxmdy - fmdxmdy) / (2 * i_d)) / (2 * i_d); //d^2f/dxdy = d^2f/dydx
    return result;
    }

  double
  Derivative(const std::function<double(const point&)>& i_fun, const point& i_p, const point& i_direction, double i_d)
    {
    const auto dp = Normalize(i_direction) * i_d;
    const auto f1 = i_fun(i_p);
    const auto f2 = i_fun(i_p + dp);
    return (f2 - f1) / Length(dp);
    }


  point OneDimensionalOptGolden(const std::function<double(const point&)>& i_fun, const point& i_p1, const point& i_p2,
                                double i_epsilon)
    {
    constexpr double GOLDEN = 1.6180339887498948482;
    const auto direction = Normalize(i_p2 - i_p1);

    auto df1 = Derivative(i_fun, i_p1, direction);
    auto df2 = Derivative(i_fun, i_p2, direction);

    if (df1 * df2 > 0)
      return (i_fun(i_p1) < i_fun(i_p1) ? i_p1 : i_p2);

    auto p1 = i_p1;
    auto p2 = i_p2;
    point min;
    do
      {
      min = DivideInRation(p1, p2, GOLDEN);
      df1 = Derivative(i_fun, p1, direction);
      df2 = Derivative(i_fun, min, direction);
      if (df1 * df2 > 0)
        p1 = min;
      else
        p2 = min;
      } while (Length(p1, p2) > i_epsilon);

    return min;
    }

  /*
   * TODO:i_bb is used only for one dimensional optimization and it is a bottle neck of the realisation
   * it will be great to find the limits for the optimization algorithmically
   */
  point PowellsOptimization(const std::function<double(const point&)>& i_fun, const point& i_p0, const bbox& i_bb,
                            double i_epsilon, bool i_print)
    {
    if (i_print)
      std::cout << "PowellsOptimization: " << std::endl;
    std::size_t k = 0, M = 20;
    sqmatrix A = {1., 0., 0., 1.};

    point prev, curr = i_p0;
    point prev_grad, curr_grad, curr_dirr;
    std::size_t success_counter = 0;

    curr_grad = Gradient(i_fun, curr);
    //10
    curr_dirr = A * curr_grad * (-1.);
    prev_grad = curr_grad;
    prev = curr;
    auto edge_for_opt = IntersectBBox(i_bb, curr, curr_dirr);
    curr = OneDimensionalOptGolden(i_fun, edge_for_opt.p1, edge_for_opt.p2);


    if (Length(curr - prev) < i_epsilon && std::abs(i_fun(curr) - i_fun(prev)) < i_epsilon)
      ++success_counter;

    if (i_print)
      {
      std::cout << "Start optimisation from point: " << prev;
      std::cout << " along  " << curr_dirr << " direction" << std::endl;
      std::cout << "gradient  " << prev_grad << " direction" << std::endl;
      std::cout << "Optimization from  " << edge_for_opt.p1 << " to " << edge_for_opt.p2 << std::endl;
      std::cout << std::endl;
      }

    do
      {
      //3
      curr_grad = Gradient(i_fun, curr);
      //4
      if (Length(curr_grad) < i_epsilon)
        break;
      //5
      if (k >= M)
        break;
      //6
      const auto d_grad = curr_grad - prev_grad;
      //7
      const auto d_p = curr - prev;
      //8
      const sqmatrix left = ColumnOnRow(d_p, d_p) / (d_p * d_grad);
      const sqmatrix right = ColumnOnRow((A * d_grad), (d_grad * A)) / ((d_grad * A) * d_grad);
      //9
      A += left - right;
      //10
      curr_dirr = A * curr_grad * (-1.);
      //11-12
      prev_grad = curr_grad;
      prev = curr;
      edge_for_opt = IntersectBBox(i_bb, curr, curr_dirr);
      curr = OneDimensionalOptGolden(i_fun, edge_for_opt.p1, edge_for_opt.p2);

      //13
      if (Length(curr - prev) < i_epsilon && std::abs(i_fun(curr) - i_fun(prev)) < i_epsilon)
        ++success_counter;
      else
        success_counter = 0;
      ++k;

      if (i_print)
        {
        std::cout << "Iteration: " << k << std::endl;
        std::cout << "\tOptimisation from point: " << prev;
        std::cout << " along  " << curr_dirr << " direction" << std::endl;
        std::cout << "\tgradient  " << prev_grad << " direction" << std::endl;
        std::cout << "\tOptimization from  " << edge_for_opt.p1 << " to " << edge_for_opt.p2 << std::endl;
        //std::cout << "A:\n" << A << std::endl;
        }
      } while (success_counter < 2);


    if (i_print)
      std::cout << "f(opt)=f(" << curr << ") = " << i_fun(curr) << std::endl;


    return curr;
    }

  point NewtonsOptimization(const std::function<double(const point&)>& i_fun, const point& i_p0, double i_epsilon,
                            bool i_print)
    {
    std::size_t k = 0, M = 20;
    point prev, curr = i_p0;
    do
      {
      if (k >= M)
        break;
      if (i_print)
        std::cout << "f(curr)=f(" << curr << ") = " << i_fun(curr) << std::endl;
      prev = curr;
      curr = curr - Invert(Hessian(i_fun, curr)) * Gradient(i_fun, curr);
      ++k;
      } while (Length(curr - prev) > i_epsilon || std::abs(i_fun(curr) - i_fun(prev)) > i_epsilon);

    if (i_print)
      std::cout << "f(opt)=f(" << curr << ") = " << i_fun(curr) << std::endl;


    return curr;
    }

  }