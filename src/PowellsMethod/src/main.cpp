#include "PowellsMethod.Math/point.h"
#include "PowellsMethod.Math/Optimization.h"

#include <iostream>
#include <iomanip>
#include <cmath>


namespace 
{
//function to minimize; return value of function in the i_point
double MinimizeMe(const point& i_point)
  {
  return pow(i_point[0], 3.0) + pow(i_point[1], 3.0) - 3 * i_point[0] - 2 * i_point[1] + 2;
  //return 4*i_point[0]*i_point[0] + 3*i_point[1]*i_point[1] - 4 * i_point[0]*i_point[1] - 2 * i_point[0];
  }

double MinimizeMeByNewton(const point& i_p)
  {
  return i_p[0] * i_p[0] + i_p[0] * i_p[1] + 3 * i_p[1] * i_p[1] - 12 * i_p[0] - 15 * i_p[1] + 2;
  }
}


int main()
  {
  const auto result = Optimization::PowellsOptimization(MinimizeMe, {0.41, 0.41}, {{0.4, 0.4},
                                                                                 {2.,  2.}}, 0.0001, true);
  std::cout << "Newton1" << std::endl;
  Optimization::NewtonsOptimization(MinimizeMeByNewton, {0.5, 0.5}, 0.0001, true);

  return 0;
  }
