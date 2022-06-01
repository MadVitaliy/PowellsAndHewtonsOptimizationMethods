#include "point.h"
#include "Optimization.h"

#include <iostream>
#include <iomanip>
#include <cmath>


namespace Tests
  {
  void RunTests();
  }

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


int main()
  {
  Tests::RunTests();

//  const auto result = Optimization::PowellsOptimization(MinimizeMe, {0.5, 0.5}, {{0.4, 0.4},
//                                                                                 {2.,  2.}}, 0.0001, true);


  std::cout << "Newton1" << std::endl;
  Optimization::NewtonsOptimization(MinimizeMeByNewton, {0.5, 0.5}, 0.0001, true);

  return 0;
  }

namespace Tests
  {
  void TestDivideInRation()
    {
    constexpr double GOLDEN = 1.6180339887498948482;
    std::cout << "DivideInRationTest:" << std::endl;
    {
      const point p1(1., -2.);
      const point p2(1., 2.);
      const auto golden = DivideInRation(p1, p2, GOLDEN);
      const auto ok = IsEqual(golden, {1., 0.47213595499957944});
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
    }

    {
      const point p1(1., -2.);
      const point p2(1., 2.);
      const auto golden = DivideInRation(p1, p2, 1.0);
      const auto ok = IsEqual(golden, {1., 0.0});
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
    }
    }

  void TestOneDimensionalOptGolden()
    {
    std::cout << "TestOneDimensionalOptGolden:" << std::endl;
    {
      const point p1(1., 0.);
      const point p2(1., 2.);
      const auto min = Optimization::OneDimensionalOptGolden([](const point& i_p) -> double
                                                               { return MinimizeMe(i_p); }, p1, p2, 0.0001);
      const auto ok = IsEqual(min, {1., 0.8164965848727971}, 0.0001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      //std::cout << min << std::endl;
    }
    {
      const point p1(0., 0.);
      const point p2(2., 2.);
      const auto min = Optimization::OneDimensionalOptGolden([](const point& i_p) -> double
                                                               { return MinimizeMe(i_p); }, p1, p2, 0.0001);
      const auto ok = IsEqual(min, {0.91292, 0.91292}, 0.0001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      //std::cout << min << std::endl;
    }

    }

  void TestIntersectBBox()
    {
    std::cout << "TestIntersectBBox:" << std::endl;
    {
      const bbox bb{{0., 0.},
                    {4., 4}};
      const point p(1., 2.);
      const point dir(3., 1.);
      const auto intersections = IntersectBBox(bb, p, dir);
      const auto ok = IsEqual(intersections.p1, {0.0, 1. + 2. / 3.}, 0.00001) && IsEqual(intersections.p2, {4., 3.});
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
//      std::cout << "\t"<< intersections.p1 << std::endl;
//      std::cout << "\t"<< intersections.p2 << std::endl;
    }
    {
      const bbox bb{{0., 0.},
                    {4., 4}};
      const point p(2., 2.);
      const point dir(1., 1.);
      const auto intersections = IntersectBBox(bb, p, dir);
      const auto ok = IsEqual(intersections.p1, {0., 0.}) && IsEqual(intersections.p2, {4., 4.});
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
//      std::cout << "\t"<< intersections.p1 << std::endl;
//      std::cout << "\t"<< intersections.p2 << std::endl;
    }

    {
      const bbox bb{{0., 0.},
                    {4., 4}};
      const point p(1., 2.);
      const point dir(1., -1.);
      const auto intersections = IntersectBBox(bb, p, dir);
      const auto ok = IsEqual(intersections.p1, {0., 3.}) && IsEqual(intersections.p2, {3., 0.});
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
//      std::cout << "\t"<< intersections.p1 << std::endl;
//      std::cout << "\t"<< intersections.p2 << std::endl;
    }


    }

  void TestInvert()
    {
    std::cout << "TestInvert:" << std::endl;
    {
      const sqmatrix m{1., 0., 0., 1.};
      const auto inv_m = Invert(m);
      const auto ok = IsEqual(m, inv_m);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << inv_m << std::endl;
    }

    {
      const auto inv_m = Invert({1., 2., 3., 4.});
      const sqmatrix golden_inv{-2., 1., 1.5, -0.5};
      const auto ok = IsEqual(inv_m, golden_inv);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << inv_m << std::endl;
    }
    }

  void TestGradient()
    {
    std::cout << "TestGradient:" << std::endl;
    {
      const auto grad = Optimization::Gradient([](const point& p)
                                                 {
                                                 return p[0] * p[0] + p[0] * p[1] + 3 * p[1] * p[1] - 12 * p[0] -
                                                        15 * p[1] + 2;
                                                 }, {0., 0.});
      const point golden_frad{-12., -15.};
      const auto ok = IsEqual(grad, golden_frad, 0.000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected: {" << golden_frad << "}, real: {" << grad << "}" << std::endl;
    }
    {
      const auto grad = Optimization::Gradient([](const point& p)
                                                 {
                                                 return p[0] * p[0] + p[0] * p[1] + 3 * p[1] * p[1] - 12 * p[0] -
                                                        15 * p[1] + 2;
                                                 }, {2., 2.});
      const point golden_frad{-6., -1};
      const auto ok = IsEqual(grad, golden_frad, 0.000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected: {" << golden_frad << "}, real: {" << grad << "}" << std::endl;
    }
    {
      const auto grad = Optimization::Gradient([](const point& i_p)
                                                 {
                                                 return pow(i_p[0], 3.0) + pow(i_p[1], 3.0) - 3 * i_p[0] - 2 * i_p[1] +
                                                        2;
                                                 }, {0., 0.});
      const point golden_frad{-3., -2.};
      const auto ok = IsEqual(grad, golden_frad, 0.000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected: {" << golden_frad << "}, real: {" << grad << "}" << std::endl;
    }
    {
      const auto grad = Optimization::Gradient([](const point& i_p)
                                                 {
                                                 return pow(i_p[0], 3.0) + pow(i_p[1], 3.0) - 3 * i_p[0] - 2 * i_p[1] +
                                                        2;
                                                 }, {0.9999999918566028, 0.8164965848727971}, 0.0001);
      const point golden_frad{0., 0.};
      const auto ok = IsEqual(grad, golden_frad, 0.000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected: {" << golden_frad << "}, real: {" << grad << "}" << std::endl;
    }
    }

  void TestHessian()
    {
    std::cout << "TestHessian:" << std::endl;
    {
      const auto hessian = Optimization::Hessian([](const point& p)
                                                   {
                                                   return p[0] * p[0] + p[0] * p[1] + 3 * p[1] * p[1] - 12 * p[0] -
                                                          15 * p[1] + 2;
                                                   }, {0., 0.});
      const sqmatrix golden_hessian{2., 1., 1., 6.};
      const auto ok = IsEqual(hessian, golden_hessian, 0.0000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected:\n" << std::fixed << std::setprecision(10) << golden_hessian << "real:\n" << hessian
                  << std::endl;
    }
    {
      const auto hessian = Optimization::Hessian([](const point& p)
                                                   {
                                                   return p[0] * p[0] + p[0] * p[1] + 3 * p[1] * p[1] - 12 * p[0] -
                                                          15 * p[1] + 2;
                                                   }, {-1., -1.});
      const sqmatrix golden_hessian{2., 1., 1., 6.};
      const auto ok = IsEqual(hessian, golden_hessian, 0.0000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected:\n" << std::fixed << std::setprecision(10) << golden_hessian << "real:\n" << hessian
                  << std::endl;
    }
    {
      const auto hessian = Optimization::Hessian([](const point& i_p)
                                                   {
                                                   return pow(i_p[0], 3.0) + pow(i_p[1], 3.0) - 3 * i_p[0] -
                                                          2 * i_p[1] + 2;
                                                   }, {-1., -1.});
      const sqmatrix golden_hessian{-6., 0., 0., -6.};
      const auto ok = IsEqual(hessian, golden_hessian, 0.0000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected:\n" << std::fixed << std::setprecision(10) << golden_hessian << "real:\n" << hessian
                  << std::endl;
    }
    {
      const auto hessian = Optimization::Hessian([](const point& i_p)
                                                   {
                                                   return pow(i_p[0], 3.0) + pow(i_p[1], 3.0) - 3 * i_p[0] -
                                                          2 * i_p[1] + 2;
                                                   }, {0.9999999918566028, 0.8164965848727971});
      const sqmatrix golden_hessian{5.999999951139617, 0., 0., 4.898979509236782};
      const auto ok = IsEqual(hessian, golden_hessian, 0.0000001);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected:\n" << std::fixed << std::setprecision(10) << golden_hessian << "real:\n" << hessian
                  << std::endl;
    }
    }


  void TestNewtonsMethod()
    {
    std::cout << "TestNewtonsMethod:" << std::endl;
    {
      constexpr auto precision = 0.0001;
      const auto opt = Optimization::NewtonsOptimization([](const point& i_p)
                                                           {
                                                           return pow(i_p[0], 3.0) + pow(i_p[1], 3.0) - 3 * i_p[0] -
                                                                  2 * i_p[1] + 2;
                                                           }, {0.2, 0.2}, precision, false);
      const point golden_opt{0.9999999918566028, 0.8164965848727971};
      const auto ok = IsEqual(opt, golden_opt, precision);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected: {" << golden_opt << "}, found: {" << opt << "}" << std::endl;
    }
    {
      constexpr auto precision = 0.0001;
      const auto opt = Optimization::NewtonsOptimization([](const point& p)
                                                           {
                                                           return p[0] * p[0] + p[0] * p[1] + 3 * p[1] * p[1] -
                                                                  12 * p[0] - 15 * p[1] + 2;
                                                           }, {0.0, 0.0}, precision, false);
      const point golden_opt{5.181818207077867, 1.6363636597928};
      const auto ok = IsEqual(opt, golden_opt, precision);
      std::cout << (ok ? "\tok" : "\t ne ok") << std::endl;
      if (!ok)
        std::cout << "Expected: {" << golden_opt << "}, found: {" << opt << "}" << std::endl;
    }
    }

  void RunTests()
    {
    TestDivideInRation();
    TestOneDimensionalOptGolden();
    TestIntersectBBox();
    TestInvert();
    TestGradient();
    TestHessian();
    TestNewtonsMethod();
    }

  }