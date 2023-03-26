#include <gtest/gtest.h>

#include "PowellsMethod.Math/point.h"
#include "PowellsMethod.Math/Optimization.h"

namespace {
	double QuadraticFunction(const point& i_p) {
		return i_p[0] * i_p[0] + i_p[0] * i_p[1] + 3 * i_p[1] * i_p[1] - 12 * i_p[0] - 15 * i_p[1] + 2;
	}

	double CubicFunction(const point& i_p) {
		return pow(i_p[0], 3.0) + pow(i_p[1], 3.0) - 3 * i_p[0] - 2 * i_p[1] + 2;
	}
}

///////////////////////////////////////////////////////////////////////////////

TEST(Gradient, Test1) {
	const auto grad = Optimization::Gradient(QuadraticFunction, { 0., 0. });
	const point golden_frad{ -12., -15. };
	EXPECT_TRUE(IsEqual(grad, golden_frad, 0.000001));
}

TEST(Gradient, Test2) {
	const auto grad = Optimization::Gradient(QuadraticFunction, { 0., 0. });
	const point golden_frad{ -12., -15. };
	EXPECT_TRUE(IsEqual(grad, golden_frad, 0.000001));
}

TEST(Gradient, Test4) {
	const auto grad = Optimization::Gradient(QuadraticFunction, { 2., 2. });
	const point golden_frad{ -6., -1 };
	EXPECT_TRUE(IsEqual(grad, golden_frad, 0.000001));
}

TEST(Gradient, Test5) {
	const auto grad = Optimization::Gradient(CubicFunction, { 0., 0. });
	const point golden_frad{ -3., -2. };
	EXPECT_TRUE(IsEqual(grad, golden_frad, 0.000001));
}

TEST(Gradient, Test6) {
	const auto grad = Optimization::Gradient(CubicFunction, { 0.9999999918566028, 0.8164965848727971 }, 0.0001);
	const point golden_frad{ 0., 0. };
	EXPECT_TRUE(IsEqual(grad, golden_frad, 0.000001));
}

///////////////////////////////////////////////////////////////////////////////

TEST(Hessian, Test1) {
	const auto hessian = Optimization::Hessian(QuadraticFunction, { 0., 0. });
	const sqmatrix golden_hessian{ 2., 1., 1., 6. };
	EXPECT_TRUE(IsEqual(hessian, golden_hessian, 0.0000001));
}

TEST(Hessian, Test2) {
	const auto hessian = Optimization::Hessian(QuadraticFunction, { -1., -1. });
	const sqmatrix golden_hessian{ 2., 1., 1., 6. };
	EXPECT_TRUE(IsEqual(hessian, golden_hessian, 0.0000001));
}

TEST(Hessian, Test3) {
	const auto hessian = Optimization::Hessian(CubicFunction, { -1., -1. });
	const sqmatrix golden_hessian{ -6., 0., 0., -6. };
	EXPECT_TRUE(IsEqual(hessian, golden_hessian, 0.0000001));
}

TEST(Hessian, Test4) {
	const auto hessian = Optimization::Hessian(CubicFunction, { 0.9999999918566028, 0.8164965848727971 });
	const sqmatrix golden_hessian{ 5.999999951139617, 0., 0., 4.898979509236782 };
	EXPECT_TRUE(IsEqual(hessian, golden_hessian, 0.0000001));
}

///////////////////////////////////////////////////////////////////////////////

TEST(NewtonsMethod, Test1) {
	constexpr auto precision = 0.0001;
	const auto opt = Optimization::NewtonsOptimization(CubicFunction, { 0.2, 0.2 }, precision, false);
	const point golden_opt{ 0.9999999918566028, 0.8164965848727971 };
	EXPECT_TRUE(IsEqual(opt, golden_opt, precision));
}

TEST(NewtonsMethod, Test2) {
	constexpr auto precision = 0.0001;
	const auto opt = Optimization::NewtonsOptimization(QuadraticFunction, { 0.0, 0.0 }, precision, false);
	const point golden_opt{ 5.181818207077867, 1.6363636597928 };
	EXPECT_TRUE(IsEqual(opt, golden_opt, precision));
}

///////////////////////////////////////////////////////////////////////////////

TEST(OneDimensionalOptGolden, Test1) {
	const point p1(1., 0.);
	const point p2(1., 2.);
	const auto min = Optimization::OneDimensionalOptGolden(CubicFunction, p1, p2, 0.0001);
	const auto ok = IsEqual(min, { 1., 0.8164965848727971 }, 0.0001);
}

TEST(OneDimensionalOptGolden, Test2) {
	const point p1(0., 0.);
	const point p2(2., 2.);
	const auto min = Optimization::OneDimensionalOptGolden(CubicFunction, p1, p2, 0.0001);
	const auto ok = IsEqual(min, { 0.91292, 0.91292 }, 0.0001);
}
