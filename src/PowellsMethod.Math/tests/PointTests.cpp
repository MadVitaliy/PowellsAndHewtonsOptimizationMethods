
#include <gtest/gtest.h>

#include "PowellsMethod.Math/point.h"

///////////////////////////////////////////////////////////////////////////////

TEST(DivideInRation, Test1) {
	constexpr double GOLDEN = 1.6180339887498948482;
	const point p1(1., -2.);
	const point p2(1., 2.);
	const auto golden = DivideInRation(p1, p2, GOLDEN);
	EXPECT_TRUE(IsEqual(golden, { 1., 0.47213595499957944 }));
}

TEST(DivideInRation, Test2) {
	const point p1(1., -2.);
	const point p2(1., 2.);
	const auto golden = DivideInRation(p1, p2, 1.0);
	EXPECT_TRUE(IsEqual(golden, { 1., 0.0 }));
}

///////////////////////////////////////////////////////////////////////////////

TEST(IntersectBBox, Test1) {
	const bbox bb{ {0., 0.}, {4., 4} };
	const point p(1., 2.);
	const point dir(3., 1.);
	const auto intersections = IntersectBBox(bb, p, dir);
	EXPECT_TRUE(IsEqual(intersections.p1, { 0.0, 1. + 2. / 3. }, 0.00001));
	EXPECT_TRUE(IsEqual(intersections.p2, { 4., 3. }));
}

TEST(IntersectBBox, Test2) {
	const bbox bb{ {0., 0.}, {4., 4} };
	const point p(2., 2.);
	const point dir(1., 1.);
	const auto intersections = IntersectBBox(bb, p, dir);
	EXPECT_TRUE(IsEqual(intersections.p1, { 0., 0. }) && IsEqual(intersections.p2, { 4., 4. }));
}

TEST(IntersectBBox, Test3) {
	const bbox bb{ {0., 0.}, {4., 4} };
	const point p(1., 2.);
	const point dir(1., -1.);
	const auto intersections = IntersectBBox(bb, p, dir);
	EXPECT_TRUE(IsEqual(intersections.p1, { 0., 3. }) && IsEqual(intersections.p2, { 3., 0. }));
}

///////////////////////////////////////////////////////////////////////////////

TEST(Invert, Test1_IdentityMatrixIsNotChanged) {
	const sqmatrix m{ 1., 0., 0., 1. };
	EXPECT_TRUE(IsEqual(m, Invert(m)));
}

TEST(Invert, Test2) {
	const auto inv_m = Invert({ 1., 2., 3., 4. });
	const sqmatrix golden_inv{ -2., 1., 1.5, -0.5 };
	EXPECT_TRUE(IsEqual(inv_m, golden_inv));
}
