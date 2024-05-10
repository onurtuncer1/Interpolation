#if !defined(SIM_MATH_UTILS_H)

#include "sim_math_constants.h"

// NOTE(batuhan): 0 is considered positive
#define SIM_SIGN(x)					(((x) >= 0) ? 1 : -1)
#define SIM_SQUARE(a)				((a) * (a))
#define SIM_CUBE(a)					(SIM_SQUARE((a)) * (a))
#define SIM_FOURTH_POWER(a)			(SIM_CUBE((a)) * (a))
#define SIM_FIFTH_POWER(a)			(SIM_FOURTH_POWER((a)) * (a))
#define SIM_DEG_2_RAD(x)			((x) * SIM_DEG_2_RAD_FAC)
#define SIM_RAD_2_DEG(x)			((x) * SIM_RAD_2_DEG_FAC)

// NOTE(batuhan): The Art of Computer Programming by Donald Knuth, Volume II: Chapter IV
inline static
bool SIM_approximately_equal(double x, double y /*, double margin = 1 */)
{
	//return (fabs(x - y) <= (SIM_MAX(fabs(x), fabs(y)) * margin));
	return x == y;
}

inline static
bool SIM_essentially_equal(double x, double y/* , double margin = 1 */)
{
	//return (fabs(x - y) <= (SIM_MIN(fabs(x), fabs(y)) * margin));
	return x == y;
}

// x > y
inline static bool SIM_definitely_greater_than(double x, double y/* , double margin = 1 */) {
	//return (x - y > (SIM_MAX(fabs(x), fabs(y)) * margin));
	return x > y;
}

// x < y
inline static bool SIM_definitely_less_than(double x, double y/* , double margin = 1 */)
{
	//return (y - x > (SIM_MAX(fabs(x), fabs(y)) * margin));
	return x < y;
}

inline static double SIM_normalize_degree(double degree)
{
	double result = degree;

	/*while (result < 0.0)		result += 360;
	while (result > 360.0)		result -= 360;*/

	while (SIM_definitely_less_than(result, 0.0))		result += 360.0;
	while (SIM_definitely_greater_than(result, 360.0))	result -= 360.0;

	return result;
}

inline static double SIM_normalize_radian(double radian) {
	double result = radian;

	while (SIM_definitely_less_than(result, 0.0))			result += SIM_TWO_PI;
	while (SIM_definitely_greater_than(result, SIM_TWO_PI))	result -= SIM_TWO_PI;

	return result;
}

#define SIM_MATH_UTILS_H
#endif
