#if !defined(SIM_FUNCTIONS_H)

#include "sim_utils.h"
//#include "sim_math_data_types.h"
#include <glm/glm.hpp>
#include <tuple>
#include <deque>
#include <vector>
#include <array>

enum Ref_Plane_Type
{
	PLANE_XY = 0,
	PLANE_YZ = 1,
	PLANE_XZ = 2
};

enum  Arc_Segment_Direction
{
	DIRECTION_CW = -1,
	DIRECTION_CCW = 1
};

struct Pose
{
	glm::dvec3 position; 	// workpiece coordinate, position
	glm::dvec3 orientation{0,0,1};	// workpiece coordinate, orientation
};

struct Interpolation_Input{
	Pose						start_pose; 			// in workpiece frame
	Pose						end_pose;				// in workpiece frame
	glm::dvec3					center_point;
	glm::dvec3					normal_vec;
	std::array<float,6>			start_joint;			// in machine frame
	std::array<float,6>			end_joint; 				// in machine frame
	double						raw_feedrate;
	double						feedrate;
	double						start_length;
	double						acceleration;
	double						deceleration;
	double						initial_feedrate = 0;
	std::array<double, 6> 		initial_feedrate_ruckig{0,0,0,0,0,0};
	double						final_feedrate = 0;
	double						helix;
	unsigned long long 			number_of_turns;
	Ref_Plane_Type				plane;
	Arc_Segment_Direction		direction;
	double						jerk_limit;
	double						time_step;
	double						one_over_time_step;
	bool 						requiresTransform = false;	// if true, calculate orientations with rodrigues formula
	bool 						isContouring = false;
};

struct Interpolation_Output
{
	std::deque<glm::dvec3>			tool_path;
	std::deque<glm::dvec3>			jerks; // haha
	std::deque<glm::dvec3>			velocities;
	std::deque<glm::dvec3>			accelerations;
	std::deque<glm::dvec3>			orientations;
	std::deque<glm::dvec3>			orientation_velocities;
	std::deque<glm::dvec3>			orientation_accelerations;
	std::deque<glm::dvec3>			orientation_jerks;

	// here for rapid traverse outputs
	std::deque<std::array<double, 6>> joint_positions;
	std::deque<std::array<double, 6>> joint_velocities;
	std::deque<std::array<double, 6>> joint_accelerations;
	std::deque<std::array<double, 6>> joint_jerks;

	// TODO(batuhan): This is here for debugging,
	// generate_spline function should return its velocity with end points
	double							NextInitialFeedRate;
	std::array<double,6>			NextInitialFeedRate_ruckig{0,0,0,0,0,0};	 	// necessary for ruckig driven segments. the size is equal to number of axes
	glm::dvec3						EndPoint;
};

struct Trapezoidal_Acceleration_Output {
	std::vector<double>		displacement_variation;
	std::vector<double>		velocity_variation;
	std::vector<double>		acceleration_variation;
	std::vector<double>		jerk_variation;
	std::vector<double>		delta;
};

struct Circle_Intersection_Output {
	bool					valid;
	double					x3[2];
	double					y3[2];
};

struct Spline_Generation_Output {
	bool					valid;
	double					spline_initial_speed;
	double					spline_final_speed;
};

unsigned long long  SIM_linear_interpolation(const Interpolation_Input& input, Interpolation_Output* result);
unsigned long long  SIM_linear_interpolation(const Interpolation_Input& input, Interpolation_Output* result);
void Rapid_Interpolation(const Interpolation_Input& input, Interpolation_Output* result);

unsigned long long  SIM_circular_interpolation(const Interpolation_Input& input, Interpolation_Output* result);
std::tuple<unsigned long long, double>  SIM_trajectory_with_spline(const Interpolation_Input& input, int input_type, Interpolation_Input& next_input, int next_input_type, Interpolation_Output* result, const double cornering_tolerance);

glm::dvec3 SIM_find_spline_endpoint_linesegment(const Interpolation_Input& input, bool end, const double cornering_tolerance);
glm::dvec3 SIM_find_spline_endpoint_arcsegment(const Interpolation_Input& input, bool ELEVEN, const double cornering_tolerance);
void SIM_convert_2D_to_3D(const glm::dvec3& start_point, const glm::dvec3& end_point, const glm::dvec3& corner_point, double* sp_x, double* sp_y, double* ep_x, double* ep_y, double* angle_before, double* angle_after, glm::dmat3 *R);
Circle_Intersection_Output SIM_circcirc(double x0, double y0, double r0, double x1, double y1, double r1);

#define SIM_FUNCTIONS_H
#endif // SIM_FUNCTIONS_H