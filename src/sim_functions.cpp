#include <vector>
#include <algorithm>
#include <limits>
#include "sim_functions.h"
#include "sim_math_utils.h"
#include "sim_utils.h"
#include <ruckig/ruckig.hpp>

#define ONE_OVER_SIX		0.166666667

static void SIM_trapezoidal_acceleration(const Interpolation_Input& input, const double travel_distance_of_tool, Trapezoidal_Acceleration_Output* result)
{

	auto jerk = SIM_MIN(input.acceleration * input.one_over_time_step, input.deceleration * input.one_over_time_step);
	jerk = SIM_MIN(jerk, input.jerk_limit);
	auto jerk_abs = fabs(jerk);

	//#pragma region ACCELERATION
	auto initial_feedrate_diff = (input.feedrate - input.initial_feedrate);
	auto sign_of_acc = SIM_SIGN(initial_feedrate_diff);
	auto acceleration = sign_of_acc * fabs(input.acceleration);
	auto acceleration_sign = SIM_SIGN(acceleration);
	auto J1 = acceleration_sign * jerk_abs;
	auto J3 = J1;
	auto one_over_J1 = 1.0 / J1;
	auto one_over_J3 = 1.0 / J3;
	auto T1 = 0.0;
	auto T2 = 0.0;
	auto T3 = 0.0;
	auto N1 = 0.0;
	auto N2 = 0.0;
	auto N3 = 0.0;
	auto acceleration_result = 0.0;
	auto acceleration_result_sign = 1;
	// TODO(batuhan): This is not going to work -- double comparison

	if (input.feedrate != input.initial_feedrate)
	{
		T1 = acceleration * one_over_J1;
		T3 = acceleration * one_over_J3;
		T2 = ((initial_feedrate_diff) / acceleration) - (T1);

		if (T2 < 0.0)
		{
			T2 = 0.0;
			acceleration = acceleration_sign * sqrt(J1 * initial_feedrate_diff);
			T1 = acceleration * one_over_J1;
			T3 = acceleration * one_over_J3;
		}
		acceleration_result = acceleration;
		acceleration_result_sign = SIM_SIGN(acceleration_result);

		N1 = round(fabs(T1) * input.one_over_time_step);
		N3 = round(fabs(T3) * input.one_over_time_step);
		N2 = round(T2 * input.one_over_time_step);
	}

	auto final_feedrate_diff = (input.feedrate - input.final_feedrate);
	auto sign_of_dec = SIM_SIGN(final_feedrate_diff);
	auto deceleration = sign_of_dec * fabs(input.deceleration);
	auto deceleration_sign = SIM_SIGN(deceleration);
	auto J5 = deceleration_sign * jerk_abs;
	auto J7 = J5;
	auto one_over_J5 = 1.0 / J5;
	auto one_over_J7 = 1.0 / J7;
	auto N5 = 0.0;
	auto N6 = 0.0;
	auto N7 = 0.0;
	auto T5 = 0.0;
	auto T6 = 0.0;
	auto T7 = 0.0;
	auto deceleration_result = 0.0;
	auto deceleration_result_sign = 1;
	// TODO(batuhan): This is not going to work -- double comparison
	if (input.feedrate != input.final_feedrate)
	{
		T5 = deceleration * one_over_J5;
		T7 = deceleration * one_over_J7;
		T6 = (final_feedrate_diff / deceleration) - (T5);

		if (T6 < 0)
		{
			T6 = 0.0;
			deceleration = deceleration_sign * sqrt(J5 * (final_feedrate_diff));
			T5 = deceleration * one_over_J5;
			T7 = deceleration * one_over_J7;
		}
		deceleration_result = deceleration;
		deceleration_result_sign = SIM_SIGN(deceleration_result); 				// TODO enes : There is "else". I did not change yet

		N5 = round(fabs(T5) * input.one_over_time_step);
		N6 = round(T6 * input.one_over_time_step);
		N7 = round(fabs(T7) * input.one_over_time_step);
	}

	auto a1 = 0.0;
	auto a2 = 0.0;
	auto a3 = 0.0;

	if (N1 != 0)
	{
		auto two_J1 = 2 * J1;
		a1 = SIM_ONE_HALF / acceleration_result;
		a2 = acceleration_result / two_J1;
		a3 = (acceleration_result * input.initial_feedrate / (two_J1)) - (SIM_SQUARE(input.initial_feedrate) / (2 * acceleration_result));
	}

	if (N5 != 0) {
		auto two_J5 = 2 * J5;
		a1 += SIM_ONE_HALF / deceleration_result;
		a2 += deceleration_result / two_J5;
		a3 += (deceleration_result * input.final_feedrate / (two_J5)) - (SIM_SQUARE(input.final_feedrate) / (2 * deceleration_result));
	}

	auto feedrate = input.feedrate;
	auto T4 = (1.0 / input.feedrate) * (travel_distance_of_tool - (a1 * SIM_SQUARE(input.feedrate) + a2 * input.feedrate + a3));
	auto N4 = 0;
	if (T4 <= 0)
	{
		T4 = 0;
		feedrate = (-a2 + sqrt(SIM_SQUARE(a2) - (4.0 * a1 * (a3 - travel_distance_of_tool)))) / (2.0 * a1);
		initial_feedrate_diff = feedrate - input.initial_feedrate;

		T2 = (initial_feedrate_diff / acceleration_result) - (acceleration_result * one_over_J1);
		if (T2 < 0)
		{
			T2 = 0;
			acceleration_result = acceleration_result_sign * sqrt(fabs(J1 * initial_feedrate_diff));
			T1 = acceleration_result * one_over_J1;
			T3 = T1;
		}
		N1 = round(fabs(T1) * input.one_over_time_step);
		N3 = N1;
		N2 = round(T2 * input.one_over_time_step);

		final_feedrate_diff = feedrate - input.final_feedrate;
		T6 = (final_feedrate_diff / deceleration_result) - (deceleration_result * one_over_J5);
		if (T6 < 0)
		{
			T6 = 0;
			deceleration_result = deceleration_result_sign * sqrt(fabs(J5 * final_feedrate_diff));
			T5 = deceleration_result * one_over_J5;
			T7 = T5;
		}
		N5 = round(fabs(T5) * input.one_over_time_step);
		N7 = N5;
		N6 = round(T6 * input.one_over_time_step);
	}
	N4 = round(T4 * input.one_over_time_step);

	auto f1s = input.initial_feedrate;
	auto l1s = input.start_length;
	for (int64_t k = 1; k <= N1; ++k)
	{
		const auto k_times_timestep = k * input.time_step;
		l1s = input.start_length + input.initial_feedrate * k_times_timestep + (J1 * SIM_CUBE(k_times_timestep)) * ONE_OVER_SIX;
		result->displacement_variation.push_back(l1s);

		f1s = input.initial_feedrate + (J1 * SIM_SQUARE(k_times_timestep)) * SIM_ONE_HALF;
		result->velocity_variation.push_back(f1s);

		result->acceleration_variation.push_back(J1 * k_times_timestep);
		result->jerk_variation.push_back(J1);
	}

	auto f2s = f1s;
	auto l2s = l1s;
	for (int64_t k = 1; k <= N2; ++k)
	{
		const auto k_times_timestep = k * input.time_step;
		l2s = l1s + f1s * k_times_timestep + acceleration_result * SIM_SQUARE(k_times_timestep) * SIM_ONE_HALF;
		result->displacement_variation.push_back(l2s);

		f2s = f1s + acceleration_result * k_times_timestep;
		result->velocity_variation.push_back(f2s);

		result->acceleration_variation.push_back(acceleration_result);
		result->jerk_variation.push_back(0);
	}

	auto f3s = f2s;
	auto l3s = l2s;
	for (int64_t k = 1; k <= N3; ++k)
	{
		const auto k_times_timestep = k * input.time_step;
		l3s = l2s + f2s * k_times_timestep + acceleration_result * SIM_SQUARE(k_times_timestep) * SIM_ONE_HALF - (J3 * SIM_CUBE(k_times_timestep)) * ONE_OVER_SIX;
		result->displacement_variation.push_back(l3s);

		f3s = f2s + acceleration_result * k_times_timestep - (J3 * SIM_SQUARE(k_times_timestep)) * SIM_ONE_HALF;
		result->velocity_variation.push_back(f3s);

		result->acceleration_variation.push_back(acceleration_result - J3 * k_times_timestep);
		result->jerk_variation.push_back(-J3);
	}

	auto f4s = f3s;
	auto l4s = l3s;
	for (int64_t k = 1; k <= N4; ++k)
	{
		const auto k_times_timestep = k * input.time_step;
		l4s = l3s + f3s * k_times_timestep;
		result->displacement_variation.push_back(l4s);

		f4s = f3s;
		result->velocity_variation.push_back(f4s);

		result->acceleration_variation.push_back(0);
		result->jerk_variation.push_back(0);
	}

	auto f5s = f4s;
	auto l5s = l4s;
	for (int64_t k = 1; k <= N5; ++k)
	{
		const auto k_times_timestep = k * input.time_step;
		l5s = l4s + f4s * k_times_timestep - (J5 * SIM_CUBE(k_times_timestep)) * ONE_OVER_SIX;
		result->displacement_variation.push_back(l5s);

		f5s = f4s - (J5 * SIM_SQUARE(k_times_timestep) * SIM_ONE_HALF);
		result->velocity_variation.push_back(f5s);

		result->acceleration_variation.push_back(-J5 * k_times_timestep);
		result->jerk_variation.push_back(-J5);
	}

	auto f6s = f5s;
	auto l6s = l5s;
	for (int64_t k = 1; k <= N6; ++k) {

		const auto k_times_timestep = k * input.time_step;
		l6s = l5s + f5s * k_times_timestep - deceleration_result * SIM_SQUARE(k_times_timestep) * SIM_ONE_HALF;
		result->displacement_variation.push_back(l6s);

		f6s = f5s - deceleration_result * k_times_timestep;
		result->velocity_variation.push_back(f6s);

		result->acceleration_variation.push_back(-deceleration_result);
		result->jerk_variation.push_back(0);
	}

	auto f7s = f6s;
	auto l7s = l6s;
	for (int64_t k = 1; k <= N7; ++k)
	{
		const auto k_times_timestep = k * input.time_step;
		l7s = l6s + f6s * k_times_timestep - (deceleration_result * (SIM_SQUARE(k_times_timestep)) * SIM_ONE_HALF) + (J7 * SIM_CUBE(k_times_timestep) * ONE_OVER_SIX);
		result->displacement_variation.push_back(l7s);

		f7s = f6s - deceleration_result * k_times_timestep + (J7 * SIM_SQUARE(k_times_timestep)) * SIM_ONE_HALF; 	// enes : deceleration -> deceleration_result
		result->velocity_variation.push_back(f7s);

		result->acceleration_variation.push_back(-deceleration_result + J7 * k_times_timestep);
		result->jerk_variation.push_back(J7);
	}

	auto kkkk = 1.0;
	if (result->displacement_variation.back() != travel_distance_of_tool) 												// Throws exception here if the start and end tool tip positions are equal !!!!
	{
		kkkk = (result->displacement_variation.back() - input.start_length) / travel_distance_of_tool;
	}

	for (auto& it : result->displacement_variation)
	{
		double delta = it - input.start_length;
		it = (delta / kkkk) + input.start_length;
	}

	double start_length = input.start_length;
	for (const auto& it : result->displacement_variation)
	{
		double delta = it - start_length;
		start_length = it;
		result->delta.push_back(delta);
	}

	kkkk = 1.0;
	// TODO(batuhan): This is not going to work -- double comparison / what is 
	if (result->velocity_variation.back() != input.final_feedrate)
	{
		// TODO(batuhan): division by 0?
		kkkk = (result->velocity_variation.back() / input.final_feedrate);
	}

	if (input.final_feedrate == 0)
	{
		kkkk = 1;
		result->velocity_variation.back() = 0;
	}

	for (int i = 0; i < int(result->velocity_variation.size()); ++i)
	{
		result->velocity_variation[i] = result->velocity_variation[i] / kkkk;
	}
}

void SIM_convert_2D_to_3D(const glm::dvec3& start_point, const glm::dvec3& end_point, const glm::dvec3& corner_point, double* sp_x, double* sp_y, double* ep_x, double* ep_y, double* angle_before, double* angle_after, glm::dmat3* R) {
	
	glm::dvec3 V1 = start_point - corner_point;
	glm::dvec3 V2 = end_point - corner_point;
	// isn't this just length(V1)?
	double L1 = glm::length(V1);
	double L2 = glm::length(V2);
	glm::dvec3 UV1 = V1 / L1;
	glm::dvec3 UV2 = V2 / L2;

	glm::dvec3 V2_cross_V1 = glm::cross(V2, V1);
	// NOTE(batuhan): divisor is V1_cross_V2 in the matlab
	glm::dvec3 normal_vec = V2_cross_V1 / (length(V2_cross_V1) + std::numeric_limits<double>::epsilon());

	if (V2_cross_V1.x == 0 &&	V2_cross_V1.y == 0 &&		V2_cross_V1.z == 0)	{
		
		normal_vec.z = 1;
	}

	UV2 = -1.0* glm::cross(UV1, normal_vec);

	*R = { UV1, UV2, normal_vec };
	glm::dmat3 R_transpose = transpose(*R);

	glm::dvec3 planar_A = R_transpose * V1;
	glm::dvec3 planar_B = R_transpose * V2;

	double theta_start = atan2(planar_A[1], planar_A[0]);
	double theta_end = atan2(planar_B[1], planar_B[0]);
	theta_start = SIM_normalize_radian(theta_start);
	theta_end = SIM_normalize_radian(theta_end);

	double xs = L1 * cos(theta_start);
	double ys = L1 * sin(theta_start);
	double xe = L2 * cos(theta_end);
	double ye = L2 * sin(theta_end);

	// TODO(batuhan): group these
	*angle_before = atan2(-ys, -xs);
	*angle_after = atan2(ye, xe);
	*sp_x = xs;
	*sp_y = ys;
	*ep_x = xe;
	*ep_y = ye;
}

// Intersection of two circles: http://paulbourke.net/geometry/circlesphere/
Circle_Intersection_Output SIM_circcirc(double x0, double y0, double r0, double x1, double y1, double r1) {

	Circle_Intersection_Output result;
	result.valid = false;

	double d2 = (SIM_SQUARE(x0 - x1) + SIM_SQUARE(y0 - y1));
	// d = a + b
	double d = sqrt(d2);

	if ((d > r0 + r1) || d < fabs(r0 - r1)) return result;

	double a = (SIM_SQUARE(r0) - SIM_SQUARE(r1) + (d2)) / (2.0 * d);
	double h = sqrt(SIM_SQUARE(r0) - SIM_SQUARE(a));

	double x2 = x0 + a * (x1 - x0) / d;
	double y2 = y0 + a * (y1 - y0) / d;

	result.valid = true;
	result.x3[0] = x2 - h * (y1 - y0) / d;
	result.y3[0] = y2 + h * (x1 - x0) / d;
	result.x3[1] = x2 + h * (y1 - y0) / d;
	result.y3[1] = y2 - h * (x1 - x0) / d;

	return result;
}

struct orientation_profile {
    // Member variables
    std::deque<glm::dvec3> orientations;
    std::deque<glm::dvec3> orientation_velocities;
    std::deque<glm::dvec3> orientation_accelerations;
    std::deque<glm::dvec3> orientation_jerks;
    std::deque<double> theta;
    std::deque<double> theta_velocity;
    std::deque<double> theta_acceleration;
    glm::dvec3 normal_vector_of_rotation{};

    // Constructor
    orientation_profile(const std::deque<glm::dvec3>& ori,
                        const std::deque<glm::dvec3>& vel,
                        const std::deque<glm::dvec3>& acc,
                        const std::deque<glm::dvec3>& jer,
                        const std::deque<double>& th,
                        const std::deque<double>& th_vel,
                        const std::deque<double>& th_acc,
                        const glm::vec3& normal = glm::vec3(0.0f))
        : orientations(ori),
          orientation_velocities(vel),
          orientation_accelerations(acc),
          orientation_jerks(jer),
          theta(th),
          theta_velocity(th_vel),
          theta_acceleration(th_acc),
          normal_vector_of_rotation(normal) {/*empty*/}
	
}; // struct orientation_profile


double get_angle_between_vectors(const Interpolation_Input& input) // calculates angle between orientations (rad)
{
	glm::vec3 normalized_start_orient = glm::normalize(input.start_pose.orientation);
	glm::vec3 normalized_end_orient = glm::normalize(input.end_pose.orientation);
	return glm::acos(glm::dot(normalized_start_orient, normalized_end_orient));
}

orientation_profile get_mid_orientation_with_rodrigues(const Interpolation_Input& input,const Trapezoidal_Acceleration_Output& mid_res)
{
	
	if (input.start_pose.orientation == input.end_pose.orientation) 	// TODO check this double comparison
	{
		std::deque<glm::dvec3> orientations(mid_res.displacement_variation.size(), input.start_pose.orientation); 		// If there is no change in orientation, equalize it
		
		std::deque<glm::dvec3> orientation_vels(mid_res.displacement_variation.size(), {{0}, {0}, {0}});
		std::deque<glm::dvec3> orientation_accs(orientation_vels);
		std::deque<glm::dvec3> orientation_jerks(orientation_vels);

		std::deque<double> theta(mid_res.displacement_variation.size(), 0.0);
		std::deque<double> theta_vel(theta);
		std::deque<double> theta_acc(theta);

		return orientation_profile{orientations, orientation_vels, orientation_accs, orientation_jerks, theta, theta_vel, theta_acc};
	}
	else
	{
		glm::dvec3 normal_vector = glm::cross(input.start_pose.orientation, input.end_pose.orientation);
		// glm::normalize(normal_vector);
		normal_vector = normal_vector / glm::length(normal_vector);

		double angle_between_orientations = get_angle_between_vectors(input); // get angle between end and start orientations

		double start_length_previous = input.start_length;

		std::deque<double> theta(1, 0.0);
		std::deque<double> theta_vel(1, 0.0);
		std::deque<double> theta_acc(1, 0.0);
		std::deque<double> theta_jerk(1, 0.0);

		double theta_0 = 0.0;
		
		for (int i = 0; i < int(mid_res.displacement_variation.size()); i++) // Fill theta vector
		{
			auto theta_temp = (mid_res.displacement_variation[i] - start_length_previous) / (mid_res.displacement_variation.back() - input.start_length) * angle_between_orientations + theta_0;
			theta_0 = theta_temp;
			theta.emplace_back(theta_temp);
			start_length_previous = mid_res.displacement_variation[i];
		}

		for (int i = 0; i < int(mid_res.displacement_variation.size()); i++) // Fill theta velocity vector
		{
			auto theta_vel_temp = (mid_res.velocity_variation[i]) / (mid_res.displacement_variation.back() - input.start_length) * angle_between_orientations;
			theta_vel.emplace_back(theta_vel_temp);
		}

		for (int i = 0; i < int(mid_res.displacement_variation.size()); i++) // Fill theta acc vector
		{
			auto theta_acc_temp = (mid_res.acceleration_variation[i]) / (mid_res.displacement_variation.back() - input.start_length) * angle_between_orientations;
			theta_acc.emplace_back(theta_acc_temp);
		}

		for (int i = 0; i < int(mid_res.displacement_variation.size()); i++) // Fill theta jerk vector
		{
			auto theta_jerk_temp = (mid_res.jerk_variation[i]) / (mid_res.displacement_variation.back() - input.start_length) * angle_between_orientations;
			theta_jerk.emplace_back(theta_jerk_temp);
		}

		std::deque<glm::dvec3> orientations;
		std::deque<glm::dvec3> orientation_vels;
		std::deque<glm::dvec3> orientation_accs;
		std::deque<glm::dvec3> orientation_jerks;

		for (int i = 0; i < int(theta.size()); i++) // Fill orientation
		{
			auto orientation_temp = input.start_pose.orientation * cos(theta[i]) 
												+ glm::cross(normal_vector, input.start_pose.orientation) * sin(theta[i]) 
												+ normal_vector * (glm::dot(normal_vector, input.start_pose.orientation)) * (1 - cos(theta[i]));
			
			orientations.emplace_back(orientation_temp);
		}

		for (int i = 0; i < int(theta.size()); i++) // Fill orientation vel
		{
			auto orientation_vel_temp = -input.start_pose.orientation * sin(theta[i]) * theta_vel[i]
                + glm::cross(normal_vector, input.start_pose.orientation) * cos(theta[i]) * theta_vel[i]
                + normal_vector * (glm::dot(normal_vector, input.start_pose.orientation)) * (sin(theta[i]) * theta_vel[i]);
			
			orientation_vels.emplace_back(orientation_vel_temp);
		}

		for (int i = 0; i < int(theta.size()); i++) // Fill orientation acc
		{
			auto orientation_acc_temp = - input.start_pose.orientation * cos(theta[i]) * SIM_SQUARE(theta_vel[i])
											- input.start_pose.orientation * sin(theta[i]) * theta_acc[i]
											- cross(normal_vector, input.start_pose.orientation) * (sin(theta[i]) * SIM_SQUARE(theta_vel[i]) - cos(theta[i])*theta_acc[i])
											+ normal_vector * (glm::dot(normal_vector,input.start_pose.orientation )) * (cos(theta[i]) * SIM_SQUARE(theta_vel[i]) + sin(theta[i]) * theta_acc[i]);
			
			orientation_accs.emplace_back(orientation_acc_temp);
		}

		for (int i = 0; i < int(theta.size()); i++) // Fill orientation jerk
		{
			auto orientation_jerk_temp = input.start_pose.orientation * sin(theta[i]) * SIM_CUBE(theta_vel[i]) - input.start_pose.orientation * cos(theta[i]) * 2.0 * theta_acc[i]
											- input.start_pose.orientation * cos(theta[i]) * theta_vel[i] * theta_acc[i] - input.start_pose.orientation * sin(theta[i]) * 2.0 * theta_jerk[i]
											+ glm::cross(normal_vector, input.start_pose.orientation) * (-cos(theta[i]) * SIM_CUBE(theta_vel[i]) - sin(theta[i]) * 2.0 * theta_acc[i])
											+ glm::cross(normal_vector,input.start_pose.orientation) * (-sin(theta[i]) * theta_vel[i] * theta_acc[i] + cos(theta[i]) * theta_jerk[i])
											+ normal_vector * (glm::dot(normal_vector, input.start_pose.orientation)) * (-sin(theta[i]) * SIM_CUBE(theta_vel[i]) + cos(theta[i]) * 2.0 * theta_acc[i])
											+ normal_vector * (glm::dot(normal_vector, input.start_pose.orientation)) * (cos(theta[i]) * theta_vel[i] * theta_acc[i] + sin(theta[i]) * theta_jerk[i]);
			
			orientation_jerks.emplace_back(orientation_jerk_temp);
		}
		
		return orientation_profile{orientations, orientation_vels, orientation_accs, orientation_jerks, theta, theta_vel, theta_acc};
	}

}

unsigned long long  SIM_linear_interpolation(const Interpolation_Input& input, Interpolation_Output* result) {
	
	const auto& sp = input.start_pose.position;
	const auto& ep = input.end_pose.position;

	glm::dvec3 delta = ep - sp;
	double travel_distance_of_tool = glm::length(delta); 			// if this is zero throws exception

	double epsilon = 0.003; 			// TODO take this to Parameter Manager
	if ((travel_distance_of_tool < epsilon && !input.requiresTransform) || input.feedrate == 0.0) 			// if travel distance is zero or feed rate is zero, generate no trajectory
	{
		result->tool_path.resize(1);
		result->tool_path[0] = glm::dvec3{input.end_pose.position[0], input.end_pose.position[1], input.end_pose.position[2]};
		return std::size_t(1);
	}
	else if ((travel_distance_of_tool < epsilon) 
				&& input.end_pose.orientation != input.start_pose.orientation 
				&& input.requiresTransform) 								// if tool tip is stationary with respoect to workpiece and tool just orientates, interpolate tool orientation using circular interp
	{
		Interpolation_Input input_arc;
		Interpolation_Output result_arc;
		input_arc.acceleration = input.acceleration;
		input_arc.deceleration = input.deceleration;
		input_arc.jerk_limit = input.jerk_limit;
		input_arc.time_step = input.time_step;
		input_arc.one_over_time_step = 1.0 / input_arc.time_step;
		input_arc.normal_vec = glm::cross(input.start_pose.orientation, input.end_pose.orientation);
		input_arc.normal_vec = input_arc.normal_vec / glm::length(input_arc.normal_vec);
		input_arc.center_point = glm::dvec3{0, 0, 0};
		input_arc.end_pose.position = input.end_pose.orientation;
		input_arc.start_pose.position = input.start_pose.orientation;
		input_arc.feedrate = input.feedrate / 2;
		SIM_circular_interpolation(input_arc, &result_arc);
		result->tool_path.clear();		
		result->tool_path.resize(result_arc.tool_path.size());
		for (size_t i = 0; i < result_arc.tool_path.size(); i++)
		{
			result->tool_path[i] = input.end_pose.position;
		}
		result->orientations.clear();
		result->orientations.resize(result_arc.tool_path.size());
		result->orientations = result_arc.tool_path;
		return result->orientations.size(); 
	}
	else if ((travel_distance_of_tool > epsilon)
				&& input.end_pose.orientation == input.start_pose.orientation
				&& input.requiresTransform) 									// if tool tip position changes but tool orientation does not change
	{
		Trapezoidal_Acceleration_Output mid_res = {};
		SIM_trapezoidal_acceleration(input, travel_distance_of_tool, &mid_res);

		for (int i = 0; i < int(mid_res.displacement_variation.size()); ++i)	{

			result->tool_path.emplace_back(sp + delta * (mid_res.displacement_variation[i] - input.start_length) / travel_distance_of_tool);
			result->jerks.emplace_back(delta * (mid_res.jerk_variation[i] / travel_distance_of_tool));
			result->velocities.emplace_back(delta * (mid_res.velocity_variation[i] / travel_distance_of_tool));
			result->accelerations.emplace_back(delta * (mid_res.acceleration_variation[i] / travel_distance_of_tool));
			result->NextInitialFeedRate = mid_res.velocity_variation.back();
		}
		result->orientations.clear();
		result->orientations.resize(result->tool_path.size());

		for (size_t i = 0; i < result->orientations.size(); i++)
		{
			result->orientations[i] = input.end_pose.orientation;
		}


		return mid_res.displacement_variation.size();

	}
	else
	{
		Trapezoidal_Acceleration_Output mid_res = {};
		SIM_trapezoidal_acceleration(input, travel_distance_of_tool, &mid_res);

		for (int i = 0; i < int(mid_res.displacement_variation.size()); ++i)	{

			result->tool_path.emplace_back(sp + delta * (mid_res.displacement_variation[i] - input.start_length) / travel_distance_of_tool);
			result->jerks.emplace_back(delta * (mid_res.jerk_variation[i] / travel_distance_of_tool));
			result->velocities.emplace_back(delta * (mid_res.velocity_variation[i] / travel_distance_of_tool));
			result->accelerations.emplace_back(delta * (mid_res.acceleration_variation[i] / travel_distance_of_tool));
			result->NextInitialFeedRate = mid_res.velocity_variation.back();
		}

		if (input.requiresTransform)
		{
			orientation_profile tempOrientationProfile = get_mid_orientation_with_rodrigues(input, mid_res);
			result->orientations 				= tempOrientationProfile.orientations;
			result->orientation_velocities 		= tempOrientationProfile.orientation_velocities;
			result->orientation_accelerations 	= tempOrientationProfile.orientation_accelerations;
			result->orientation_jerks 			= tempOrientationProfile.orientation_jerks;
		}
		return mid_res.displacement_variation.size();
	}

	
}
void Rapid_Interpolation(const Interpolation_Input& input, Interpolation_Output* result)
{
	using namespace ruckig;

    Ruckig<6> otg {input.time_step};
    InputParameter<6> input_ruckig;
    OutputParameter<6> output_ruckig;

    input_ruckig.current_position = {	input.start_joint[0], 
										input.start_joint[1],
										input.start_joint[2],
										input.start_joint[3],
										input.start_joint[4],
										input.start_joint[5]};

    input_ruckig.current_velocity = input.initial_feedrate_ruckig;
    input_ruckig.current_acceleration = {0, 0, 0, 0, 0, 0}; 		// TODO also equalize this 

    input_ruckig.target_position = {	input.end_joint[0], 
										input.end_joint[1],
										input.end_joint[2],
										input.end_joint[3],
										input.end_joint[4],
										input.end_joint[5]};

	if (input_ruckig.current_position == input_ruckig.target_position) 		// if current pos and final pos are equal return single tick
	{
		result->joint_positions.push_back(std::array<double, 6>{input.start_joint[0],input.start_joint[1],input.start_joint[2],input.start_joint[3],input.start_joint[4],input.start_joint[5]});
		result->joint_velocities.push_back(std::array<double, 6>{0,0,0,0,0,0});
		result->joint_accelerations.push_back(std::array<double, 6>{0,0,0,0,0,0});
		result->joint_jerks.push_back(std::array<double, 6>{0,0,0,0,0,0});
		return;
	}

	if (input.isContouring)
	{
    	input_ruckig.target_velocity = {input.final_feedrate, input.final_feedrate, input.final_feedrate, input.final_feedrate, input.final_feedrate, input.final_feedrate};
	}
	else
	{
    	input_ruckig.target_velocity = {0,0,0,0,0,0};
	}

    input_ruckig.max_velocity = {input.feedrate, input.feedrate, input.feedrate, input.feedrate, input.feedrate, input.feedrate};
    input_ruckig.max_acceleration = {input.acceleration,input.acceleration,input.acceleration,input.acceleration,input.acceleration,input.acceleration};
    input_ruckig.max_jerk = {input.jerk_limit,input.jerk_limit,input.jerk_limit,input.jerk_limit,input.jerk_limit,input.jerk_limit};

    while (otg.update(input_ruckig, output_ruckig) == Result::Working) 
    {
        auto& p = output_ruckig.new_position;
        auto& v = output_ruckig.new_velocity;
        auto& a = output_ruckig.new_acceleration;
        auto& j = output_ruckig.new_jerk;

		result->joint_positions.push_back(std::array<double, 6>{p[0],p[1],p[2],p[3],p[4],p[5]});
		result->joint_velocities.push_back(std::array<double, 6>{v[0],v[1],v[2],v[3],v[4],v[5]});
		result->NextInitialFeedRate_ruckig = std::array<double, 6>{{v[0],v[1],v[2],v[3],v[4],v[5]}};
		result->joint_accelerations.push_back(std::array<double, 6>{a[0],a[1],a[2],a[3],a[4],a[5]});
		result->joint_jerks.push_back(std::array<double, 6>{j[0],j[1],j[2],j[3],j[4],j[5]});

        output_ruckig.pass_to_input(input_ruckig);
    }

}

unsigned long long  SIM_circular_interpolation(const Interpolation_Input& input, Interpolation_Output* result) {

	glm::dvec3 sp = input.start_pose.position;
	glm::dvec3 ep = input.end_pose.position;
	glm::dvec3 cp = input.center_point;
	glm::dvec3 nv = input.normal_vec;	

	glm::dvec3 Oprime_sp = sp - cp;
	glm::dvec3 Oprime_ep = ep - cp;

	double R_sp = glm::length(Oprime_sp);
	double R_ep = glm::length(Oprime_ep);
	double k = R_ep / R_sp;

	glm::dvec3 UA = Oprime_sp / R_sp;
	glm::dvec3 UB = Oprime_ep / R_ep;
	Oprime_ep /= k;
	
	UA = glm::normalize(Oprime_sp);

	UB = -1.0 * glm::cross(UA, nv);
	
	int8_t dir = 1;

	if (glm::length(glm::cross(UA, UB) - nv) > 0.001)
	{
		UB *= - 1;
		dir = -1;
	}

	glm::dmat3 R = { UA, UB, nv };
	glm::dmat3 R_transpose = glm::transpose(R);
	glm::vec3 planar_A = R_transpose * Oprime_sp;
	glm::vec3 planar_B = R_transpose * Oprime_ep;
	double theta_start = atan2(planar_A[1], planar_A[0]);
	double theta_end = atan2(planar_B[1], planar_B[0]);

	theta_start = 180.0 / SIM_PI * theta_start;
	theta_end = 180.0 / SIM_PI * theta_end;

	if (-0.01 < theta_start && theta_start < 0.01)
	{
		theta_start = 0.0;
	}

	if (-0.01 < theta_end && theta_end < 0.01)
	{
		theta_end = 0.0;
	}

	if (dir > 0)
	{
		if (theta_start < 0)
		{
			theta_start = 360 + theta_start;
		}
		if (theta_end < 0)
		{
			theta_end = 360 + theta_end;
		}
		if (theta_end < theta_start)
		{
			theta_end = 360 + theta_end;
		}
		
	}

	if (dir < 0)
	{
		double temp = theta_start;
		theta_start = theta_end;
		theta_end = temp;
		if (theta_start < 0)
		{
			theta_start = 360 + theta_start;
		}
		if (theta_end < 0)
		{
			theta_end = 360 + theta_end;
		}
		if (theta_end < theta_start)
		{
			theta_end = 360 + theta_end;
		}

	}

	if (theta_end == 0.0)
	{
		theta_end = 360;
	}

	if (theta_start == 360.0)
	{
		theta_start = 0.0;
	}

	theta_start = theta_start * SIM_PI / 180.0;
	theta_end = theta_end * SIM_PI / 180.0;

	double xs = R_sp * cos(theta_start);
	double ys = R_sp * sin(theta_start);

	double travel_distance_of_tool = R_sp * (theta_end - theta_start);

	Trapezoidal_Acceleration_Output mid_res = {};
	SIM_trapezoidal_acceleration(input, travel_distance_of_tool, &mid_res);

	for (auto& it : mid_res.delta){
		it /= R_sp;
	}

	std::vector<glm::vec3> new_points;
	new_points.push_back({ xs, ys, 0 });
	for (int i = 0; i < int(mid_res.delta.size()); ++i)
	{
		new_points.push_back(
			{
				(new_points[i].x * cos(mid_res.delta[i])) - (sin(mid_res.delta[i]) * new_points[i].y),
				(new_points[i].y * cos(mid_res.delta[i])) + (sin(mid_res.delta[i]) * new_points[i].x),
				0
			});
	}
	if (dir < 1)
	{
		std::reverse(new_points.begin(), new_points.end());		
	}


	for (int i = 1; i < int(new_points.size()); ++i)
	{
		result->tool_path.push_back(cp + R * new_points[i]);
	}

	double angle = theta_start;

	dir = 1;
	for (int i = 0; i < int(mid_res.displacement_variation.size()); ++i)	{

		double pos_dv = angle + ((mid_res.displacement_variation[i] - input.start_length) * double(dir) / R_sp);
		double pos_dv_cos = R_sp * cos(pos_dv);
		double pos_dv_sin = R_sp * sin(pos_dv);

		double vel_dv = mid_res.velocity_variation[i] / R_sp * double(dir);
		double acc_dv = mid_res.acceleration_variation[i] / R_sp * double(dir);
		double jerk_dv = mid_res.jerk_variation[i] / R_sp * double(dir);

		auto velocity = R * glm::dvec3(-pos_dv_sin * vel_dv, +pos_dv_cos * vel_dv, 0.0 );
		result->velocities.push_back(velocity);

		auto acceleration = R * glm::dvec3(-result->velocities.back().y * vel_dv - pos_dv_sin * acc_dv,
											+result->velocities.back().x * vel_dv + pos_dv_cos * acc_dv,
											0.0 );
		result->accelerations.push_back(acceleration);

		auto jerk = R * glm::dvec3(-result->accelerations.back().y * vel_dv - 2.0 * result->velocities.back().y * acc_dv - pos_dv_sin * jerk_dv,
									+result->accelerations.back().x * vel_dv + 2.0 * result->velocities.back().x * acc_dv + pos_dv_cos * jerk_dv,
									0.0);
		result->jerks.push_back(jerk);

		result->NextInitialFeedRate = mid_res.velocity_variation.back();
	}

	orientation_profile tempOrientationProfile = get_mid_orientation_with_rodrigues(input, mid_res);
	result->orientations 				= tempOrientationProfile.orientations;
	result->orientation_velocities 		= tempOrientationProfile.orientation_velocities;
	result->orientation_accelerations 	= tempOrientationProfile.orientation_accelerations;
	result->orientation_jerks 			= tempOrientationProfile.orientation_jerks;

	return mid_res.displacement_variation.size();
}

//void SIM_generate_spline(const v3& start_point, const v3& corner_point, const v3& end_point)
//{
//	constexpr double spline_travel_time_increment = 0.001;
//	bool spline_fit = false;
//	Interpolation_Output mid_res;
//
//	double spline_sp_x, spline_sp_y;
//	double spline_ep_x, spline_ep_y;
//	double angle_before, angle_after;
//	mat3 R;
//	SIM_convert_2D_to_3D(start_point, end_point, corner_point,
//						 &spline_sp_x, &spline_sp_y, &spline_ep_x, &spline_ep_y,
//						 &angle_before, &angle_after, &R);
//
//	// TODO(batuhan): -v- v0, v1
//	double v0 = 0.0;
//	double v1 = 0.0;
//	double WTF_x = 0;
//	double WTF_y = 0;
//
//	for (double k = 1.0; !spline_fit; k = k * 0.999)
//	{
//		v0 = k * v0;
//		v1 = k * v1;
//		double check_two = 0.01;
//		std::vector<double> check_one(11, 1.0);
//
//		int iteration_number = 0;
//		double min_ttt = std::numeric_limits<double>::infinity();
//		double current_spline_travel_time = spline_travel_time_increment;
//		// spline fit
//		while ((*std::max_element(check_one.begin(), check_one.end()) > 0.0) ||
//			   check_two > 0.01 && min_ttt < 0.001)
//		{
//			auto q0x = spline_sp_x;
//			auto q0y = spline_sp_y;
//
//			// calculate spline end point
//			auto q1x = spline_ep_x;
//			auto q1y = spline_ep_y;
//
//			// spline start speed
//			auto v0x = v0 * cos(angle_before);
//			auto v0y = v1 * sin(angle_before);
//
//			// spline end speed
//			auto v1x = v0 * cos(angle_after);
//			auto v1y = v1 * sin(angle_after);
//
//			// spline start acceleration (it should be zero if it's line to line spline)
//			auto a0x = 0.0;
//			auto a0y = 0.0;
//
//			// spline end acceleration (it should be zero if it's line to line spline)
//			auto a1x = 0.0;
//			auto a1y = 0.0;
//
//			// travel time
//			auto Tx = current_spline_travel_time;
//			auto Ty = current_spline_travel_time;
//
//			int maximum_time_index = current_spline_travel_time * 1000.0;
//
//			mid_res = {};
//
//			// begin
//			for (int time_index = 0; time_index <= maximum_time_index; ++time_index)
//			{
//				// time
//				double t = (double)time_index / 1000.0;
//				// positions
//				auto xp = (q0x + t * v0x + (a0x * SIM_SQUARE(t)) / 2.0 -
//						   (SIM_FIFTH_POWER(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x + (SIM_SQUARE(Tx) * a0x) / 2.0 -
//						   (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) - (SIM_CUBE(t) * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x +
//						   (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx) + (SIM_FOURTH_POWER(t) *
//						   (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x + (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - SIM_SQUARE(Tx) * a1x)) / SIM_FOURTH_POWER(Tx));
//
//				auto yp = (q0y + t * v0y + (a0y * SIM_SQUARE(t)) / 2.0 -
//						   (SIM_FIFTH_POWER(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y + (SIM_SQUARE(Ty) * a0y) / 2.0 -
//						   (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) - (SIM_CUBE(t) * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y +
//						   (3.0 * SIM_SQUARE(t) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty) + (SIM_FOURTH_POWER(t) *
//						   (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y + (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - SIM_SQUARE(Ty) * a1y)) / SIM_FOURTH_POWER(Ty));
//
//				// velocities
//				auto xv = (v0x + a0x * t - (5.0 * SIM_FOURTH_POWER(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x +
//						   (SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) - (3.0 * SIM_SQUARE(t) * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x +
//						   (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx) + (4.0 * SIM_CUBE(t) * (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x +
//						   (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - SIM_SQUARE(Tx) * a1x)) / SIM_FOURTH_POWER(Tx));
//
//				auto yv = (v0y + a0y * t - (5.0 * SIM_FOURTH_POWER(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y +
//						   (SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) - (3.0 * SIM_SQUARE(t) * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y +
//						   (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty) + (4.0 * SIM_CUBE(t) * (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y +
//						   (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - SIM_SQUARE(Ty) * a1y)) / SIM_FOURTH_POWER(Ty));
//
//				// accelerations
//				auto xa = (a0x -
//						   (20.0 * SIM_CUBE(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x + ((SIM_SQUARE(Tx) * a0x) / 2.0) - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) +
//						   (12.0 * SIM_SQUARE(t) * (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x + ((3.0 * SIM_SQUARE(Tx) * a0x) / 2.0) - (SIM_SQUARE(Tx) * a1x))) / SIM_FOURTH_POWER(Tx) -
//						   (6.0 * t * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x + ((3.0 * SIM_SQUARE(Tx) * a0x) / 2.0) - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx));
//
//				auto ya = (a0y -
//						   (20.0 * SIM_CUBE(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y + ((SIM_SQUARE(Ty) * a0y) / 2.0) - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) +
//						   (12.0 * SIM_SQUARE(t) * (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y + ((3.0 * SIM_SQUARE(Ty) * a0y) / 2.0) - (SIM_SQUARE(Ty) * a1y))) / SIM_FOURTH_POWER(Ty) -
//						   (6.0 * t * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y + ((3.0 * SIM_SQUARE(Ty) * a0y) / 2.0) - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty));
//
//				// jerks
//				auto xj = ((24.0 * t * (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x + (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - SIM_SQUARE(Tx) * a1x)) / SIM_FOURTH_POWER(Tx) -
//						   (60.0 * SIM_SQUARE(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x + (SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) -
//						   (6.0 * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x + (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx));
//
//				auto yj = ((24.0 * t * (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y + (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - SIM_SQUARE(Ty) * a1y)) / SIM_FOURTH_POWER(Ty) -
//						   (60.0 * SIM_SQUARE(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y + (SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) -
//						   (6.0 * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y + (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty));
//
//				auto ttt_candidate = sqrt(SIM_SQUARE(xp - WTF_x) + SIM_SQUARE(yp - WTF_y));
//				if (ttt_candidate < min_ttt)
//				{
//					min_ttt = ttt_candidate;
//				}
//				// input.end_pose.position.z
//				mid_res.tool_path.push_back({ xp, yp, 0 });
//				mid_res.velocities.push_back({ xv, yv, 0.0 });
//				mid_res.accelerations.push_back({ xa, ya, 0.0 });
//				mid_res.jerks.push_back({ xj, yj, 0.0 });
//			} // end for X(t), X'(t), X''(t), X'''(t)
//
//			double sum = 0.0;
//
//			for (int i = 0; i <= maximum_time_index; ++i)
//			{
//				// TODO(batuhan): Why did the areas change, wtf?
//				double Area_PP1P2 = SIM_ONE_HALF * fabs(determinant({ mid_res.tool_path[i].x, mid_res.tool_path[i].y, 1, q0x, q0y, 1, WTF_x, WTF_y, 1 }));
//				double Area_PP2P3 = SIM_ONE_HALF * fabs(determinant({ mid_res.tool_path[i].x, mid_res.tool_path[i].y, 1, WTF_x, WTF_y, 1, q1x, q1y, 1 }));
//				double Area_PP3P1 = SIM_ONE_HALF * fabs(determinant({ mid_res.tool_path[i].x, mid_res.tool_path[i].y, 1, q1x, q1y, 1, q0x, q0y, 1 }));
//				double Area_P1P2P3 = SIM_ONE_HALF * fabs(determinant({ q0x, q0y, 1, WTF_x, WTF_y, 1, q1x, q1y, 1 }));
//				sum += ((Area_PP1P2 + Area_PP2P3 + Area_PP3P1) - Area_P1P2P3) / Area_P1P2P3;
//			} // end for sum
//
//			double vx_min = std::min_element(
//				mid_res.velocities.begin(), mid_res.velocities.end(),
//				[](auto const& a, auto const& b) { return a.x < b.x; })->x;
//
//			double vx_max = std::max_element(
//				mid_res.velocities.begin(), mid_res.velocities.end(),
//				[](auto const& a, auto const& b) { return a.x < b.x; })->x;
//
//			double vy_min = std::min_element(
//				mid_res.velocities.begin(), mid_res.velocities.end(),
//				[](auto const& a, auto const& b) { return a.y < b.y; })->y;
//
//			double vy_max = std::max_element(
//				mid_res.velocities.begin(), mid_res.velocities.end(),
//				[](auto const& a, auto const& b) { return a.y < b.y; })->y;
//
//			double ax_min = std::min_element(
//				mid_res.accelerations.begin(), mid_res.accelerations.end(),
//				[](auto const& a, auto const& b) { return a.x < b.x; })->x;
//
//			double ax_max = std::max_element(
//				mid_res.accelerations.begin(), mid_res.accelerations.end(),
//				[](auto const& a, auto const& b) { return a.x < b.x; })->x;
//
//			double ay_min = std::min_element(
//				mid_res.accelerations.begin(), mid_res.accelerations.end(),
//				[](auto const& a, auto const& b) { return a.y < b.y; })->y;
//
//			double ay_max = std::max_element(
//				mid_res.accelerations.begin(), mid_res.accelerations.end(),
//				[](auto const& a, auto const& b) { return a.y < b.y; })->y;
//
//			double jx_min = std::min_element(
//				mid_res.jerks.begin(), mid_res.jerks.end(),
//				[](auto const& a, auto const& b) { return a.x < b.x; })->x;
//
//			double jx_max = std::max_element(
//				mid_res.jerks.begin(), mid_res.jerks.end(),
//				[](auto const& a, auto const& b) { return a.x < b.x; })->x;
//
//			double jy_min = std::min_element(
//				mid_res.jerks.begin(), mid_res.jerks.end(),
//				[](auto const& a, auto const& b) { return a.y < b.y; })->y;
//
//			double jy_max = std::max_element(
//				mid_res.jerks.begin(), mid_res.jerks.end(),
//				[](auto const& a, auto const& b) { return a.y < b.y; })->y;
//
//			// TODO(batuhan): Is there a bug here, it actually should be maximum_time_index + 1;
//			// but this works somehow.
//			check_two = sum / (maximum_time_index);
//			check_one =
//			{
//				(vx_max - vmax) / vmax,
//				(vy_max - vmax) / vmax,
//				(-vx_min - vmax) / vmax,
//				(-vy_min - vmax) / vmax,
//				(ax_max - amax) / amax,
//				(ay_max - amax) / amax,
//				(-ax_min - amax) / amax,
//				(-ay_min - amax) / amax,
//				(jx_max - jmax) / jmax,
//				(jy_max - jmax) / jmax,
//				(-jx_min - jmax) / jmax,
//				(-jy_min - jmax) / jmax
//			};
//
//			auto check_one_max = *std::max_element(check_one.begin(), check_one.end());
//			auto check_one_min_spec = *std::min_element(check_one.begin() + 1, check_one.end(), [](auto const& a, auto const& b) { return fabs(a) < fabs(b); });
//
//			if (check_two < 1e-2 &&
//				check_one_max <= 0.0 &&
//				check_one_min_spec <= 0.3)
//			{
//				spline_fit = true;
//				break;
//			}
//			//} // end for ctf iteration
//			current_spline_travel_time += spline_travel_time_increment;
//
//			if (++iteration_number > 100)
//			{
//				break;
//			}
//
//			// TODO(batuhan): a0, a1 always 0?
//		} // end while(check)
//	}
//}
bool can_fit_spline(const Interpolation_Input& input, int input_type, const Interpolation_Input& next_input, int next_input_type, double cornering_tolerance)
{
	auto arc_length = [](Interpolation_Input _input)
	{
		glm::dvec3 start_vector = _input.start_pose.position - _input.center_point;
		glm::dvec3 end_vector = _input.end_pose.position - _input.center_point;
		double inner_product = glm::dot(start_vector, end_vector);
		double radius = glm::length(start_vector);
		double angle_rad = 0;	
		angle_rad = acos(inner_product / SIM_SQUARE(radius));
		return radius * abs(angle_rad);
	};

	double first_segment_length;
	double second_segment_length = 0;
	if(input_type == 0) 	// if first segment is line
	{
		first_segment_length = glm::length(input.end_pose.position - input.start_pose.position);
	}
	else if(input_type == 1) 	// if first segment is arc
	{
		first_segment_length = arc_length(input);
	}
	if(next_input_type == 0) 	// if second segment is line
	{
		second_segment_length = glm::length(next_input.end_pose.position - next_input.start_pose.position);
	}
	else if(next_input_type == 1) 	// if second segment is arc
	{
		second_segment_length = arc_length(next_input);
	}

	if (((first_segment_length) <= cornering_tolerance) 
			|| ((second_segment_length) <= cornering_tolerance))
	{
		return false;
	}
	return true;
}

std::tuple<unsigned long long, double>  SIM_trajectory_with_spline(const Interpolation_Input& input, int input_type, Interpolation_Input& next_input, int next_input_type, Interpolation_Output* result, const double cornering_tolerance)
{
	constexpr double spline_travel_time_increment = 0.001; // TODO enes: this was different
	constexpr double max_spline_travel_time = 1.0;

	// constexpr double cornering_tolerance_factor_increment = 0.2;
	constexpr auto JERK_LIMIT = 50000; 				// TODO enes: there is a jerk limit that overrides upper jerk limits !!! 
	unsigned long long  point_count = 0;

	if (!can_fit_spline(input, input_type, next_input, next_input_type, cornering_tolerance))
	{
		Interpolation_Input first_path_input{ input };
		first_path_input.final_feedrate = input.feedrate;

		std::cout << "cornering tolerance is larger than segment length!" << std::endl;
		if (input_type == 0)
		{
			point_count += SIM_linear_interpolation(first_path_input, result);
		}
		else if (input_type == 1)
		{
			point_count += SIM_circular_interpolation(first_path_input, result);
		}

		// result->NextInitialFeedRate = result->NextInitialFeedRate;  			//------>>>>>> Next initial feedrate is moved out  
		result->EndPoint = input.end_pose.position;
		return std::make_pair(point_count, 0); // do not fit spline and do not stop
	}

	glm::dvec3 spline_sp;
	glm::dvec3 spline_ep;

	if (input_type == 0) {
		spline_sp = SIM_find_spline_endpoint_linesegment(input, true, cornering_tolerance);
	}
	else if (input_type == 1)
	{
		spline_sp = SIM_find_spline_endpoint_arcsegment(input, false, cornering_tolerance);
	}

	if (next_input_type == 0)
	{
		spline_ep = SIM_find_spline_endpoint_linesegment(next_input, false, cornering_tolerance);
	}
	else if (next_input_type == 1)
	{
		spline_ep = SIM_find_spline_endpoint_arcsegment(next_input, true, cornering_tolerance);
	}

	glm::dvec3 corner_point = input.end_pose.position;

	double spline_sp_x, spline_sp_y;
	double spline_ep_x, spline_ep_y;
	double angle_before, angle_after;

	glm::dmat3 R;

	SIM_convert_2D_to_3D(spline_sp, spline_ep, corner_point, &spline_sp_x, &spline_sp_y, &spline_ep_x, &spline_ep_y, &angle_before, &angle_after, &R);

	glm::dvec3 n = glm::cross(spline_sp - corner_point, spline_ep - corner_point);
	glm::dvec3 n_spline_plane = n / glm::length(n);

	double vmax = input.feedrate;
	double amax = input.acceleration;
	double jmax = JERK_LIMIT;
	double v0 = vmax;
	double v1 = vmax;
	double V_proj_mag = 0.0;
	double Radius;

	const Interpolation_Input* ref_inp = nullptr;
	glm::dvec3* ref_point = nullptr;

	if (input_type == 1 && next_input_type == 0) // arc to line
	{
		ref_inp = &input;
		ref_point = &spline_sp;
	}
	else if (input_type == 0 && next_input_type == 1) // line to arc
	{
		ref_inp = &next_input;
		ref_point = &spline_ep;
	}

	if (ref_inp)
	{
		glm::dvec3 Oprime_start = ref_inp->start_pose.position - ref_inp->center_point;
		Radius = sqrt(dot(Oprime_start, Oprime_start));
		glm::dvec3 V1 = *ref_point - ref_inp->center_point;

		glm::dvec3 arc_segment_normal = { 0, 0, 1 };
		if (ref_inp->plane == Ref_Plane_Type::PLANE_XZ)
		{
			arc_segment_normal = { 0, 1, 0 };
		}
		else if (ref_inp->plane == Ref_Plane_Type::PLANE_YZ)
		{
			arc_segment_normal = { 1, 0, 0 };
		}
		glm::dvec3 V2 = glm::cross(V1/glm::length(V1), arc_segment_normal);
		if (ref_inp->direction == Arc_Segment_Direction::DIRECTION_CCW)
		{
			V2 = -V2;
		}

		glm::dvec3 V_proj = V2 - glm::dot((V2 - (*ref_point)), n_spline_plane) * n_spline_plane;
		V_proj_mag = sqrt(glm::dot(V_proj, V_proj));

		// TODO(batuhan): a0, a1 always 0?
	}

	Interpolation_Output mid_res;
	bool spline_fit = false;
	// TODO(batuhan): idk why x and y are always 0.
	double WTF_x = 0;
	double WTF_y = 0;

	int iteration_number = 0;
	double current_spline_travel_time = spline_travel_time_increment;
	double ki_bisection = 0.0;
	double kf_bisection = 1.0;
	double k_bisection = 1.0;
	std::vector<double> k_valid_iterations = {k_bisection};
	double k_error_bisection = 1.0;

	while (k_error_bisection > 0.1) 	// find the spline for optimum velocity // this is "dum" matlab
	{
		v0 = k_bisection * vmax;
		v1 = k_bisection * vmax;

		if(v0 < 0.3) 		// TODO enes: pay attention here. I think that this part is extra
		{
			break;
		}

		double check_two = 0.001; 		// TODO enes: this is 0.001 in matlab but was 0.01 here
		std::vector<double> check_one(11, 1.0);

		double min_ttt = std::numeric_limits<double>::infinity();
		current_spline_travel_time = spline_travel_time_increment;
		spline_fit = false;

		// spline fit
		while ((*std::max_element(check_one.begin(), check_one.end()) > 0.0) ||
			   ((check_two > 0.001) && (min_ttt < 0.001))) 		// TODO enes: here the check two tolerance is larger than matlab's
		{													// while start : find the spline for a velocity
			auto q0x = spline_sp_x;
			auto q0y = spline_sp_y;

			// calculate spline end point
			auto q1x = spline_ep_x;
			auto q1y = spline_ep_y;

			// spline start speed
			auto v0x = v0 * cos(angle_before);
			auto v0y = v1 * sin(angle_before);

			// spline end speed
			auto v1x = v0 * cos(angle_after);
			auto v1y = v1 * sin(angle_after);

			auto a0 = 0.0;
			auto a1 = 0.0;

			if (input_type == 0 && next_input_type == 0) // both segments are line
			{
				a0 = 0;
				a1 = 0;
			}
			else if (input_type == 1) // first segment is arc
			{
				a0 = (SIM_SQUARE(V_proj_mag * v0)) / Radius;
				a1 = 0;
				
			}
			else if (next_input_type == 1) // second segment is arc
			{
				a0 = 0;
				a1 = (SIM_SQUARE(V_proj_mag * v1)) / Radius;
			}
			// spline start acceleration (it should be zero if it's line to line spline)
			auto a0x = a0 * cos(angle_before);
			auto a0y = a0 * sin(angle_before);
			// spline end acceleration (it should be zero if it's line to line spline)
			auto a1x = a1 * cos(angle_after);
			auto a1y = a1 * sin(angle_after);

			// travel time
			auto Tx = current_spline_travel_time;
			auto Ty = current_spline_travel_time;

			int maximum_time_index = current_spline_travel_time * 1000.0;

			mid_res = {};

			// begin
			for (int time_index = 0; time_index <= maximum_time_index; ++time_index)
			{
				// time
				double t = (double)time_index / 1000.0;
				// positions
				auto xp = (q0x + t * v0x + (a0x * SIM_SQUARE(t)) / 2.0 -
						   (SIM_FIFTH_POWER(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x + (SIM_SQUARE(Tx) * a0x) / 2.0 -
						   (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) - (SIM_CUBE(t) * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x +
						   (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx) + (SIM_FOURTH_POWER(t) *
						   (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x + (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - SIM_SQUARE(Tx) * a1x)) / SIM_FOURTH_POWER(Tx));

				auto yp = (q0y + t * v0y + (a0y * SIM_SQUARE(t)) / 2.0 -
						   (SIM_FIFTH_POWER(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y + (SIM_SQUARE(Ty) * a0y) / 2.0 -
						   (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) - (SIM_CUBE(t) * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y +
						   (3.0 * SIM_SQUARE(t) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty) + (SIM_FOURTH_POWER(t) *
						   (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y + (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - SIM_SQUARE(Ty) * a1y)) / SIM_FOURTH_POWER(Ty));

				// velocities
				auto xv = (v0x + a0x * t - (5.0 * SIM_FOURTH_POWER(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x +
						   (SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) - (3.0 * SIM_SQUARE(t) * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x +
						   (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx) + (4.0 * SIM_CUBE(t) * (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x +
						   (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - SIM_SQUARE(Tx) * a1x)) / SIM_FOURTH_POWER(Tx));

				auto yv = (v0y + a0y * t - (5.0 * SIM_FOURTH_POWER(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y +
						   (SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) - (3.0 * SIM_SQUARE(t) * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y +
						   (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty) + (4.0 * SIM_CUBE(t) * (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y +
						   (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - SIM_SQUARE(Ty) * a1y)) / SIM_FOURTH_POWER(Ty));

				// accelerations
				auto xa = (a0x -
						   (20.0 * SIM_CUBE(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x + ((SIM_SQUARE(Tx) * a0x) / 2.0) - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) +
						   (12.0 * SIM_SQUARE(t) * (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x + ((3.0 * SIM_SQUARE(Tx) * a0x) / 2.0) - (SIM_SQUARE(Tx) * a1x))) / SIM_FOURTH_POWER(Tx) -
						   (6.0 * t * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x + ((3.0 * SIM_SQUARE(Tx) * a0x) / 2.0) - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx));

				auto ya = (a0y -
						   (20.0 * SIM_CUBE(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y + ((SIM_SQUARE(Ty) * a0y) / 2.0) - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) +
						   (12.0 * SIM_SQUARE(t) * (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y + ((3.0 * SIM_SQUARE(Ty) * a0y) / 2.0) - (SIM_SQUARE(Ty) * a1y))) / SIM_FOURTH_POWER(Ty) -
						   (6.0 * t * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y + ((3.0 * SIM_SQUARE(Ty) * a0y) / 2.0) - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty));

				// jerks
				auto xj = ((24.0 * t * (15.0 * q0x - 15.0 * q1x + 8.0 * Tx * v0x + 7.0 * Tx * v1x + (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - SIM_SQUARE(Tx) * a1x)) / SIM_FOURTH_POWER(Tx) -
						   (60.0 * SIM_SQUARE(t) * (6.0 * q0x - 6.0 * q1x + 3.0 * Tx * v0x + 3.0 * Tx * v1x + (SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_FIFTH_POWER(Tx) -
						   (6.0 * (10.0 * q0x - 10.0 * q1x + 6.0 * Tx * v0x + 4.0 * Tx * v1x + (3.0 * SIM_SQUARE(Tx) * a0x) / 2.0 - (SIM_SQUARE(Tx) * a1x) / 2.0)) / SIM_CUBE(Tx));

				auto yj = ((24.0 * t * (15.0 * q0y - 15.0 * q1y + 8.0 * Ty * v0y + 7.0 * Ty * v1y + (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - SIM_SQUARE(Ty) * a1y)) / SIM_FOURTH_POWER(Ty) -
						   (60.0 * SIM_SQUARE(t) * (6.0 * q0y - 6.0 * q1y + 3.0 * Ty * v0y + 3.0 * Ty * v1y + (SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_FIFTH_POWER(Ty) -
						   (6.0 * (10.0 * q0y - 10.0 * q1y + 6.0 * Ty * v0y + 4.0 * Ty * v1y + (3.0 * SIM_SQUARE(Ty) * a0y) / 2.0 - (SIM_SQUARE(Ty) * a1y) / 2.0)) / SIM_CUBE(Ty));

				auto ttt_candidate = sqrt(SIM_SQUARE(xp - WTF_x) + SIM_SQUARE(yp - WTF_y));
				if (ttt_candidate < min_ttt)
				{
					min_ttt = ttt_candidate;
				}
				// input.end_pose.position.z
				mid_res.tool_path.push_back({ xp, yp, 0 });
				mid_res.velocities.push_back({ xv, yv, 0.0 });
				mid_res.accelerations.push_back({ xa, ya, 0.0 });
				mid_res.jerks.push_back({ xj, yj, 0.0 });
			} // end for X(t), X'(t), X''(t), X'''(t)

			double sum = 0.0;

			for (int i = 0; i <= maximum_time_index; ++i)
			{
				// TODO(batuhan): Why did the areas change, wtf?
				double Area_PP1P2 = SIM_ONE_HALF * fabs(glm::determinant(glm::dmat3{ mid_res.tool_path[i].x, mid_res.tool_path[i].y, 1, q0x, q0y, 1, WTF_x, WTF_y, 1 }));
				double Area_PP2P3 = SIM_ONE_HALF * fabs(glm::determinant(glm::dmat3{ mid_res.tool_path[i].x, mid_res.tool_path[i].y, 1, WTF_x, WTF_y, 1, q1x, q1y, 1 }));
				double Area_PP3P1 = SIM_ONE_HALF * fabs(glm::determinant(glm::dmat3{ mid_res.tool_path[i].x, mid_res.tool_path[i].y, 1, q1x, q1y, 1, q0x, q0y, 1 }));
				double Area_P1P2P3 = SIM_ONE_HALF * fabs(glm::determinant(glm::dmat3{ q0x, q0y, 1, WTF_x, WTF_y, 1, q1x, q1y, 1 }));
				sum += ((Area_PP1P2 + Area_PP2P3 + Area_PP3P1) - Area_P1P2P3) / Area_P1P2P3;
			} // end for sum

			double vx_min = std::min_element(
				mid_res.velocities.begin(), mid_res.velocities.end(),
				[](auto const& a, auto const& b) { return a.x < b.x; })->x;

			double vx_max = std::max_element(
				mid_res.velocities.begin(), mid_res.velocities.end(),
				[](auto const& a, auto const& b) { return a.x < b.x; })->x;

			double vy_min = std::min_element(
				mid_res.velocities.begin(), mid_res.velocities.end(),
				[](auto const& a, auto const& b) { return a.y < b.y; })->y;

			double vy_max = std::max_element(
				mid_res.velocities.begin(), mid_res.velocities.end(),
				[](auto const& a, auto const& b) { return a.y < b.y; })->y;

			double ax_min = std::min_element(
				mid_res.accelerations.begin(), mid_res.accelerations.end(),
				[](auto const& a, auto const& b) { return a.x < b.x; })->x;

			double ax_max = std::max_element(
				mid_res.accelerations.begin(), mid_res.accelerations.end(),
				[](auto const& a, auto const& b) { return a.x < b.x; })->x;

			double ay_min = std::min_element(
				mid_res.accelerations.begin(), mid_res.accelerations.end(),
				[](auto const& a, auto const& b) { return a.y < b.y; })->y;

			double ay_max = std::max_element(
				mid_res.accelerations.begin(), mid_res.accelerations.end(),
				[](auto const& a, auto const& b) { return a.y < b.y; })->y;

			double jx_min = std::min_element(
				mid_res.jerks.begin(), mid_res.jerks.end(),
				[](auto const& a, auto const& b) { return a.x < b.x; })->x;

			double jx_max = std::max_element(
				mid_res.jerks.begin(), mid_res.jerks.end(),
				[](auto const& a, auto const& b) { return a.x < b.x; })->x;

			double jy_min = std::min_element(
				mid_res.jerks.begin(), mid_res.jerks.end(),
				[](auto const& a, auto const& b) { return a.y < b.y; })->y;

			double jy_max = std::max_element(
				mid_res.jerks.begin(), mid_res.jerks.end(),
				[](auto const& a, auto const& b) { return a.y < b.y; })->y;

			// TODO(batuhan): Is there a bug here, it actually should be maximum_time_index + 1;
			// but this works somehow.
			check_two = sum / (maximum_time_index);
			check_one =
			{
				(vx_max - vmax) / vmax,
				(vy_max - vmax) / vmax,
				(-vx_min - vmax) / vmax,
				(-vy_min - vmax) / vmax,
				(ax_max - amax) / amax,
				(ay_max - amax) / amax,
				(-ax_min - amax) / amax,
				(-ay_min - amax) / amax,
				(jx_max - jmax) / jmax,
				(jy_max - jmax) / jmax,
				(-jx_min - jmax) / jmax,
				(-jy_min - jmax) / jmax
			};

			auto check_one_max = *std::max_element(check_one.begin(), check_one.end());
			auto check_one_min_spec = *std::min_element(check_one.begin() + 1, check_one.end(), [](auto const& a, auto const& b) { return fabs(a) < fabs(b); });

			if (check_two < 1e-2 && 					// TODO enes: here is different from matlab
				check_one_max <= 0.0 &&
				check_one_min_spec <= 0.3)
			{
				spline_fit = true;
				break;
			}
			else
			{
				if (current_spline_travel_time > max_spline_travel_time)
				{
					spline_fit = false;
					break;
				}
			}
			//} // end for ctf iteration
			current_spline_travel_time += spline_travel_time_increment;


			// TODO(batuhan): a0, a1 always 0?
		} // while end : find the spline for a velocity

		if (spline_fit)
		{
			ki_bisection = k_bisection;
			k_error_bisection = fabs(k_valid_iterations.back() - k_bisection);
			k_valid_iterations.push_back(k_bisection);
		}
		else
		{
			kf_bisection = (k_bisection + kf_bisection) / 2;
		}

		k_bisection = (ki_bisection + kf_bisection) / 2;

		if (++iteration_number > 50)
		{
			break;
		}
	} // while end : find the spline for optimum velocity


	if (spline_fit)
	{
		Interpolation_Input linear_path{ input };
		linear_path.end_pose.position = spline_sp;
		linear_path.final_feedrate = v0;
		next_input.start_pose.position = spline_ep;
		
		if (input_type == 0) {

			point_count += SIM_linear_interpolation(linear_path, result);
		}
		else if (input_type == 1)
		{
			point_count += SIM_circular_interpolation(linear_path, result);
		}

		for (size_t i = 0; i < mid_res.tool_path.size(); ++i) {

			mid_res.tool_path[i] = corner_point + R * mid_res.tool_path[i];
			mid_res.velocities[i] = R * mid_res.velocities[i];
			mid_res.accelerations[i] = R * mid_res.accelerations[i];
			mid_res.jerks[i] = R * mid_res.jerks[i];
		}

		// next_input.start_pose.position = spline_ep;
		result->NextInitialFeedRate = SIM_MIN(v0, result->NextInitialFeedRate);
		result->EndPoint = spline_ep;

		// not pushing the first and last point of spline to prevent duplication ( [begin() + 1] to [end()] )
		result->tool_path.insert(result->tool_path.end(), mid_res.tool_path.begin() + 1, mid_res.tool_path.end());
		result->velocities.insert(result->velocities.end(), mid_res.velocities.begin() + 1, mid_res.velocities.end());
		result->accelerations.insert(result->accelerations.end(), mid_res.accelerations.begin() + 1, mid_res.accelerations.end());
		result->jerks.insert(result->jerks.end(), mid_res.jerks.begin() + 1, mid_res.jerks.end());

		point_count += mid_res.tool_path.size() - 1;
		std::cout << "Spline fit" << std::endl;
	}
	else // spline not fit
	{
		std::cout << "Spline could not fit !" << std::endl;
		if (input_type == 0)
		{
			point_count += SIM_linear_interpolation(input, result);
		}
		else if (input_type == 1) {

			point_count += SIM_circular_interpolation(input, result);
		}

		// result->NextInitialFeedRate = 0.0; 			//------>>>>>> Next initial feedrate is moved out
		result->EndPoint = input.end_pose.position;
	}

	return std::make_tuple(point_count, current_spline_travel_time);
}

glm::dvec3 SIM_find_spline_endpoint_linesegment(const Interpolation_Input& input, bool end, const double cornering_tolerance) {

	glm::dvec3 result;

	glm::dvec3 V1 = input.end_pose.position - input.start_pose.position;
	// isn't this just length(V1)?
	double L1 = sqrt(glm::dot(V1, V1));
	glm::dvec3 UV1 = V1 / L1;

	if (end){

		result = input.end_pose.position - UV1 * cornering_tolerance;
	}
	else {
		result = input.start_pose.position + UV1 * cornering_tolerance;
	}

	return result;
}

glm::dvec3 SIM_find_spline_endpoint_arcsegment(const Interpolation_Input& input, bool ELEVEN, const double cornering_tolerance)
{
	glm::dvec3 result;
	
	glm::dvec3 V1 = input.start_pose.position - input.center_point;
	glm::dvec3 V2 = input.end_pose.position - input.center_point;
	// isn't this just length(V1)?
	double L1 = sqrt(dot(V1, V1));
	double L2 = sqrt(dot(V2, V2));
	glm::dvec3 UV1 = V1 / L1;
	glm::dvec3 UV2 = V2 / L2;

	glm::dvec3 V2_cross_V1 = cross(V2, V1);
	//v3 V1_cross_V2 = cross(V1, V2);
	/* double theta = acos(dot(V2, V1)) / (length(V2) * length(V1)); */
	// NOTE(batuhan): divisor is V1_cross_V2 in the matlab
	glm::dvec3 normal_vec = V2_cross_V1 / (length(V2_cross_V1) + std::numeric_limits<double>::epsilon());

	if (V2_cross_V1.x == 0 &&
		V2_cross_V1.y == 0 &&
		V2_cross_V1.z == 0)
	{
		normal_vec.z = 1;
	}
	UV2 = -1.0 * glm::cross(UV1, normal_vec);

	glm::dmat3 R = { UV1, UV2, normal_vec };
	glm::dmat3 R_transpose = transpose(R);

	glm::dvec3 planar_A = R_transpose * V1;
	glm::dvec3 planar_B = R_transpose * V2;

	double theta_start = atan2(planar_A[1], planar_A[0]);
	double theta_end = atan2(planar_B[1], planar_B[0]);
	theta_start = SIM_normalize_radian(theta_start);
	theta_end = SIM_normalize_radian(theta_end);

	double xs, ys, xe, ye;
	double angle_before, angle_after;
	glm::dmat3 _unused;
	SIM_convert_2D_to_3D(input.start_pose.position, input.end_pose.position, input.center_point, &xs, &ys, &xe, &ye, &angle_before, &angle_after, &_unused);

	Circle_Intersection_Output xyout;
	if (ELEVEN)
	{
		xyout = SIM_circcirc(0, 0, L1, xs, ys, cornering_tolerance);
	}
	else
	{
		xyout = SIM_circcirc(0, 0, L1, xe, ye, cornering_tolerance);
	}

	// TODO(batuhan): this is 3x2 matrix in matlab
	glm::dmat3 intersection_points{ xyout.x3[0], xyout.y3[0], 0, xyout.x3[1], xyout.y3[1], 0, 0, 0, 0 };
	intersection_points = R * intersection_points;

	intersection_points[0] += input.center_point;
	intersection_points[1] += input.center_point;

	Interpolation_Output out;
	Interpolation_Input in{ input };
	in.initial_feedrate = 20; //TODO why these magic numbers?
	in.final_feedrate = 20;
	SIM_circular_interpolation(in, &out);

	SIM_ASSERT(out.tool_path.size() > 1);

	// TODO(batuhan): Ask about index 11, @splineend_pointFinder5 : 178, 179
	glm::dvec3 d1, d2;

	if (ELEVEN)	{
		SIM_ASSERT(out.tool_path.size() > 10);
		d1 = out.tool_path[10] - intersection_points[0];
		d2 = out.tool_path[10] - intersection_points[1];
	}
	else
	{
		SIM_ASSERT(out.tool_path.size() > 10);
		d1 = out.tool_path[out.tool_path.size() - 10] - intersection_points[0];
		d2 = out.tool_path[out.tool_path.size() - 10] - intersection_points[1];
	}

	double d1_l = glm::length(d1);
	double d2_l = glm::length(d2);

	if (d1_l < d2_l)
	{
		result = intersection_points[0];
	}
	else
	{
		result = intersection_points[1];
	}

	return result;
}