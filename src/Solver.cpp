#include "Solver.h"

#include <boost/graph/graphviz.hpp>
#include <fstream>

std::vector<std::string> show_edges = { "Left of", "Right of", "Front of", "Behind", "Above", "Under", "Close by", "Align with" };
std::vector<std::string> show_orientations = { "up", "down", "left", "right", "front", "back" };

Solver::Solver() : env(), model(env) {
    // Initialize solver-related data if needed
    hyperparameters = {1, 1, 1, 1};
	scalingFactor = 3;
}

Solver::~Solver() {}

bool Solver::has_path(const SceneGraph& g, VertexDescriptor start, VertexDescriptor target)
{
    std::vector<bool> visited_1(boost::num_vertices(g), false), visited_2(boost::num_vertices(g), false),
		visited_3(boost::num_vertices(g), false), visited_4(boost::num_vertices(g), false),
		visited_5(boost::num_vertices(g), false), visited_6(boost::num_vertices(g), false);
	return (dfs_check_path(g, start, target, LeftOf, visited_1) || dfs_check_path(g, start, target, RightOf, visited_1) || 
		dfs_check_path(g, start, target, FrontOf, visited_3) || dfs_check_path(g, start, target, Behind, visited_4) || 
		dfs_check_path(g, start, target, Above, visited_5) || dfs_check_path(g, start, target, Under, visited_6));
}

bool Solver::dfs_check_path(const SceneGraph& g, VertexDescriptor u, VertexDescriptor target, EdgeType required_type, std::vector<bool>& visited)
{
	if (g[u].id == g[target].id) return true;
	visited[g[u].id] = true;

	for (const auto& edge : boost::make_iterator_range(boost::out_edges(u, g))) {
		VertexDescriptor v = boost::target(edge, g);
		if (!visited[g[v].id] && g[edge].type == required_type) {
			if (dfs_check_path(g, v, target, required_type, visited)) {
				return true;
			}
		}
	}
	return false;
}

void Solver::addConstraints()
{
	int num_vertices = boost::num_vertices(g);
	int num_obstacles = obstacles.size();
	double M = boundary.size[0] + boundary.size[1] + boundary.size[2];
	std::vector<GRBVar> x_i(num_vertices), y_i(num_vertices), z_i(num_vertices),
		l_i(num_vertices), w_i(num_vertices), h_i(num_vertices);
	std::vector<std::vector<GRBVar>> sigma_L(num_vertices, std::vector<GRBVar>(num_vertices)),
		sigma_R(num_vertices, std::vector<GRBVar>(num_vertices)),
		sigma_F(num_vertices, std::vector<GRBVar>(num_vertices)),
		sigma_B(num_vertices, std::vector<GRBVar>(num_vertices)),
		sigma_U(num_vertices, std::vector<GRBVar>(num_vertices)),
		sigma_D(num_vertices, std::vector<GRBVar>(num_vertices)),
		L(num_vertices, std::vector<GRBVar>(num_vertices)),
		R(num_vertices, std::vector<GRBVar>(num_vertices)),
		F(num_vertices, std::vector<GRBVar>(num_vertices)),
		B(num_vertices, std::vector<GRBVar>(num_vertices)),
		sigma_oL(num_vertices, std::vector<GRBVar>(num_obstacles)),
		sigma_oR(num_vertices, std::vector<GRBVar>(num_obstacles)),
		sigma_oF(num_vertices, std::vector<GRBVar>(num_obstacles)),
		sigma_oB(num_vertices, std::vector<GRBVar>(num_obstacles)),
		sigma_oU(num_vertices, std::vector<GRBVar>(num_obstacles)),
		sigma_oD(num_vertices, std::vector<GRBVar>(num_obstacles));
	for (int i = 0; i < num_vertices; ++i) {
		x_i[i] = model.addVar(boundary.origin_pos[0], boundary.origin_pos[0] + boundary.size[0], 0.0, GRB_CONTINUOUS, "x_" + std::to_string(i));
		y_i[i] = model.addVar(boundary.origin_pos[1], boundary.origin_pos[1] + boundary.size[1], 0.0, GRB_CONTINUOUS, "y_" + std::to_string(i));
		l_i[i] = model.addVar(0.0, boundary.size[0], 0.0, GRB_CONTINUOUS, "l_" + std::to_string(i));
		w_i[i] = model.addVar(0.0, boundary.size[1], 0.0, GRB_CONTINUOUS, "w_" + std::to_string(i));
		if (!floorplan) {
			z_i[i] = model.addVar(boundary.origin_pos[2], boundary.origin_pos[2] + boundary.size[2], 0.0, GRB_CONTINUOUS, "z_" + std::to_string(i));
			h_i[i] = model.addVar(0.0, boundary.size[2], 0.0, GRB_CONTINUOUS, "h_" + std::to_string(i));
		}
	}
	// Inside Constraints & tolerance Constraint
	VertexIterator vi, vi_end;
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		model.addConstr(x_i[g[*vi].id] - l_i[g[*vi].id] / 2 >= boundary.origin_pos[0], "Inside_Object_" + std::to_string(g[*vi].id) + "_x_left");
		model.addConstr(x_i[g[*vi].id] + l_i[g[*vi].id] / 2 <= boundary.origin_pos[0] + boundary.size[0], "Inside_Object_" + std::to_string(g[*vi].id) + "_x_right");
		model.addConstr(y_i[g[*vi].id] - w_i[g[*vi].id] / 2 >= boundary.origin_pos[1], "Inside_Object_" + std::to_string(g[*vi].id) + "_y_back");
		model.addConstr(y_i[g[*vi].id] + w_i[g[*vi].id] / 2 <= boundary.origin_pos[1] + boundary.size[1], "Inside_Object_" + std::to_string(g[*vi].id) + "_y_front");
		if (!floorplan) {
			model.addConstr(z_i[g[*vi].id] - h_i[g[*vi].id] / 2 >= boundary.origin_pos[2], "Inside_Object_" + std::to_string(g[*vi].id) + "_z_bottom");
			model.addConstr(z_i[g[*vi].id] + h_i[g[*vi].id] / 2 <= boundary.origin_pos[2] + boundary.size[2], "Inside_Object_" + std::to_string(g[*vi].id) + "_z_top");
		}
		if (!g[*vi].pos_tolerance.empty() && !g[*vi].target_pos.empty()) {
			model.addConstr(x_i[g[*vi].id] >= g[*vi].target_pos[0] - g[*vi].pos_tolerance[0], "Pos_Tolerance_Object_" + std::to_string(g[*vi].id) + "_x_left");
			model.addConstr(x_i[g[*vi].id] <= g[*vi].target_pos[0] + g[*vi].pos_tolerance[0], "Pos_Tolerance_Object_" + std::to_string(g[*vi].id) + "_x_right");
			model.addConstr(y_i[g[*vi].id] >= g[*vi].target_pos[1] - g[*vi].pos_tolerance[1], "Pos_Tolerance_Object_" + std::to_string(g[*vi].id) + "_y_back");
			model.addConstr(y_i[g[*vi].id] <= g[*vi].target_pos[1] + g[*vi].pos_tolerance[1], "Pos_Tolerance_Object_" + std::to_string(g[*vi].id) + "_y_front");
			if (!floorplan) {
				model.addConstr(z_i[g[*vi].id] >= g[*vi].target_pos[2] - g[*vi].pos_tolerance[2], "Pos_Tolerance_Object_" + std::to_string(g[*vi].id) + "_z_bottom");
				model.addConstr(z_i[g[*vi].id] <= g[*vi].target_pos[2] + g[*vi].pos_tolerance[2], "Pos_Tolerance_Object_" + std::to_string(g[*vi].id) + "_z_top");
			}
		}
		if (!g[*vi].size_tolerance.empty() && !g[*vi].target_size.empty()) {
			model.addConstr(l_i[g[*vi].id] >= g[*vi].target_size[0] - g[*vi].size_tolerance[0], "Size_Tolerance_Object_" + std::to_string(g[*vi].id) + "_l_min");
			model.addConstr(l_i[g[*vi].id] <= g[*vi].target_size[0] + g[*vi].size_tolerance[0], "Size_Tolerance_Object_" + std::to_string(g[*vi].id) + "_l_max");
			model.addConstr(w_i[g[*vi].id] >= g[*vi].target_size[1] - g[*vi].size_tolerance[1], "Size_Tolerance_Object_" + std::to_string(g[*vi].id) + "_w_min");
			model.addConstr(w_i[g[*vi].id] <= g[*vi].target_size[1] + g[*vi].size_tolerance[1], "Size_Tolerance_Object_" + std::to_string(g[*vi].id) + "_w_max");
			if (!floorplan) {
				model.addConstr(h_i[g[*vi].id] >= g[*vi].target_size[2] - g[*vi].size_tolerance[2], "Size_Tolerance_Object_" + std::to_string(g[*vi].id) + "_h_min");
				model.addConstr(h_i[g[*vi].id] <= g[*vi].target_size[2] + g[*vi].size_tolerance[2], "Size_Tolerance_Object_" + std::to_string(g[*vi].id) + "_h_max");
			}
		}
	}
	// On floor Constraints
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		if (g[*vi].on_floor && !floorplan) {
			model.addConstr(z_i[g[*vi].id] == boundary.origin_pos[2] + h_i[g[*vi].id] / 2, "On_Floor_Object_" + std::to_string(g[*vi].id));
		}
	}
	// Hanging Constraints
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		if (g[*vi].hanging && !floorplan) {
			model.addConstr(z_i[g[*vi].id] == boundary.origin_pos[2] + boundary.size[2] - h_i[g[*vi].id] / 2, "Hanging_Object_" + std::to_string(g[*vi].id));
		}
	}
	// Adjacency Constraints
	EdgeIterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
		VertexDescriptor source = boost::source(*ei, g);
		VertexDescriptor target = boost::target(*ei, g);
		int ids = g[source].id, idt = g[target].id;
		switch(g[*ei].type)
		{
		case LeftOf:
			if (g[*ei].distance >= 0)
				model.addConstr(x_i[ids] + l_i[ids] / 2 <= x_i[idt] - l_i[idt] / 2, "Object_" + std::to_string(ids) + "_LeftOf_Object_" + std::to_string(idt));
			else
				model.addConstr(x_i[ids] + l_i[ids] / 2 == x_i[idt] - l_i[idt] / 2, "Object_" + std::to_string(ids) + "_LeftOf_Object_" + std::to_string(idt));
			break;
		case RightOf:
			if (g[*ei].distance >= 0)
				model.addConstr(x_i[ids] - l_i[ids] / 2 >= x_i[idt] + l_i[idt] / 2, "Object_" + std::to_string(ids) + "_RightOf_Object_" + std::to_string(idt));
			else
				model.addConstr(x_i[ids] - l_i[ids] / 2 == x_i[idt] + l_i[idt] / 2, "Object_" + std::to_string(ids) + "_RightOf_Object_" + std::to_string(idt));
			break;
		case Behind:
			if (g[*ei].distance >= 0)
				model.addConstr(y_i[ids] + w_i[ids] / 2 <= y_i[idt] - w_i[idt] / 2, "Object_" + std::to_string(ids) + "_Behind_Object_" + std::to_string(idt));
			else
				model.addConstr(y_i[ids] + w_i[ids] / 2 == y_i[idt] - w_i[idt] / 2, "Object_" + std::to_string(ids) + "_Behind_Object_" + std::to_string(idt));
			break;
		case FrontOf:
			if (g[*ei].distance >= 0)
				model.addConstr(y_i[ids] - w_i[ids] / 2 >= y_i[idt] + w_i[idt] / 2, "Object_" + std::to_string(ids) + "_FrontOf_Object_" + std::to_string(idt));
			else
				model.addConstr(y_i[ids] - w_i[ids] / 2 == y_i[idt] + w_i[idt] / 2, "Object_" + std::to_string(ids) + "_FrontOf_Object_" + std::to_string(idt));
			break;
		case Under:
			model.addConstr(z_i[ids] + h_i[ids] / 2 == z_i[idt] - h_i[idt] / 2, "Object_" + std::to_string(ids) + "_Under_Object_" + std::to_string(idt));
			break;
		case Above:
			model.addConstr(z_i[ids] - h_i[ids] / 2 == z_i[idt] + h_i[idt] / 2, "Object_" + std::to_string(ids) + "_Above_Object_" + std::to_string(idt));
			break;
		case CloseBy:
			L[ids][idt] = model.addVar(0, 1, 0, GRB_BINARY);
			R[ids][idt] = model.addVar(0, 1, 0, GRB_BINARY);
			F[ids][idt] = model.addVar(0, 1, 0, GRB_BINARY);
			B[ids][idt] = model.addVar(0, 1, 0, GRB_BINARY);
			model.addConstr(x_i[ids] - l_i[ids] / 2 <= x_i[idt] + l_i[idt] / 2 + M * (1 - R[ids][idt]), "Object_" + std::to_string(ids) + "_CloseBy_Object_" + std::to_string(idt) + "ieqa");
			model.addConstr(x_i[ids] + l_i[ids] / 2 >= x_i[idt] - l_i[idt] / 2 - M * (1 - L[ids][idt]), "Object_" + std::to_string(ids) + "_CloseBy_Object_" + std::to_string(idt) + "ieqb");
			model.addConstr(y_i[ids] - w_i[ids] / 2 <= y_i[idt] + w_i[idt] / 2 + M * (1 - F[ids][idt]), "Object_" + std::to_string(ids) + "_CloseBy_Object_" + std::to_string(idt) + "ieqc");
			model.addConstr(y_i[ids] + w_i[ids] / 2 >= y_i[idt] - w_i[idt] / 2 - M * (1 - B[ids][idt]), "Object_" + std::to_string(ids) + "_CloseBy_Object_" + std::to_string(idt) + "ieqd");
			model.addConstr(L[ids][idt] + R[ids][idt] + F[ids][idt] + B[ids][idt] <= 1, "Object_" + std::to_string(ids) + "_CloseBy_Object_" + std::to_string(idt) + "ieqe");
			break;
		case AlignWith:
			switch (g[*ei].align_edge)
			{
			case 0:
				model.addConstr(y_i[ids] - w_i[ids] / 2 == y_i[idt] - w_i[idt] / 2, "Object_" + std::to_string(ids) + "_AlignWith_Object_" + std::to_string(idt));
				break;
			case 1:
				model.addConstr(x_i[ids] + l_i[ids] / 2 == x_i[idt] + l_i[idt] / 2, "Object_" + std::to_string(ids) + "_AlignWith_Object_" + std::to_string(idt));
				break;
			case 2:
				model.addConstr(y_i[ids] + w_i[ids] / 2 == y_i[idt] + w_i[idt] / 2, "Object_" + std::to_string(ids) + "_AlignWith_Object_" + std::to_string(idt));
				break;
			case 3:
				model.addConstr(x_i[ids] - l_i[ids] / 2 == x_i[idt] - l_i[idt] / 2, "Object_" + std::to_string(ids) + "_AlignWith_Object_" + std::to_string(idt));
				break;
			case 4:
				model.addConstr(z_i[ids] - h_i[ids] / 2 == z_i[idt] + h_i[idt] / 2, "Object_" + std::to_string(ids) + "_AlignWith_Object_" + std::to_string(idt));
				break;
			case 5:
				model.addConstr(z_i[ids] + h_i[ids] / 2 == z_i[idt] - h_i[idt] / 2, "Object_" + std::to_string(ids) + "_AlignWith_Object_" + std::to_string(idt));
				break;
			default:break;
			}
			break;
		default:break;
		}
	}
	// Non overlap Constraints
	auto vi_start_end = boost::vertices(g);
	VertexIterator vj;
	for (vi = vi_start_end.first; vi != vi_start_end.second; ++vi) {
		for (vj = vi_start_end.first; vj != vi_start_end.second; ++vj) {
			if (g[*vi].id < g[*vj].id && !has_path(g, *vi, *vj) && !has_path(g, *vj, *vi)) {
				sigma_R[g[*vi].id][g[*vj].id] = model.addVar(0, 1, 0, GRB_BINARY);
				sigma_L[g[*vi].id][g[*vj].id] = model.addVar(0, 1, 0, GRB_BINARY);
				sigma_F[g[*vi].id][g[*vj].id] = model.addVar(0, 1, 0, GRB_BINARY);
				sigma_B[g[*vi].id][g[*vj].id] = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(x_i[g[*vi].id] - l_i[g[*vi].id] / 2 >= x_i[g[*vj].id] + l_i[g[*vj].id] / 2 -
					M * (1 - sigma_R[g[*vi].id][g[*vj].id]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id) + "R");
				model.addConstr(x_i[g[*vi].id] + l_i[g[*vi].id] / 2 <= x_i[g[*vj].id] - l_i[g[*vj].id] / 2 +
					M * (1 - sigma_L[g[*vi].id][g[*vj].id]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id) + "L");
				model.addConstr(y_i[g[*vi].id] - w_i[g[*vi].id] / 2 >= y_i[g[*vj].id] + w_i[g[*vj].id] / 2 - 
					M * (1 - sigma_F[g[*vi].id][g[*vj].id]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id) + "F");
				model.addConstr(y_i[g[*vi].id] + w_i[g[*vi].id] / 2 <= y_i[g[*vj].id] - w_i[g[*vj].id] / 2 + 
					M * (1 - sigma_B[g[*vi].id][g[*vj].id]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id) + "B");

				if (!floorplan) {
					sigma_U[g[*vi].id][g[*vj].id] = model.addVar(0, 1, 0, GRB_BINARY);
					sigma_D[g[*vi].id][g[*vj].id] = model.addVar(0, 1, 0, GRB_BINARY);
					model.addConstr(z_i[g[*vi].id] - h_i[g[*vi].id] / 2 >= z_i[g[*vj].id] + h_i[g[*vj].id] / 2 -
						M * (1 - sigma_U[g[*vi].id][g[*vj].id]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id) + "U");
					model.addConstr(z_i[g[*vi].id] + h_i[g[*vi].id] / 2 <= z_i[g[*vj].id] - h_i[g[*vj].id] / 2 +
						M * (1 - sigma_D[g[*vi].id][g[*vj].id]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id) + "D");
					model.addConstr(sigma_L[g[*vi].id][g[*vj].id] + sigma_R[g[*vi].id][g[*vj].id] +
						sigma_F[g[*vi].id][g[*vj].id] + sigma_B[g[*vi].id][g[*vj].id] +
						sigma_U[g[*vi].id][g[*vj].id] + sigma_D[g[*vi].id][g[*vj].id] >= 1, "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id));
				}
				else {
					model.addConstr(sigma_L[g[*vi].id][g[*vj].id] + sigma_R[g[*vi].id][g[*vj].id] +
						sigma_F[g[*vi].id][g[*vj].id] + sigma_B[g[*vi].id][g[*vj].id] >= 1, "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Object_" + std::to_string(g[*vj].id));
				}
			}
		}
	}
	// Obstacle Constraints
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		for (int i = 0; i < num_obstacles; ++i) {
			GRBLinExpr sigma_o = 0;
			sigma_oL[g[*vi].id][i] = model.addVar(0, 1, 0, GRB_BINARY);
			model.addConstr(x_i[g[*vi].id] + l_i[g[*vi].id] / 2 <= obstacles[i].pos[0] - obstacles[i].size[0] / 2 +
				M * (1 - sigma_oL[g[*vi].id][i]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i) + "L");
			sigma_o += sigma_oL[g[*vi].id][i];

			sigma_oR[g[*vi].id][i] = model.addVar(0, 1, 0, GRB_BINARY);
			model.addConstr(x_i[g[*vi].id] - l_i[g[*vi].id] / 2 >= obstacles[i].pos[0] + obstacles[i].size[0] / 2 -
				M * (1 - sigma_oR[g[*vi].id][i]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i) + "R");
			sigma_o += sigma_oR[g[*vi].id][i];

			sigma_oB[g[*vi].id][i] = model.addVar(0, 1, 0, GRB_BINARY);
			model.addConstr(y_i[g[*vi].id] + w_i[g[*vi].id] / 2 <= obstacles[i].pos[1] - obstacles[i].size[1] / 2 +
				M * (1 - sigma_oB[g[*vi].id][i]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i) + "B");	
			sigma_o += sigma_oB[g[*vi].id][i];
			sigma_oF[g[*vi].id][i] = model.addVar(0, 1, 0, GRB_BINARY);

			model.addConstr(y_i[g[*vi].id] - w_i[g[*vi].id] / 2 >= obstacles[i].pos[1] + obstacles[i].size[1] / 2 -
				M * (1 - sigma_oF[g[*vi].id][i]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i) + "F");
			sigma_o += sigma_oF[g[*vi].id][i];

			if (!floorplan) {
				sigma_oD[g[*vi].id][i] = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(z_i[g[*vi].id] + h_i[g[*vi].id] / 2 <= obstacles[i].pos[2] - obstacles[i].size[2] / 2 +
					M * (1 - sigma_oD[g[*vi].id][i]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i) + "D");
				sigma_o += sigma_oD[g[*vi].id][i];

				sigma_oU[g[*vi].id][i] = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(z_i[g[*vi].id] - h_i[g[*vi].id] / 2 >= obstacles[i].pos[2] + obstacles[i].size[2] / 2 -
					M * (1 - sigma_oU[g[*vi].id][i]), "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i) + "U");
				sigma_o += sigma_oU[g[*vi].id][i];
			}
			model.addConstr(sigma_o >= 1, "NonOverlap_Object_" + std::to_string(g[*vi].id) + "and_Obstacle_" + std::to_string(i));	
		}
	}
	// Boundary Constraints
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		if (g[*vi].boundary >= 0) {
			double x1 = boundary.points[g[*vi].boundary][0], x2 = boundary.points[(g[*vi].boundary + 1) % boundary.Orientations.size()][0];
			double y1 = boundary.points[g[*vi].boundary][1], y2 = boundary.points[(g[*vi].boundary + 1) % boundary.Orientations.size()][1];
			double x1_ = std::min(x1, x2), x2_ = std::max(x1, x2);
			double y1_ = std::min(y1, y2), y2_ = std::max(y1, y2);
			switch (boundary.Orientations[g[*vi].boundary])
			{
			case LEFT:
				model.addConstr(x_i[g[*vi].id] - l_i[g[*vi].id] / 2 == x1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Left_eq");
				// Here we assume that on boundary means at least half of length is on the wall
				model.addConstr(y_i[g[*vi].id] >= y1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Left_ieq");
				model.addConstr(y_i[g[*vi].id] <= y2_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Left_ieqq");
				break;
			case RIGHT:
				model.addConstr(x_i[g[*vi].id] + l_i[g[*vi].id] / 2 == x1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Right_eq");
				model.addConstr(y_i[g[*vi].id] >= y1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Right_ieq");
				model.addConstr(y_i[g[*vi].id] <= y2_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Right_ieqq");
				break;
			case FRONT:
				model.addConstr(y_i[g[*vi].id] + w_i[g[*vi].id] / 2 == y1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Front_eq");
				model.addConstr(x_i[g[*vi].id] >= x1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Front_ieq");
				model.addConstr(x_i[g[*vi].id] <= x2_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Front_ieqq");
				break;
			case BACK:
				model.addConstr(y_i[g[*vi].id] - w_i[g[*vi].id] / 2 == y1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Back_eq");
				model.addConstr(x_i[g[*vi].id] >= x1_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Back_ieq");
				model.addConstr(x_i[g[*vi].id] <= x2_, "Boundary_Object_" + std::to_string(g[*vi].id) + "_Back_ieqq");
				break;
			default:break;
			}
		}
	}
	// Floor Plan Constraints :AREA
	if (floorplan)
	{
		GRBQuadExpr total_area = 0;
		double unuse_area = boundary.size[0] * boundary.size[1];
		for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
			total_area += l_i[g[*vi].id] * w_i[g[*vi].id];
		}
		for (int i = 0; i < num_obstacles; ++i) {
			unuse_area -= obstacles[i].size[0] * obstacles[i].size[1];
		}
		model.addQConstr(total_area == unuse_area, "Area_Constraint_for_FloorPlan");
	}
	// Corner Constraints
	std::vector<std::vector<GRBVar>> cor(num_vertices);
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		switch (g[*vi].corner)
		{
			case BOTTOMLEFT:
			{
				int num = boundary.BLcorner.size();
				cor[g[*vi].id].resize(num);
				GRBLinExpr posx, posy, cors;
				for (int i = 0; i < num; ++i) {
					cor[g[*vi].id][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					posx += boundary.points[boundary.BLcorner[i]][0] * cor[g[*vi].id][i];
					posy += boundary.points[boundary.BLcorner[i]][1] * cor[g[*vi].id][i];
					cors += cor[g[*vi].id][i];
				}
				model.addConstr(cors == 1, "BottomLeft_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqa");
				model.addConstr(x_i[g[*vi].id] - l_i[g[*vi].id] / 2 == posx, "BottomLeft_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqb");
				model.addConstr(y_i[g[*vi].id] - w_i[g[*vi].id] / 2 == posy, "BottomLeft_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqc");
				break;
			}
			case BOTTOMRIGHT:
			{
				int num = boundary.BRcorner.size();
				cor[g[*vi].id].resize(num);
				GRBLinExpr posx, posy, cors;
				for (int i = 0; i < num; ++i) {
					cor[g[*vi].id][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					posx += boundary.points[boundary.BRcorner[i]][0] * cor[g[*vi].id][i];
					posy += boundary.points[boundary.BRcorner[i]][1] * cor[g[*vi].id][i];
					cors += cor[g[*vi].id][i];
				}
				model.addConstr(cors == 1, "BottomRight_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqa");
				model.addConstr(x_i[g[*vi].id] + l_i[g[*vi].id] / 2 == posx, "BottomRight_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqb");
				model.addConstr(y_i[g[*vi].id] - w_i[g[*vi].id] / 2 == posy, "BottomRight_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqc");
				break;
			}
			case TOPLEFT:
			{
				int num = boundary.TLcorner.size();
				cor[g[*vi].id].resize(num);
				GRBLinExpr posx, posy, cors;
				for (int i = 0; i < num; ++i) {
					cor[g[*vi].id][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					posx += boundary.points[boundary.TLcorner[i]][0] * cor[g[*vi].id][i];
					posy += boundary.points[boundary.TLcorner[i]][1] * cor[g[*vi].id][i];
					cors += cor[g[*vi].id][i];
				}
				model.addConstr(cors == 1, "TopLeft_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqa");
				model.addConstr(x_i[g[*vi].id] - l_i[g[*vi].id] / 2 == posx, "TopLeft_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqb");
				model.addConstr(y_i[g[*vi].id] + w_i[g[*vi].id] / 2 == posy, "TopLeft_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqc");
				break;
			}
			case TOPRIGHT:
			{
				int num = boundary.TRcorner.size();
				cor[g[*vi].id].resize(num);
				GRBLinExpr posx, posy, cors;
				for (int i = 0; i < num; ++i) {
					cor[g[*vi].id][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					posx += boundary.points[boundary.TRcorner[i]][0] * cor[g[*vi].id][i];
					posy += boundary.points[boundary.TRcorner[i]][1] * cor[g[*vi].id][i];
					cors += cor[g[*vi].id][i];
				}
				model.addConstr(cors == 1, "TopRight_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqa");
				model.addConstr(x_i[g[*vi].id] + l_i[g[*vi].id] / 2 == posx, "TopRight_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqb");
				model.addConstr(y_i[g[*vi].id] + w_i[g[*vi].id] / 2 == posy, "TopRight_Corner_of_Object_" + std::to_string(g[*vi].id) + "eqc");
				break;
			}
			default: break;
		}
	}
	
	
	// Objective Function
	// Notice that hyperparameters are the weights of area, size error, position error, adjacency error.
	GRBQuadExpr obj = hyperparameters[0];
	for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
		bool area_flag = true;
		boost::graph_traits<SceneGraph>::out_edge_iterator e_out, e_end;
		for (boost::tie(e_out, e_end) = boost::out_edges(*vi, g); e_out != e_end; ++e_out) {
			if (g[*e_out].type == Above || g[*e_out].type == Under) {
				area_flag = false;
				break;
			}
		}
		if (area_flag)
			obj -= hyperparameters[0] * l_i[g[*vi].id] * w_i[g[*vi].id] / boundary.size[0] / boundary.size[1];
		if (!g[*vi].target_size.empty()) {
			obj += hyperparameters[1] * (l_i[g[*vi].id] - g[*vi].target_size[0]) * (l_i[g[*vi].id] - g[*vi].target_size[0]) / boundary.size[0] / boundary.size[0];
			obj += hyperparameters[1] * (w_i[g[*vi].id] - g[*vi].target_size[1]) * (w_i[g[*vi].id] - g[*vi].target_size[1]) / boundary.size[1] / boundary.size[1];
			if (!floorplan)
				obj += hyperparameters[1] * (h_i[g[*vi].id] - g[*vi].target_size[2]) * (h_i[g[*vi].id] - g[*vi].target_size[2]) / boundary.size[2] / boundary.size[2];
		}
		if (!g[*vi].target_pos.empty()) {
			obj += hyperparameters[2] * (x_i[g[*vi].id] - g[*vi].target_pos[0]) * (x_i[g[*vi].id] - g[*vi].target_pos[0]) / boundary.size[0] / boundary.size[0];
			obj += hyperparameters[2] * (y_i[g[*vi].id] - g[*vi].target_pos[1]) * (y_i[g[*vi].id] - g[*vi].target_pos[1]) / boundary.size[1] / boundary.size[1];
			if (!floorplan)
				obj += hyperparameters[2] * (z_i[g[*vi].id] - g[*vi].target_pos[2]) * (z_i[g[*vi].id] - g[*vi].target_pos[2]) / boundary.size[2] / boundary.size[2];
		}
	}
	for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
		VertexDescriptor source = boost::source(*ei, g);
		VertexDescriptor target = boost::target(*ei, g);
		std::vector<double> offset = g[*ei].xyoffset;
		if (offset.empty())
			offset = { 0, 0 };
		if (g[*ei].distance >= 0) {
			switch (g[*ei].type)
			{
			case LeftOf:
				obj += hyperparameters[3] * (x_i[g[target].id] - l_i[g[target].id] / 2 - x_i[g[source].id] - l_i[g[source].id] / 2 - g[*ei].distance) *
					(x_i[g[target].id] - l_i[g[target].id] / 2 - x_i[g[source].id] - l_i[g[source].id] / 2 - g[*ei].distance) / boundary.size[0] / boundary.size[0];
				break;
			case RightOf:
				obj += hyperparameters[3] * (x_i[g[source].id] - l_i[g[source].id] / 2 - x_i[g[target].id] - l_i[g[target].id] / 2 - g[*ei].distance) *
					(x_i[g[source].id] - l_i[g[source].id] / 2 - x_i[g[target].id] - l_i[g[target].id] / 2 - g[*ei].distance) / boundary.size[0] / boundary.size[0];
				break;
			case Behind:
				obj += hyperparameters[3] * (y_i[g[target].id] - w_i[g[target].id] / 2 - y_i[g[source].id] - w_i[g[source].id] / 2 - g[*ei].distance) *
					(y_i[g[target].id] - w_i[g[target].id] / 2 - y_i[g[source].id] - w_i[g[source].id] / 2 - g[*ei].distance) / boundary.size[1] / boundary.size[1];
				break;
			case FrontOf:
				obj += hyperparameters[3] * (y_i[g[source].id] - w_i[g[source].id] / 2 - y_i[g[target].id] - w_i[g[target].id] / 2 - g[*ei].distance) *
					(y_i[g[source].id] - w_i[g[source].id] / 2 - y_i[g[target].id] - w_i[g[target].id] / 2 - g[*ei].distance) / boundary.size[1] / boundary.size[1];
				break;
			default:break;
			}
		}
		if (g[*ei].type == Above || g[*ei].type == Under || g[*ei].type == CloseBy) {
			obj += hyperparameters[3] * (x_i[g[source].id] - x_i[g[target].id] - offset[0]) * (x_i[g[source].id] - x_i[g[target].id] - offset[0]) / boundary.size[0] / boundary.size[0];
			obj += hyperparameters[3] * (y_i[g[source].id] - y_i[g[target].id] - offset[1]) * (y_i[g[source].id] - y_i[g[target].id] - offset[1]) / boundary.size[1] / boundary.size[1];
		}
	}
	model.setObjective(obj, GRB_MINIMIZE);
}

void Solver::optimizeModel()
{
    try {
        model.set(GRB_DoubleParam_TimeLimit, 10);
		if (floorplan)
        	model.set(GRB_DoubleParam_MIPGap, 0.11);
		else
			model.set(GRB_DoubleParam_MIPGap, 0.01);
		model.set(GRB_IntParam_MIPFocus, 1);
        model.set(GRB_IntParam_Method, 2);
        model.set(GRB_DoubleParam_BarConvTol, 1e-4);
        model.set(GRB_IntParam_Cuts, 2);
        model.set(GRB_IntParam_Presolve, 0);
        model.optimize();
		int iter = 0;
        while (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE && iter < 3) {
            //model.computeIIS();
            //model.write("model.ilp");
            //std::cout << "Infeasible constraints written to 'model.ilp'" << std::endl;
			std::cout << "Model is infeasible. Calling IIS computation..." << std::endl;
            handleInfeasibleModel();
			iter++;
        }

		if (graphProcessor.conflict_info.empty()) {
			GRBVar* vars = model.getVars();
        	int numVars = model.get(GRB_IntAttr_NumVars);
        	for (auto i = 0; i < numVars; ++i) {
        	    std::string varName = vars[i].get(GRB_StringAttr_VarName);
        	    double varValue = vars[i].get(GRB_DoubleAttr_X);
        	    std::cout << "Variable " << varName << ": Value = " << varValue << std::endl;
        	}
        	VertexIterator vi1, vi_end1;
			for (boost::tie(vi1, vi_end1) = boost::vertices(g); vi1 != vi_end1; ++vi1) {
        	    g[*vi1].pos = { 0, 0, 0 };
        	    g[*vi1].size = { 0, 0, 0 };
				std::string varName = "x_" + std::to_string(g[*vi1].id);
				double varValue = model.getVarByName(varName).get(GRB_DoubleAttr_X);
				g[*vi1].pos[0] = varValue;

				varName = "y_" + std::to_string(g[*vi1].id);
				varValue = model.getVarByName(varName).get(GRB_DoubleAttr_X);
				g[*vi1].pos[1] = varValue;

				varName = "l_" + std::to_string(g[*vi1].id);
				varValue = model.getVarByName(varName).get(GRB_DoubleAttr_X);
				g[*vi1].size[0] = varValue;

				varName = "w_" + std::to_string(g[*vi1].id);
				varValue = model.getVarByName(varName).get(GRB_DoubleAttr_X);
				g[*vi1].size[1] = varValue;

        	    if (!floorplan) {
        	        varName = "z_" + std::to_string(g[*vi1].id);
        	        varValue = model.getVarByName(varName).get(GRB_DoubleAttr_X);
        	        g[*vi1].pos[2] = varValue;

        	        varName = "h_" + std::to_string(g[*vi1].id);
        	        varValue = model.getVarByName(varName).get(GRB_DoubleAttr_X);
        	        g[*vi1].size[2] = varValue;
        	    }
        	    else {
					g[*vi1].pos[2] = g[*vi1].target_size[2] / 2;
					g[*vi1].size[2] = g[*vi1].target_size[2];
        	    }
        	}
        	std::cout << "Value of objective function: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
		}
    }
    catch (GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...) {
        std::cout << "Exception during optimization" << std::endl;
    }
	

    VertexIterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
        std::cout << "Vertex " << g[*vi].id << " (" << g[*vi].label << ")" <<
            " Boundary Constraint: " << g[*vi].boundary <<
            " Orientation: " << graphProcessor.orientationnames[g[*vi].orientation] <<
            " Target Position: " << g[*vi].target_pos.size() <<
            " Target Size: " << g[*vi].target_size.size() << std::endl;
    }

    EdgeIterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        std::cout << "Edge (" << g[boost::source(*ei, g)].label << " -> "
            << g[boost::target(*ei, g)].label << ") "
            << ", Type: " << graphProcessor.edgenames[g[*ei].type] << std::endl;
    }
}

void Solver::saveGraph()
{
    std::ofstream file_in("graph_in.dot");
    if (!file_in.is_open()) {
        std::cerr << "Failed to open file for writing: graph_in.dot" << std::endl;
    } else {
        boost::write_graphviz(file_in, inputGraph, vertex_writer_in<SceneGraph::vertex_descriptor>(inputGraph), 
            edge_writer<SceneGraph::edge_descriptor>(inputGraph));
    }
    
    std::ofstream file_out("graph_out.dot");
    if (!file_out.is_open()) {
        std::cerr << "Failed to open file for writing: graph_out.dot" << std::endl;
    } else {
        boost::write_graphviz(file_out, g, vertex_writer_out<SceneGraph::vertex_descriptor>(g),
            edge_writer<SceneGraph::edge_descriptor>(g));
    }
	model.write("model.lp");

	try
    {
        std::ifstream ifs(inputpath);
        if (!ifs.is_open())
            std::cerr << "Failed to open input JSON file: " << inputpath << std::endl;
        
        nlohmann::json j;
        ifs >> j;
        ifs.close();

		if (!graphProcessor.conflict_info.empty()) {
			j["conflict_info"] = graphProcessor.conflict_info;
			j["plan_info"] = {};
			for (auto i = 0; i < graphProcessor.plan_info.size(); ++i)
				j["plan_info"].push_back(graphProcessor.plan_info[i]);
		}
		else {
			j["conflict_info"] = "";
			j["plan_info"] = {};
			for (auto i = 0; i < j["vertices"].size(); ++i) {
				j["vertices"][i]["position"] = {
					g[boost::vertex(i, g)].pos[0],
					g[boost::vertex(i, g)].pos[1],
					g[boost::vertex(i, g)].pos[2]
				};
				j["vertices"][i]["size"] = {
					g[boost::vertex(i, g)].size[0],
					g[boost::vertex(i, g)].size[1],
					g[boost::vertex(i, g)].size[2]
				};
			}
		}

		std::string outputpath = "output.json";
        std::ofstream ofs(outputpath);
        if (!ofs.is_open())
        {
            std::cerr << "Failed to open output JSON file: " << outputpath << std::endl;
            return;
        }
        ofs << j.dump(4) << std::endl;
        ofs.close();
        std::cout << "JSON file has been updated and saved to: " << outputpath << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception occurred: " << e.what() << std::endl;
    }
}

void Solver::solve()
{
	if (inputGraph.m_vertices.empty()) {
		std::cerr << "Scene Graph is empty!" << std::endl;
	}
	else if (!graphProcessor.conflict_info.empty()) {
		std::cout << graphProcessor.conflict_info << std::endl;
		std::cerr << "Conflict Constraints Found" << std::endl;
	}
	else {
		clearModel();
		addConstraints();
    	optimizeModel();
	}
	saveGraph();
}

void Solver::readSceneGraph(const std::string& path, float wallwidth)
{
	reset();
	inputpath = path;
	// Read JSON file
	std::ifstream file(path);
    nlohmann::json scene_graph_json;
    file >> scene_graph_json;

	floorplan = scene_graph_json["floorplan"];

    // Parse JSON to set boundary
    boundary.origin_pos = scene_graph_json["boundary"]["origin_pos"].get<std::vector<double>>();
    boundary.size = scene_graph_json["boundary"]["size"].get<std::vector<double>>();
    boundary.points = scene_graph_json["boundary"]["points"].get<std::vector<std::vector<double>>>();
	// calculate orientations
    boundary.Orientations = std::vector<Orientation>(boundary.points.size(), FRONT);
	
	if (!floorplan)
	{
		boundary.origin_pos[0] += wallwidth / 2; boundary.origin_pos[1] += wallwidth / 2;
		boundary.size[0] -= wallwidth; boundary.size[1] -= wallwidth;
		ClipperLib::Path boundaryp;
		ClipperLib::Paths solu;
		for (const auto& point : boundary.points) {
			boundaryp.push_back({ (int)(point[0] * std::pow(10, scalingFactor)), (int)(point[1] * std::pow(10, scalingFactor)) });
		}
		ClipperLib::ClipperOffset co;
		co.AddPath(boundaryp, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);
		co.Execute(solu, -(wallwidth / 2) * std::pow(10, scalingFactor));
		std::vector<std::vector<double>> new_points = {};
		for (const auto& point : solu[0]) {
			new_points.push_back({point.X / std::pow(10, scalingFactor), point.Y / std::pow(10, scalingFactor)});
		}
		int idx_ = 0;
		for (auto i = 0; i < new_points.size(); ++i)
		{
			if (std::fabs(std::fabs(new_points[i][0] - boundary.points[0][0]) - wallwidth / 2) < 1e-3 
			&& std::fabs(std::fabs(new_points[i][1] - boundary.points[0][1]) - wallwidth / 2) < 1e-3)
			{
				idx_ = i;
				break;
			}
		}
		for (auto i = 0; i < new_points.size(); ++i)
			boundary.points[i] = new_points[(i + idx_) % new_points.size()];
	}
	
	for (auto i = 0; i < boundary.points.size(); ++i) {
        // Get the current edge
        std::vector<double> p1 = boundary.points[i];
        std::vector<double> p2 = boundary.points[(i + 1) % boundary.points.size()];

        // Calculate direction vector
        std::vector<double> direction = {p2[0] - p1[0], p2[1] - p1[1]};

        // Determine normal vector (rotate direction by 90 degrees)
        std::vector<double> normal = {direction[1], -direction[0]};

        // Determine orientation
        if (normal[0] > 0) {
            boundary.Orientations[i] = RIGHT;
        } else if (normal[0] < 0) {
            boundary.Orientations[i] = LEFT;
        } else if (normal[1] > 0) {
            boundary.Orientations[i] = FRONT;
        } else {
            boundary.Orientations[i] = BACK;
        }
    }
	for (int i = 0; i < boundary.Orientations.size(); ++i) {
		Orientation o1 = boundary.Orientations[i];
		Orientation o2 = boundary.Orientations[(i + 1) % boundary.Orientations.size()];
		if (o1 == LEFT && o2 == BACK)
			boundary.BLcorner.push_back((i + 1) % boundary.Orientations.size());
		else if (o1 == BACK && o2 == RIGHT)
			boundary.BRcorner.push_back((i + 1) % boundary.Orientations.size());
		else if (o1 == RIGHT && o2 == FRONT)
			boundary.TRcorner.push_back((i + 1) % boundary.Orientations.size());
		else if (o1 == FRONT && o2 == LEFT)
			boundary.TLcorner.push_back((i + 1) % boundary.Orientations.size());
	}
    // set obstacles from boundary and doors/windows
	for (const auto& door : scene_graph_json["doors"]) {
		Doors d;
		d.orientation = door["orientation"];
		d.pos = door["pos"].get<std::vector<double>>();
		d.size = door["size"].get<std::vector<double>>();
		doors.push_back(d);
		Obstacles door_obstacle;
		switch (d.orientation) {
		case FRONT:
			door_obstacle.pos = {d.pos[0], d.pos[1] - d.size[0] / 2, d.pos[2]};
			door_obstacle.size = {d.size[0], d.size[0], d.size[2]};
			break;
		case BACK:
			door_obstacle.pos = {d.pos[0], d.pos[1] + d.size[0] / 2, d.pos[2]};
			door_obstacle.size = {d.size[0], d.size[0], d.size[2]};
			break;
		case LEFT:
			door_obstacle.pos = {d.pos[0] + d.size[1] / 2, d.pos[1], d.pos[2]};
			door_obstacle.size = {d.size[1], d.size[1], d.size[2]};
			break;
		case RIGHT:
			door_obstacle.pos = {d.pos[0] - d.size[1] / 2, d.pos[1], d.pos[2]};
			door_obstacle.size = {d.size[1], d.size[1], d.size[2]};
			break;
		default:break;
		}
		obstacles.push_back(door_obstacle);
	}
	for (const auto& window : scene_graph_json["windows"]) {
		Windows w;
		w.orientation = window["orientation"];
		w.pos = window["pos"].get<std::vector<double>>();
		w.size = window["size"].get<std::vector<double>>();
		windows.push_back(w);
		Obstacles window_obstacle;
		switch (w.orientation) {
		case FRONT:
			window_obstacle.pos = {w.pos[0], w.pos[1] - w.size[0] / 2, w.pos[2]};
			window_obstacle.size = {w.size[0], w.size[0], w.size[2]};
			break;
		case BACK:
			window_obstacle.pos = {w.pos[0], w.pos[1] + w.size[0] / 2, w.pos[2]};
			window_obstacle.size = {w.size[0], w.size[0], w.size[2]};
			break;
		case LEFT:
			window_obstacle.pos = {w.pos[0] + w.size[1] / 2, w.pos[1], w.pos[2]};
			window_obstacle.size = {w.size[1], w.size[1], w.size[2]};
			break;
		case RIGHT:
			window_obstacle.pos = {w.pos[0] - w.size[1] / 2, w.pos[1], w.pos[2]};
			window_obstacle.size = {w.size[1], w.size[1], w.size[2]};
			break;
		default:break;
		}
		obstacles.push_back(window_obstacle);
	}
    for (const auto& obs : scene_graph_json["obstacles"]) {
        Obstacles obstacle;
        obstacle.pos = obs["pos"].get<std::vector<double>>();
        obstacle.size = obs["size"].get<std::vector<double>>();
        obstacles.push_back(obstacle);
    }
	// set obstacles from boundary
	ClipperLib::Path boundary_path, bounding_box;
	bounding_box << ClipperLib::IntPoint(
    	std::round(boundary.origin_pos[0] * std::pow(10, scalingFactor)),
    	std::round(boundary.origin_pos[1] * std::pow(10, scalingFactor)))
		<< ClipperLib::IntPoint(
    	std::round((boundary.origin_pos[0] + boundary.size[0]) * std::pow(10, scalingFactor)),
    	std::round(boundary.origin_pos[1] * std::pow(10, scalingFactor)))
		<< ClipperLib::IntPoint(
		std::round((boundary.origin_pos[0] + boundary.size[0]) * std::pow(10, scalingFactor)),
		std::round((boundary.origin_pos[1] + boundary.size[1]) * std::pow(10, scalingFactor)))
		<< ClipperLib::IntPoint(
		std::round(boundary.origin_pos[0] * std::pow(10, scalingFactor)),
		std::round((boundary.origin_pos[1] + boundary.size[1]) * std::pow(10, scalingFactor)));
	for (const auto& point : boundary.points) {
		boundary_path.push_back({ (int)(point[0] * std::pow(10, scalingFactor)), (int)(point[1] * std::pow(10, scalingFactor)) });
	}
	ClipperLib::Clipper clipper;
    clipper.AddPath(bounding_box, ClipperLib::ptSubject, true);
    clipper.AddPath(boundary_path, ClipperLib::ptClip, true);
    ClipperLib::Paths solution;
    clipper.Execute(ClipperLib::ctDifference, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	for (const auto& path : solution) {
        Obstacles ob_from_boundary;
		double minX = std::min({path[0].X, path[1].X, path[2].X, path[3].X});
		double maxX = std::max({path[0].X, path[1].X, path[2].X, path[3].X});
		double minY = std::min({path[0].Y, path[1].Y, path[2].Y, path[3].Y});
		double maxY = std::max({path[0].Y, path[1].Y, path[2].Y, path[3].Y});
		ob_from_boundary.pos = {(minX + maxX) / 2.0 / std::pow(10, scalingFactor), (minY + maxY)/ 2.0 / std::pow(10, scalingFactor), boundary.size[2] / 2.0};
		ob_from_boundary.size = {(maxX - minX) / std::pow(10, scalingFactor), (maxY - minY)/ std::pow(10, scalingFactor), boundary.size[2]};
		if (!floorplan) {
			ob_from_boundary.size[0] += wallwidth;
			ob_from_boundary.size[1] += wallwidth;
		}
		obstacles.push_back(ob_from_boundary);
    }
    // Parse JSON to set vertices
	// std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_int_distribution<> dis(0, 1);
    for (const auto& vertex : scene_graph_json["vertices"]) {
        VertexProperties vp;
        vp.label = vertex["label"];
        vp.id = vertex["id"];
        vp.boundary = vertex["boundary"];
        vp.on_floor = vertex["on_floor"];
		vp.hanging = vertex["hanging"];
		vp.corner = vertex["corner"];
        vp.target_pos = vertex["target_pos"].get<std::vector<double>>();
        vp.target_size = vertex["target_size"].get<std::vector<double>>();
		vp.orientation = vertex["orientation"];
		vp.size_tolerance = {};
		vp.pos_tolerance = {};
		/* Adjust orientation based on target size (NEEDED OR NOT?)
		if (!vp.target_size.empty()) {
			vp.size_tolerance = vertex["size_tolerance"].get<std::vector<double>>();
			if (vp.target_size[0] > vp.target_size[1])
				vp.orientation = dis(gen) == 0 ? FRONT : BACK;
			else
				vp.orientation = dis(gen) == 0 ? LEFT : RIGHT;
		}
		*/
		if (!vp.target_pos.empty()) {
			vp.pos_tolerance = vertex["pos_tolerance"].get<std::vector<double>>();
		}
        auto v = add_vertex(vp, inputGraph);
    }

    // Parse JSON to set edges
    for (const auto& edge : scene_graph_json["edges"]) {
        EdgeProperties ep;
        ep.type = edge["type"];
		if (ep.type == AlignWith) {
			ep.align_edge = edge["align_edge"];
			ep.distance = -1;
			ep.xyoffset = {};
		}
		else if (ep.type == CloseBy) {
			ep.xyoffset = edge["xyoffset"].get<std::vector<double>>();
			ep.distance = -1;
			ep.align_edge = -1;
		}
		else if (ep.type == Above || ep.type == Under) {
			ep.xyoffset = edge["xyoffset"].get<std::vector<double>>();
			ep.distance = edge["distance"];
			ep.align_edge = -1;
		}
		else {
			ep.distance = edge["distance"];
			ep.align_edge = -1;
			ep.xyoffset = {};
		}
        auto source = vertex(edge["source"], inputGraph);
        auto target = vertex(edge["target"], inputGraph);
        add_edge(source, target, ep, inputGraph);
    }

    g = graphProcessor.process(inputGraph, boundary, obstacles);
	if (floorplan)
    	g = graphProcessor.splitGraph2(g, boundary);

	VertexIterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
        std::cout << "Vertex " << g[*vi].id << " (" << g[*vi].label << ")" <<
            " Boundary Constraint: " << g[*vi].boundary <<
            //" Priority :" << outputGraph[*vi].priority <<
            " Orientation: " << graphProcessor.orientationnames[g[*vi].orientation] << std::endl;
    }

    EdgeIterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        std::cout << "Edge (" << g[boost::source(*ei, g)].label << " -> "
            << g[boost::target(*ei, g)].label << ") "
            << ", Type: " << graphProcessor.edgenames[g[*ei].type] << std::endl;
    }
}

void Solver::reset()
{
	inputGraph.clear();
	g.clear();
	boundary = Boundary();
	obstacles.clear();
	doors.clear();
	windows.clear();
	graphProcessor.reset();
	clearModel();
}
	

void Solver::clearModel() {
	auto vars = model.getVars();
	for (auto i = 0; i < model.get(GRB_IntAttr_NumVars); ++i) {
		model.remove(vars[i]);
	}
	auto constrs = model.getConstrs();
	for (auto i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i) {
		model.remove(constrs[i]);
	}
	auto qconstrs = model.getQConstrs();
	for (auto i = 0; i < model.get(GRB_IntAttr_NumQConstrs); ++i) {
		model.remove(qconstrs[i]);
	}
	model.update();
}

float Solver::getboundaryMaxSize()
{
	return std::max(boundary.size[0], std::max(boundary.size[1], boundary.size[2]));
}

void Solver::handleInfeasibleModel() {
	model.computeIIS();
	graphProcessor.conflict_info = "Infeasible constraints found in IIS. List of constraints: \n";
	graphProcessor.plan_info = {};
    
    GRBConstr* constrs = model.getConstrs();
    int numConstrs = model.get(GRB_IntAttr_NumConstrs);
    
    std::vector<GRBConstr> infeasibleConstraints;
    for (int i = 0; i < numConstrs; ++i) {
        if (constrs[i].get(GRB_IntAttr_IISConstr) == 1) {
            infeasibleConstraints.push_back(constrs[i]);
        }
    }
	for (int i = 0; i < infeasibleConstraints.size(); ++i) {
        std::string constrName = infeasibleConstraints[i].get(GRB_StringAttr_ConstrName);
		graphProcessor.plan_info.push_back("Constraint " + std::to_string(i) + ": " + constrName + "\n");
        //std::cout << "Constraint " << i << ": " << constrName << std::endl;
    }
	model.optimize();
}

void Solver::removeIIS(std::string name) {
	// Unused Now
	GRBConstr* constrs = model.getConstrs();
	int numConstrs = model.get(GRB_IntAttr_NumConstrs);
	size_t pos1 = name.find_first_of("0123456789");
	size_t pos2 = name.find_last_of("0123456789");
	//std::string prefix = name.substr(0, pos1) + name.substr(pos1, pos2 - pos1 + 1);
	std::string prefix = name.substr(0, pos2 + 1);
	if (name.find("Inside") != std::string::npos || name.find("NonOverlap") != std::string::npos) {
		// TODO
		std::cout << "It's not recommended to remove Inside/NonOverlap constraint." << std::endl;
	}
	else if (name.find("Pos_Tolerance") != std::string::npos) {
		model.remove(model.getConstrByName(prefix + "_x_left"));
		model.remove(model.getConstrByName(prefix + "_x_right"));
		model.remove(model.getConstrByName(prefix + "_y_front"));
		model.remove(model.getConstrByName(prefix + "_y_back"));
		if (!floorplan) {
			model.remove(model.getConstrByName(prefix + "_z_bottom"));
			model.remove(model.getConstrByName(prefix + "_z_top"));
		}
		model.update();
	}
	else if (name.find("Size_Tolerance") != std::string::npos) {
		model.remove(model.getConstrByName(prefix + "_l_min"));
		model.remove(model.getConstrByName(prefix + "_l_max"));
		model.remove(model.getConstrByName(prefix + "_w_min"));
		model.remove(model.getConstrByName(prefix + "_w_max"));
		if (!floorplan) {
			model.remove(model.getConstrByName(prefix + "_h_min"));
			model.remove(model.getConstrByName(prefix + "_h_max"));
		}
		model.update();
	}
	else if (name.find("Boundary") != std::string::npos) {
		if (name.find("Left") != std::string::npos) {
			model.remove(model.getConstrByName(prefix + "_Left_eq"));
			model.remove(model.getConstrByName(prefix + "_Left_ieq"));
			model.remove(model.getConstrByName(prefix + "_Left_ieqq"));
		}
		else if (name.find("Right") != std::string::npos) {
			model.remove(model.getConstrByName(prefix + "_Right_eq"));
			model.remove(model.getConstrByName(prefix + "_Right_ieq"));
			model.remove(model.getConstrByName(prefix + "_Right_ieqq"));
		}
		else if (name.find("Front") != std::string::npos) {
			model.remove(model.getConstrByName(prefix + "_Front_eq"));
			model.remove(model.getConstrByName(prefix + "_Front_ieq"));
			model.remove(model.getConstrByName(prefix + "_Front_ieqq"));
		}
		else if (name.find("Back") != std::string::npos) {
			model.remove(model.getConstrByName(prefix + "_Back_eq"));
			model.remove(model.getConstrByName(prefix + "_Back_ieq"));
			model.remove(model.getConstrByName(prefix + "_Back_ieqq"));
		}
		model.update();
	}
	else if (name.find("CloseBy") != std::string::npos) {
		model.remove(model.getConstrByName(prefix + "ieqa"));
		model.remove(model.getConstrByName(prefix + "ieqb"));
		model.remove(model.getConstrByName(prefix + "ieqc"));
		model.remove(model.getConstrByName(prefix + "ieqd"));
		model.remove(model.getConstrByName(prefix + "ieqe"));
		model.update();
	}
	else if (name.find("Corner") != std::string::npos) {
		model.remove(model.getConstrByName(prefix + "eqa"));
		model.remove(model.getConstrByName(prefix + "eqb"));
		model.remove(model.getConstrByName(prefix + "eqc"));
		model.update();
	}
	else {
		model.remove(model.getConstrByName(name));
		model.update();
	}
}