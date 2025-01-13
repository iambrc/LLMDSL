/*Here we define input conditions: Boundary*/
#pragma once
#include <vector>
#include "SceneGraph.h"

struct Boundary {
	std::vector<std::vector<double>> points;
	std::vector<double> origin_pos, size;
	// Notice that orientation indicates the direction in which the outward normal direction of wall segment
	std::vector<Orientation> Orientations;
	std::vector<int> BLcorner, BRcorner, TLcorner, TRcorner;
};