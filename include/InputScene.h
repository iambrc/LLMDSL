/*Here we define input conditions: Boundary, Doors, Windows, Obstacles*/
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

struct Obstacles {
	std::vector<double> pos, size;
};

struct Doors {
	std::vector<double> pos, size;
	// Notice that orientation indicates the direction of doors'/windows' outward normal direction
	Orientation orientation;
};

struct Windows {
	std::vector<double> pos, size;
	Orientation orientation;
};