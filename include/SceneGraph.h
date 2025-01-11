/*Here we define the SceneGraph data structure, which is a directed graph that represents the scene.*/
#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>


enum Orientation { UP, DOWN, LEFT, RIGHT, FRONT, BACK };
// Notice that Frontof means coordinates y is increasing, and Behind means coordinates y is decreasing.
enum EdgeType { LeftOf, RightOf, FrontOf, Behind, Above, Under, CloseBy, AlignWith};

enum CornerType { TOPLEFT, TOPRIGHT, BOTTOMLEFT, BOTTOMRIGHT };

struct VertexProperties {
	std::string label;
	int id, boundary;
	CornerType corner;
	std::vector<double> target_pos, target_size, pos, size, pos_tolerance, size_tolerance;
	Orientation orientation;
	bool on_floor, hanging;
};

struct EdgeProperties {
	// Notice that distance means the distance between two edges, not the center of rectangles.
	double distance;
	std::vector<double> xyoffset;
	// Notice that align_edge = {0, 1, 2, 3}, each number represents the alignment of bottom/right/up/left ,respectively.
	int align_edge;
	EdgeType type;
};

typedef boost::adjacency_list<
	boost::vecS,
	boost::vecS,
	boost::directedS,
	VertexProperties,
	EdgeProperties
> SceneGraph;

typedef boost::graph_traits<SceneGraph>::edge_descriptor EdgeDescriptor;
typedef boost::graph_traits<SceneGraph>::vertex_descriptor VertexDescriptor;
typedef boost::graph_traits<SceneGraph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<SceneGraph>::vertex_iterator VertexIterator;

