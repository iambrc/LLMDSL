#pragma once
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/strong_components.hpp>
#include "SceneGraph.h"
#include "InputScene.h"
#include <map>
#include <random>

struct EdgeTypeFilter {
    EdgeTypeFilter() : g(nullptr), type(LeftOf) {}
    EdgeTypeFilter(const SceneGraph& g, EdgeType type) : g(&g), type(type) {}

    bool operator()(const EdgeDescriptor& e) const {
        return (*g)[e].type == type;
    }

    const SceneGraph* g;
    EdgeType type;
};

class GraphProcessor {
public:
    GraphProcessor();
    ~GraphProcessor();

    SceneGraph process(const SceneGraph& inputGraph, const Boundary& boundary);
    void reset();

    std::vector<EdgeType> edgetypes;
    std::vector<std::string> edgenames;
    std::vector<Orientation> orientations;
    std::vector<std::string> orientationnames;
    std::string conflict_info;
    std::vector<std::string> plan_info;

private:
    void removeCycles(SceneGraph& g, EdgeType edge_type);
    Orientation oppositeOrientation(Orientation o);
    EdgeType oppositeEdgeType(EdgeType e);
};