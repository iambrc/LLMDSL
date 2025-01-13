#include "GraphProcessor.h"
#include <iostream>

GraphProcessor::GraphProcessor() {
    edgetypes = { LeftOf, RightOf, FrontOf, Behind, Above, Under, CloseBy, AlignWith };
    edgenames = { "LeftOf", "RightOf", "FrontOf", "Behind", "Above", "Under", "CloseBy", "AlignWith" };
    orientations = { UP, DOWN, LEFT, RIGHT, FRONT, BACK };
    orientationnames = { "UP", "DOWN", "LEFT", "RIGHT", "FRONT", "BACK" };
}

GraphProcessor::~GraphProcessor() {}

Orientation GraphProcessor::oppositeOrientation(Orientation o) {
    switch (o) {
    case UP: return DOWN;
    case DOWN: return UP;
    case LEFT: return RIGHT;
    case RIGHT: return LEFT;
    case FRONT: return BACK;
    case BACK: return FRONT;
    }
    return o;
}

EdgeType GraphProcessor::oppositeEdgeType(EdgeType e) {
    switch (e) {
    case LeftOf: return RightOf;
    case RightOf: return LeftOf;
    case FrontOf: return Behind;
    case Behind: return FrontOf;
    case Above: return Under;
    case Under: return Above;
    }
    return e;
}

void GraphProcessor::removeCycles(SceneGraph& g, EdgeType edge_type) {
    EdgeTypeFilter edge_filter(g, edge_type);
    boost::filtered_graph<SceneGraph, EdgeTypeFilter> filtered_g(g, edge_filter);

    // Find the strongly connected component
    std::vector<int> component(num_vertices(filtered_g));
    int num = boost::strong_components(filtered_g, &component[0]);

    // Find rings
    std::vector<std::vector<EdgeDescriptor>> cycles_edges;
    for (int i = 0; i < num; ++i) {
        std::vector<EdgeDescriptor> edges_to_remove;
        for (const auto &e : boost::make_iterator_range(edges(filtered_g))) {
            if (component[source(e, filtered_g)] == i && component[target(e, filtered_g)] == i) {
                edges_to_remove.push_back(e);
            }
        }
        if (!edges_to_remove.empty()) {
            cycles_edges.push_back(edges_to_remove);
        }
    }

    // Report cycles
    if (!cycles_edges.empty()) {
        conflict_info = "Cycles found, please choose to remove conflict constraints: \n";
        plan_info = {};
        
        for (const auto &edges : cycles_edges) {
            for (auto i = 0; i < edges.size(); ++i) {
                plan_info.push_back("Plan " + std::to_string(i) + ": Remove edge " + g[source(edges[i], g)].label + " " + edgenames[g[edges[i]].type] + " " + g[target(edges[i], g)].label + "\n");
            }
        }
    }
    
}

SceneGraph GraphProcessor::process(const SceneGraph& inputGraph, const Boundary& boundary)
{
    SceneGraph outputGraph = inputGraph;
    // Find and work with rings in each type of edge
    std::vector<std::pair<VertexDescriptor, VertexDescriptor>> edges_to_reverse, edges_to_remove;
	std::vector<EdgeProperties> new_edge_properties, removed_edge_properties;
    EdgeIterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(outputGraph); ei != ei_end; ++ei) {
		VertexDescriptor v1 = boost::source(*ei, outputGraph), v2 = boost::target(*ei, outputGraph);
        EdgeProperties ep = outputGraph[*ei];
        if (ep.type == RightOf) {
            ep.type = LeftOf;
        }
        else if (ep.type == Behind) {
            ep.type = FrontOf;
        }
        else if (ep.type == Under) {
            ep.type = Above;
        }
        else {
            continue;
        }
        edges_to_reverse.push_back(std::make_pair(v1, v2));
		new_edge_properties.push_back(ep);
    }
    for (size_t i = 0; i < edges_to_reverse.size(); ++i) {
        boost::graph_traits<SceneGraph>::out_edge_iterator oe_i, oe_end;
        for (boost::tie(oe_i, oe_end) = boost::out_edges(edges_to_reverse[i].first, outputGraph); oe_i != oe_end; ++oe_i) {
            if (boost::target(*oe_i, outputGraph) == edges_to_reverse[i].second && outputGraph[*oe_i].type == oppositeEdgeType(new_edge_properties[i].type)) { 
                boost::remove_edge(*oe_i, outputGraph);
                break;
            }
        }
        //boost::remove_edge(edges_to_reverse[i].first, edges_to_reverse[i].second, outputGraph);
        boost::add_edge(edges_to_reverse[i].second, edges_to_reverse[i].first, new_edge_properties[i], outputGraph);
    }
    for (EdgeType edgetype : {LeftOf, FrontOf, Above}) {
        removeCycles(outputGraph, edgetype);
    }
    // Find the contradiction between boundary constraints and position/orientation constraints
    VertexIterator vi, vi_end;
    std::vector<VertexDescriptor> vertices_to_remove;
    for (boost::tie(vi, vi_end) = boost::vertices(outputGraph); vi != vi_end; ++vi) {
        Orientation o_vi = outputGraph[*vi].orientation;
        if (outputGraph[*vi].boundary >= 0 && boundary.Orientations[outputGraph[*vi].boundary] == o_vi) {
            conflict_info = "Conflict found: Object " + outputGraph[*vi].label + " face the wall, please select a plan: \n";
            plan_info = {};
            plan_info.push_back("Plan 0: Adjust orientation\n");
            plan_info.push_back("Plan 1: Remove boundary constraint\n");
            plan_info.push_back("Plan 2: Remove object: " + outputGraph[*vi].label + "\n");
        }
        // check for contradictions between position/size constraints and on-floor constarints
        if (!outputGraph[*vi].target_pos.empty() && !outputGraph[*vi].pos_tolerance.empty() &&
            !outputGraph[*vi].target_size.empty() && !outputGraph[*vi].size_tolerance.empty() &&
            outputGraph[*vi].target_pos[2] - outputGraph[*vi].pos_tolerance[2] >
            outputGraph[*vi].target_size[2] / 2 + outputGraph[*vi].size_tolerance[2] / 2 && outputGraph[*vi].on_floor)
            {
                conflict_info = "Contradiction found between position/size constraints and on-floor constraints of Object "+ outputGraph[*vi].label + ", please select a plan: \n";
                plan_info = {};
                plan_info.push_back("Plan 0: Remove on-floor constraint\n");
                plan_info.push_back("Plan 1: Adjust position/size constraints\n");
                plan_info.push_back("Plan 2: Remove object: " + outputGraph[*vi].label + "\n");
            }
    }
    return outputGraph;
}

void GraphProcessor::reset()
{
    conflict_info = "";
    plan_info = {};
}