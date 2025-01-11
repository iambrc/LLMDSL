#pragma once

#include "GraphProcessor.h"
#include <boost/graph/graphviz.hpp>
#include <fstream>
#include <gurobi_c++.h>
#include <nlohmann/json.hpp>
#include <polyclipping/clipper.hpp>

extern std::vector<std::string> show_edges;
extern std::vector<std::string> show_orientations;

template <class Vertex>
struct vertex_writer_out {
    vertex_writer_out(const SceneGraph& g) : g(g) {}
    void operator()(std::ostream& out, const Vertex& v) const {
        const VertexProperties& vp = g[v];
		double l = !vp.size.empty() ? vp.size[0] : 0.1;
		double w = !vp.size.empty() ? vp.size[1] : 0.1;
        double show_l = !vp.size.empty() ? vp.size[0] : -1;
        double show_w = !vp.size.empty() ? vp.size[1] : -1;
        out << "[label=\"" << vp.label << "\\n"
            << "Size: " << show_l << " " << show_w << "\\n"
            //<< "Priority: " << vp.priority << "\\n"
            << "Orientation: " << show_orientations[vp.orientation] << "\", "
            << "shape=rect, style=filled, fillcolor=\"lightblue\", "
            << "width=" << l << ", height=" << w << "]";
    }
    const SceneGraph& g;
};
template <class Vertex>
struct vertex_writer_in {
    vertex_writer_in(const SceneGraph& g) : g(g) {}
    void operator()(std::ostream& out, const Vertex& v) const {
        const VertexProperties& vp = g[v];
        double l = !vp.target_size.empty() ? vp.target_size[0] : 0.1;
        double w = !vp.target_size.empty() ? vp.target_size[1] : 0.1;
        double show_l = !vp.target_size.empty() ? vp.target_size[0] : -1;
        double show_w = !vp.target_size.empty() ? vp.target_size[1] : -1;
        out << "[label=\"" << vp.label << "\\n"
            << "Size: " << show_l << " " << show_w << "\\n"
            //<< "Priority: " << vp.priority << "\\n"
            << "Orientation: " << show_orientations[vp.orientation] << "\", "
            << "shape=rect, style=filled, fillcolor=\"lightblue\", "
            << "width=" << l << ", height=" << w << "]";
    }
    const SceneGraph& g;
};

template <class Edge>
struct edge_writer {
    edge_writer(const SceneGraph& g) : g(g) {}
    void operator()(std::ostream& out, const Edge& e) const {
        const EdgeProperties& ep = g[e];
        std::string color;
        switch (ep.type) {
        case LeftOf: color = "red"; break;
        case RightOf: color = "blue"; break;
        case FrontOf: color = "green"; break;
        case Behind: color = "orange"; break;
        case Above: color = "purple"; break;
        case Under: color = "brown"; break;
        case CloseBy: color = "cyan"; break;
        case AlignWith: color = "magenta"; break;
        }
        out << "[label=\"";
        out << show_edges[ep.type] << "\\n";

        if (ep.type == AlignWith) {
            out << "Align edge: " << ep.align_edge << "\\n";
        }
        else if (ep.type != CloseBy && ep.distance >= 0) {
            out << "Distance: " << ep.distance << "\\n";
        }
        out << "\", color=\"" << color << "\"]";
    }
    const SceneGraph& g;
};

class Solver {
public:
    Solver();
    ~Solver();

    void solve();
    void saveGraph();
    void readSceneGraph(const std::string& path, float wallwidth);
    void reset();
    SceneGraph getsolution() { return g; }
    float getboundaryMaxSize();
    Boundary getboundary() { return boundary; }

    bool floorplan;
    std::vector<double> hyperparameters;
    int scalingFactor;
private:
    bool has_path(const SceneGraph& g, VertexDescriptor start, VertexDescriptor target);
    bool dfs_check_path(const SceneGraph& g, VertexDescriptor u, VertexDescriptor target, EdgeType required_type, std::vector<bool>& visited);
    void addConstraints();
    void optimizeModel();
    void handleInfeasibleModel();
    void removeIIS(std::string name);
    void clearModel();

    SceneGraph inputGraph, g;
    Boundary boundary;
    std::vector<Obstacles> obstacles;
    std::vector<Doors> doors;
    std::vector<Windows> windows;
    GraphProcessor graphProcessor;

    GRBEnv env;
    GRBModel model;

    std::string inputpath;
};