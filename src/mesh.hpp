// Mesh.h
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <map>

//===========================================================
// zero‐based triangle connectivity
struct Triangle {
    int v[3];
    int physical;    // physical region tag
};

// zero‐based line (edge) connectivity for boundary markers
struct Line {
    int v[2];
    int physical;    // physical boundary tag
};

class Mesh {
public:
    Mesh() = default;

    // Read Gmsh v2 ASCII file. Returns true on success.
    bool readGmsh(const std::string& filename);

    size_t numNodes()    const { return nodes_.size(); }
    size_t numElements() const { return elems_.size(); }
    size_t numLines()    const { return lines_.size(); }

    // Access to node coordinates and element connectivity
    const std::vector<Eigen::Vector2d>& nodes()    const { return nodes_; }
    const std::vector<Triangle>&         elements() const { return elems_; }
    const std::vector<Line>&             lines()    const { return lines_; }

    // Markers
    // physical tag of node (from boundary lines), 0 if none
    int nodeMarker(size_t node)     const { return nodeMarkers_.at(node); }
    // physical tag of element (volume),  0 if none
    int elementMarker(size_t elem)  const { return elems_.at(elem).physical; }

    // Mapping from human-readable names to tags (populated from $PhysicalNames)
    int boundaryTag(const std::string& name) const {
        auto it = boundaryTags_.find(name);
        return (it != boundaryTags_.end() ? it->second : 0);
    }
    int regionTag(const std::string& name) const {
        auto it = regionTags_.find(name);
        return (it != regionTags_.end() ? it->second : 0);
    }

private:
    std::vector<Eigen::Vector2d> nodes_;
    std::vector<int>             nodeMarkers_;
    std::vector<Triangle>        elems_;
    std::vector<Line>            lines_;

    // user-supplied name→tag maps
    std::map<std::string,int>    regionTags_;
    std::map<std::string,int>    boundaryTags_;
};