// Mesh.cpp
#include <mesh.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

bool Mesh::readGmsh(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return false;
    }
    std::cout << "Reading Gmsh file: " << filename << "\n";

    std::string line;
    // 1) Look for $PhysicalNames and parse name→tag
    while (std::getline(in, line)) {
        if (line == "$PhysicalNames") {
            int nNames;
            in >> nNames;
            std::getline(in, line); // consume end of that line
            for (int i = 0; i < nNames; ++i) {
                std::getline(in, line);
                std::istringstream iss(line);
                int dim, tag;
                iss >> dim >> tag;
                iss >> std::ws;
                std::string name;
                if (iss.peek() == '\"') {
                    iss.get();                     // skip opening "
                    std::getline(iss, name, '\"'); // read until closing "
                } else {
                    iss >> name;
                }
                if (dim == 1) {
                    boundaryTags_[name] = tag;
                } else if (dim == 2) {
                    regionTags_[name] = tag;
                }
            }
            // skip the $EndPhysicalNames line
            std::getline(in, line);
            break;
        }
    }

    // rewind to start so we parse nodes & elements as before
    in.clear();
    in.seekg(0, std::ios::beg);

    // 2) Read $Nodes
    while (std::getline(in, line)) {
        if (line == "$Nodes") {
            std::size_t N; in >> N;
            nodes_.resize(N);
            nodeMarkers_.assign(N, 0);
            std::getline(in, line); // end of line
            for (std::size_t i = 0; i < N; ++i) {
                int id;
                double x,y,z;
                in >> id >> x >> y >> z;
                nodes_[id-1] = Eigen::Vector2d(x,y);
            }
            std::getline(in, line); // $EndNodes
        }
        else if (line == "$Elements") {
            std::size_t M; in >> M;
            elems_.clear();
            lines_.clear();
            std::getline(in, line); // end of line
            for (std::size_t e = 0; e < M; ++e) {
                std::getline(in, line);
                std::istringstream iss(line);
                int id, type, numTags;
                iss >> id >> type >> numTags;
                std::vector<int> tags(numTags);
                for (int t = 0; t < numTags; ++t) iss >> tags[t];
                int phys = tags.empty() ? 0 : tags[0];

                if (type == 2) {
                    // triangle
                    int n1,n2,n3;
                    iss >> n1 >> n2 >> n3;
                    Triangle tri;
                    tri.v[0] = n1-1; tri.v[1] = n2-1; tri.v[2] = n3-1;
                    tri.physical = phys;
                    elems_.push_back(tri);
                }
                else if (type == 1) {
                    // line
                    int n1,n2;
                    iss >> n1 >> n2;
                    Line ln;
                    ln.v[0] = n1-1; ln.v[1] = n2-1;
                    ln.physical = phys;
                    lines_.push_back(ln);
                    // tag the two nodes for boundary markers
                    nodeMarkers_[n1-1] = phys;
                    nodeMarkers_[n2-1] = phys;
                }
            }
        }
    }

    std::cout << "Read " << nodes_.size()
              << " nodes, " << elems_.size()
              << " triangles, " << lines_.size()
              << " boundary lines.\n";

    return !nodes_.empty() && !elems_.empty();
}