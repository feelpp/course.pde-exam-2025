// tests.cpp
#include <mesh.hpp>
#include <basis.hpp>
#include <geomap.hpp>
#include <quadrature.hpp>
#include <functionspace.hpp>
#include <fem.hpp>
#include <error.hpp>
#include <boost/ut.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <string>

using namespace boost::ut;

int main(int argc, char* argv[]) {
    // Directory containing .msh files, passed as first argument (default "data/")
    std::string dataDir = (argc > 1 ? argv[1] : "data/");
    if (!dataDir.empty() && dataDir.back() != '/') 
        dataDir += '/';

    "Mesh read simple"_test = [&]() {
        Mesh M;
        expect(M.readGmsh(dataDir + "square.msh"));
        expect(M.numNodes() >= 4u);
    };

    "Mesh has Dirichlet boundary marker"_test = [&]() {
        Mesh M;
        expect(M.readGmsh(dataDir + "square.msh"));

        int dirTag = M.boundaryTag("Dirichlet");
        expect(dirTag != 0) << "No 'Dirichlet' boundary tag defined";

        bool foundDir = false;
        for (const auto& ln : M.lines()) {
            if (ln.physical == dirTag) {
                foundDir = true;
                break;
            }
        }
        expect(foundDir) << "No line with 'Dirichlet' boundary marker found";
    };

    "Mesh read triangle"_test = [&]() {
        Mesh M;
        expect(M.readGmsh(dataDir + "triangle2.msh"));
        expect(M.numElements() >= 1u);
    };

    "P1 basis raw"_test = [](){
        Basis B;
        double dr, ds;
        double sum = B.phi(0,0.3,0.4)+B.phi(1,0.3,0.4)+B.phi(2,0.3,0.4);
        expect(std::abs(sum-1.0)<1e-12);
        B.dphi(0,dr,ds);
        expect(dr==-1.0 && ds==-1.0);
    };

    "P1 basis eigen overloads"_test = [](){
        Basis B;
        Eigen::Vector2d p(0.2,0.3);
        double v0=B.phi(0,p), v1=B.phi(1,p), v2=B.phi(2,p);
        expect(std::abs(v0+v1+v2-1.0)<1e-12);
        double dr1=B.dphi(1,0,p), ds1=B.dphi(1,1,p);
        expect(dr1==1.0 && ds1==0.0);
    };

    "Quadrature integrates x^2 exactly on reference triangle"_test = []() {
        // Exact value: ∫_0^1 ∫_0^{1−x} x^2 dy dx = ∫_0^1 x^2(1−x) dx = 1/12
        constexpr double exact = 1.0/12.0;

        // Use degree‐2 rule to be safe
        double val = integrateTriangle(
            [](double x, double /*y*/) { return x*x; },
            2
        );
        expect(std::abs(val - exact) < 1e-8) << "Got " << val << ", expected " << exact << "\n";

        // Use degree‐4 rule to be safe
        val = integrateTriangle(
            [](double x, double /*y*/) { return x*x; },
            4
        );
        expect(std::abs(val - exact) < 1e-8) << "Got " << val << ", expected " << exact << "\n";
};
    "FunctionSpace P1 local→global"_test = [&]() {
        // single‐triangle mesh
        Mesh M; 
        expect(M.readGmsh(dataDir + "triangle2.msh"));
        FunctionSpace V1(M, FEType::P1);
        // element 0 has vertices v[0], v[1], v[2]
        const auto& tri = M.elements().at(0);
        for (int lv = 0; lv < 3; ++lv) {
            int gdof = V1.localToGlobal(0, lv);
            expect(gdof == tri.v[lv]);
        }
    };

    "FunctionSpace P0 local→global"_test = [&]() {
        Mesh M;
        expect(M.readGmsh(dataDir + "triangle2.msh"));
        FunctionSpace V0(M, FEType::P0);
        // only one DOF per element, so localVertex must be 0
        int gdof = V0.localToGlobal_P0(0, 0);
        expect(gdof == 0);
    };

    "GeometricTransform simple triangle"_test = [&]() {
        // Load single right triangle [0,0]-[1,0]-[0,1]
        Mesh M; 
        expect(M.readGmsh(dataDir + "triangle2.msh"));

        GeometricTransform G;
        // loop over all element and check that for an element the ref to real mapping gives the real points
        for (int e = 0; e < int(M.numElements()); ++e) {
            G.update(M, e);
            const auto& tri = M.elements().at(e);
            Eigen::Vector2d p0 = G.map(0.0, 0.0);
            Eigen::Vector2d p1 = G.map(1.0, 0.0);
            Eigen::Vector2d p2 = G.map(0.0, 1.0);
            expect(p0 == M.nodes().at(tri.v[0]));
            expect(p1 == M.nodes().at(tri.v[1]));
            expect(p2 == M.nodes().at(tri.v[2]));
        }
    };
#if 0    
    "GeometricTransform triangle area"_test = [&]() {
        for( auto & s : std::map<std::string,double>{{"triangle2.msh", 0.5}, {"square.msh", 1.}} ) {
            Mesh M; 
            expect(M.readGmsh(dataDir + s.first));
            GeometricTransform G;
            // loop over all element and check the area of the domain (a triangle) is s.second
            double area = 0.0;
            for (int e = 0; e < int(M.numElements()); ++e) {
                G.update(M, e);
                area += G.detJ()/2;
            }
            // area of a triangle is s.second
            expect(std::abs(area - s.second) < 1e-12) << "Area mismatch: " << area << " != " << s.second;
        }
    };
#endif
    // Dirichlet BC exact solutions for various functions
    "Dirichlet BC exact solutions"_test = [&]() {

        for( auto & s : std::map<std::string,double>{{"triangle2.msh", 0.5}, {"square.msh", 1.}} ) {
            std::cout << "testing  messh " << s.first << "\n";
            Mesh M; expect(M.readGmsh(dataDir + s.first));
            FunctionSpace V(M, FEType::P1);

            // Ensure we have a “Dirichlet” boundary definition
            int dirTag = M.boundaryTag("Dirichlet");
            expect(dirTag != 0) << "Mesh missing 'Dirichlet' boundary tag";

            // Zero right-hand side
            auto zero_f = [](double,double){ return 0.0; };
            LinearForm  L(V, zero_f);
            BilinearForm A(V);

            // Candidate exact solutions (all harmonic): constant, x, y, x+2y
            std::vector<std::pair<std::string, std::function<double(double,double)>>> funcs = {
                {"constant", [](double,double){           return 5.0; }},
                {"x",        [](double x,double){         return x;   }},
                {"y",        [](double,double y){         return y;   }},
                {"x+2y",     [](double x,double y){       return x + 2.0*y; }}
            };

            for (auto& [name, g] : funcs) {
                std::cout << "testing " << name << "\n";
                auto ue = V.interpolate(g);
                std::map<std::string,std::function<double(double,double)>> bc;
                bc["Dirichlet"] = g;

                auto sol = A.solve(L, bc);
                //std::cout << "ue=" << ue.transpose() << "\n";
                for (int i = 0; i < int(sol.size()); ++i) {
                    auto p = M.nodes()[i];
                    double uex = g(p.x(), p.y());
                    expect(std::abs(sol(i) - uex) < 1e-12)
                        << "Function '" << name << "' mismatch at node " << i;
                }
                
        
                // Check L2 error
                double err = computeL2Error(V, sol, g);
                expect(err < 1e-12) << name << " L2 error: " << err;
            }
        }
    };



    #if 0
    "L2 error zero"_test = [&](){
        Mesh M; expect(M.readGmsh(dataDir+"triangle2.msh"));
        FunctionSpace V(M);
        Eigen::VectorXd sol = Eigen::VectorXd::Zero(V.numDofs());
        auto exact = [](double,double){return 0.0;};
        double err = computeL2Error(V,sol,exact);
        expect(err==0.0);
    };
#endif    

    return 0;
}