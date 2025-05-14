// FunctionSpace.h
#pragma once

#include <mesh.hpp>
#include <vector>
#include <array>
#include <stdexcept>

enum class FEType { P0, P1 };

class FunctionSpace {
public:
  // build either a P0 (one DOF per element) or P1 (one DOF per mesh node) space
  FunctionSpace(const Mesh& mesh, FEType feType);

  const Mesh& mesh() const { return mesh_; }

  // total number of global DOFs
  size_t numDofs() const { return globalDofs_; }

  // --- P1 only ----------------------------------
  // local vertex i (0,1,2) of element `elem` → global DOF index
  int localToGlobal(int elem, int localVertex) const;

  // --- P0 only ----------------------------------
  // (localVertex must be 0)
  int localToGlobal_P0(int elem, int localVertex) const;

private:
  const Mesh& mesh_;
  FEType       feType_;
  size_t       globalDofs_;
  std::vector<int> node2dof_;  // for P1: size = numNodes()
  std::vector<int> elem2dof_;  // for P0:  size = numElements()
};