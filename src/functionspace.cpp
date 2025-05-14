// FunctionSpace.cpp
#include <functionspace.hpp>

FunctionSpace::FunctionSpace(const Mesh& mesh, FEType feType)
  : mesh_(mesh), feType_(feType)
{
  if (feType_ == FEType::P1) {
    // one DOF per mesh node
    globalDofs_ = mesh_.numNodes();
    node2dof_.resize(globalDofs_);
    for (int i = 0; i < int(globalDofs_); ++i) {
      node2dof_[i] = i;
    }
  }
  else {
    // one DOF per element
    globalDofs_ = mesh_.numElements();
    elem2dof_.resize(globalDofs_);
    for (int e = 0; e < int(globalDofs_); ++e) {
      elem2dof_[e] = e;
    }
  }
}

int FunctionSpace::localToGlobal(int elem, int localVertex) const
{
  if (feType_ != FEType::P1)
    throw std::logic_error("localToGlobal(elem,i) only valid for P1");
  if (elem < 0 || elem >= int(mesh_.numElements()))
    throw std::out_of_range("element index out of range");
  if (localVertex < 0 || localVertex > 2)
    throw std::out_of_range("local vertex must be 0, 1, or 2");

  const auto& tri = mesh_.elements().at(elem);
  int nodeId = tri.v[localVertex];         // mesh node index
  return node2dof_.at(nodeId);             // global DOF
}

int FunctionSpace::localToGlobal_P0(int elem, int localVertex) const
{
  if (feType_ != FEType::P0)
    throw std::logic_error("localToGlobal_P0 only valid for P0");
  if (localVertex != 0)
    throw std::out_of_range("P0 has only one local DOF (vertex 0)");
  if (elem < 0 || elem >= int(elem2dof_.size()))
    throw std::out_of_range("element index out of range");

  return elem2dof_[elem];  // one DOF per element
}