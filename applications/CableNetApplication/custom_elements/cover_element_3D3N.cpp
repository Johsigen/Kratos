//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Johanna Sigeneger
//
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/cover_element_3D3N.hpp"
#include "includes/checks.h"
#include "includes/define.h"


namespace Kratos {
CoverElement3D3N::CoverElement3D3N(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

CoverElement3D3N::CoverElement3D3N(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
CoverElement3D3N::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_intrusive<CoverElement3D3N>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
CoverElement3D3N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_intrusive<CoverElement3D3N>(NewId, pGeom,
                                               pProperties);
}

CoverElement3D3N::~CoverElement3D3N() {}

void CoverElement3D3N::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                         VectorType &rRightHandSideVector,
                                         const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    //GetNormal();
    //KRATOS_WATCH(this->Id());
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
  rRightHandSideVector = ZeroVector(local_size);
  KRATOS_CATCH("")
}

array_1d<double,3> CoverElement3D3N::GetNormal() const
{
    array_1d<double,3> center {0.0, 0.0, 0.0};
    center = this->GetGeometry().Center();
    KRATOS_WATCH(this->GetGeometry().Normal(center))
    return this->GetGeometry().Normal(center);
}
} // namespace Kratos.