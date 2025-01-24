//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors: Johanna Sigeneger
//
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{
/**
 * @class CoverElement3D3N
 *
 * @brief This is an element to cover the surface of a ring net
 *
 * @author Johanna Sigeneger
 */

class KRATOS_API(CABLE_NET_APPLICATION) CoverElement3D3N : public Element
{
protected:

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CoverElement3D3N);


    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::DofsVectorType DofsVectorType;


    CoverElement3D3N() {};
    CoverElement3D3N(IndexType NewId,
                    GeometryType::Pointer pGeometry);
    CoverElement3D3N(IndexType NewId,
                    GeometryType::Pointer pGeometry,
                    PropertiesType::Pointer pProperties);

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    void EquationIdVector(EquationIdVectorType &rResult,
                        const ProcessInfo &rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;
            
    void GetValuesVector(Vector& rValues,int Step = 0) const override;

    void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

private:

    array_1d<double,3> GetNormal() const;
};


} // namespace Kratos

