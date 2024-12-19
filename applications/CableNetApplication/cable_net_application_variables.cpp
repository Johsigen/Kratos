//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#include "cable_net_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(Vector, SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL)
    KRATOS_CREATE_VARIABLE(double, NORMALFORCE)
    KRATOS_CREATE_VARIABLE(double, BENDING_STIFFNESS)
    KRATOS_CREATE_VARIABLE(double, MINIMAL_LENGTH)
    KRATOS_CREATE_VARIABLE(double, MINIMAL_CIRCUMFERENCE)
}
