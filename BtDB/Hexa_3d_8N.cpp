// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_elements/base_solid_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

void BaseSolidElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if( GetProperties().Has(INTEGRATION_ORDER) ) {
            const SizeType integration_order = GetProperties()[INTEGRATION_ORDER];
            switch ( integration_order )
            {
            case 1:
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                break;
            case 2:
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
                break;
            case 3:
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
                break;
            case 4:
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
                break;
            case 5:
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
                break;
            default:
                KRATOS_WARNING("BaseSolidElement") << "Integration order " << integration_order << " is not available, using default integration order for the geometry" << std::endl;
                mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
            }
        } else {
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
        }

        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
            mConstitutiveLawVector.resize( integration_points.size() );

        InitializeMaterial();

    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    // We initialize the material reponse if required
    bool required = false;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                // Call the constitutive law to update material variables
                mConstitutiveLawVector[point_number]->InitializeMaterialResponse(Values, GetStressMeasure());

                // TODO: Deprecated, remove this
                mConstitutiveLawVector[point_number]->InitializeSolutionStep( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->InitializeNonLinearIteration( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->FinalizeNonLinearIteration( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    // We finalize the material reponse if required
    bool required = false;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                // Call the constitutive law to update material variables
                mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, GetStressMeasure());

                // TODO: Deprecated, remove this
                mConstitutiveLawVector[point_number]->FinalizeSolutionStep( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeMaterial()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure BaseSolidElement::GetStressMeasure() const
{
    return ConstitutiveLaw::StressMeasure_PK2;
}

/***********************************************************************************/
/***********************************************************************************/

bool BaseSolidElement::UseElementProvidedStrain() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number )
            mConstitutiveLawVector[point_number]->ResetMaterial( r_properties,  r_geometry, row( N_values, point_number ) );
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer BaseSolidElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    KRATOS_WARNING("BaseSolidElement") << " Call BaseSolidElement (base class) Clone " << std::endl;

    BaseSolidElement::Pointer p_new_elem = Kratos::make_intrusive<BaseSolidElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = velocity[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = acceleration[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {
        VectorType element_mass_vector(mat_size);
        this->CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * dimension;

            AtomicAdd(r_geom[i].GetValue(NODAL_MASS), element_mass_vector[index]);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType element_size = dimension * number_of_nodes;

    Vector damping_residual_contribution = ZeroVector(element_size);

    // Calculate damping contribution to residual -->
    if (r_prop.Has(RAYLEIGH_ALPHA) || r_prop.Has(RAYLEIGH_BETA)) {
        Vector current_nodal_velocities = ZeroVector(element_size);
        this->GetFirstDerivativesVector(current_nodal_velocities);

        Matrix damping_matrix(element_size, element_size);
        this->CalculateDampingMatrixWithLumpedMass(damping_matrix, rCurrentProcessInfo);

        // Current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = dimension * i;

            array_1d<double, 3>& r_force_residual = r_geom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < dimension; ++j) {
                AtomicAdd(r_force_residual[j], (rRHSVector[index + j] - damping_residual_contribution[index + j]));
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    //calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                             const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = false;
    VectorType RHS;

    CalculateAll( rLeftHandSideMatrix, RHS, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // Checking if computing lumped mass matrix
    const bool compute_lumped_mass_matrix = StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(r_prop, rCurrentProcessInfo);

    // LUMPED MASS MATRIX
    if (compute_lumped_mass_matrix) {
        VectorType temp_vector(mat_size);
        this->CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);
        for (IndexType i = 0; i < mat_size; ++i)
            rMassMatrix(i, i) = temp_vector[i];
    } else { // CONSISTENT MASS
        const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
        const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

        Matrix J0(dimension, dimension);

        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            GeometryUtils::JacobianOnInitialConfiguration(
                r_geom, integration_points[point_number], J0);
            const double detJ0 = MathUtils<double>::Det(J0);
            const double integration_weight =
                GetIntegrationWeight(integration_points, point_number, detJ0) * thickness;
            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    const SizeType index_j = j * dimension;
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                    for ( IndexType k = 0; k < dimension; ++k )
                        rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const unsigned int mat_size = GetGeometry().PointsNumber() * GetGeometry().WorkingSpaceDimension();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        mat_size);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number <number_of_integration_points; ++point_number ) {
            bool value;
            mConstitutiveLawVector[point_number]->GetValue( rVariable, value);
            rOutput[point_number] = value;
        }
    } else {
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ++ii ) {
            bool solution;
            solution = mConstitutiveLawVector[ii]->CalculateValue( Values, rVariable, solution);
            rOutput[ii] = solution;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const std::size_t number_of_integration_points = integration_points.size();
    const auto& r_geometry = GetGeometry();

    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_WEIGHT) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                this_kinematic_variables.detJ0 = CalculateDerivativesOnReferenceConfiguration(this_kinematic_variables.J0,
                                                                                    this_kinematic_variables.InvJ0,
                                                                                    this_kinematic_variables.DN_DX,
                                                                                    point_number,
                                                                                    this->GetIntegrationMethod());

                double integration_weight = GetIntegrationWeight(integration_points,
                                                                    point_number,
                                                                    this_kinematic_variables.detJ0);

                if (dimension == 2 && this->GetProperties().Has(THICKNESS))
                    integration_weight *= this->GetProperties()[THICKNESS];

                rOutput[point_number] = integration_weight;
            }
        } else if ( rVariable == STRAIN_ENERGY ) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            // Reading integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                double StrainEnergy = 0.0;

                mConstitutiveLawVector[point_number]->CalculateValue(Values, STRAIN_ENERGY, StrainEnergy);

                rOutput[point_number] = StrainEnergy;
            }
        } else if ( rVariable == ERROR_INTEGRATION_POINT ) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            // Reading integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(  );

            //Calculate Cauchy Stresses from the FE solution
            std::vector<Vector> sigma_FE_solution(number_of_nodes);
            const Variable<Vector>& r_variable_stress = CAUCHY_STRESS_VECTOR;
            CalculateOnIntegrationPoints(r_variable_stress, sigma_FE_solution, rCurrentProcessInfo);

            // calculate the determinatn of the Jacobian in the current configuration
            Vector detJ(number_of_integration_points);
            detJ = r_geometry.DeterminantOfJacobian(detJ);

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            if (r_geometry[0].Has(RECOVERED_STRESS)) {
                for (IndexType point_number = 0; point_number < number_of_integration_points; point_number++) {
                    // Compute element kinematics B, F, DN_DX ...
                    CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                    double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ[point_number]);

                    if (dimension == 2 && this->GetProperties().Has(THICKNESS))
                        integration_weight *= this->GetProperties()[THICKNESS];

                    // Calculate recovered stresses at integration points
                    Vector sigma_recovered = ZeroVector(strain_size);

                    // sigma_recovered = sum(N_i * sigma_recovered_i)
                    for (IndexType node_number=0; node_number<number_of_nodes; node_number++) {
                        const auto& r_sigma_recovered_node = r_geometry[node_number].GetValue(RECOVERED_STRESS);
                        for (IndexType stress_component = 0; stress_component<strain_size; stress_component++) {
                            sigma_recovered[stress_component] += this_kinematic_variables.N[node_number] * r_sigma_recovered_node[stress_component];
                        }
                    }

                    // Calculate error_sigma
                    Vector error_sigma(strain_size);
                    error_sigma = sigma_recovered - sigma_FE_solution[point_number];

                    // For debug
                    KRATOS_TRACE("ERROR_INTEGRATION_POINT")
                    <<"sigma recovered: " << sigma_recovered << std::endl
                    <<"sigma FE: " << sigma_FE_solution[point_number] << std::endl;

                    // Calculate inverse of material matrix
                    Matrix invD(strain_size,strain_size);
                    double detD;
                    MathUtils<double>::InvertMatrix(this_constitutive_variables.D, invD,detD);

                    // Calculate error_energy
                    rOutput[point_number] = integration_weight * inner_prod(error_sigma, prod(invD, error_sigma));
                }
            } else {
                for (IndexType point_number = 0; point_number < number_of_integration_points; point_number++) {
                    rOutput[point_number] = 0.0;
                }
            }
        } else if (rVariable == VON_MISES_STRESS) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                // Compute VM stress
                if (dimension == 2 ) {
                    rOutput[point_number] = ConstitutiveLawUtilities<3>::CalculateVonMisesEquivalentStress(this_constitutive_variables.StressVector);
                } else {
                    rOutput[point_number] = ConstitutiveLawUtilities<6>::CalculateVonMisesEquivalentStress(this_constitutive_variables.StressVector);
                }
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_COORDINATES) {
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                Point global_point;
                GetGeometry().GlobalCoordinates(global_point, integration_points[point_number]);

                rOutput[point_number] = global_point.Coordinates();
            }
        } else if (rVariable == LOCAL_AXIS_1 || rVariable == LOCAL_AXIS_2 || rVariable == LOCAL_AXIS_3) {
            if (this->Has(rVariable)) {
                for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
                    rOutput[point_number] = this->GetValue(rVariable);
            } else if (rVariable == LOCAL_AXIS_3) {
                const array_1d<double, 3> r_local_axis_1 = this->GetValue(LOCAL_AXIS_1);
                const array_1d<double, 3> local_axis_2 = this->GetValue(LOCAL_AXIS_2);
                for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
                    rOutput[point_number] = MathUtils<double>::CrossProduct(r_local_axis_1, local_axis_2);
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    }  else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == INSITU_STRESS ) {
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
            Vector strain_vector( strain_size );

            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size() != strain_vector.size() )
                    rOutput[point_number].resize( strain_vector.size(), false );

                rOutput[point_number] = mConstitutiveLawVector[point_number]->GetValue( INSITU_STRESS, rOutput[point_number] );
            }
        } else if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                //call the constitutive law to update material variables
                if( rVariable == CAUCHY_STRESS_VECTOR) {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
                } else {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points,ConstitutiveLaw::StressMeasure_PK2);
                }

                if ( rOutput[point_number].size() != strain_size )
                    rOutput[point_number].resize( strain_size, false );

                rOutput[point_number] = this_constitutive_variables.StressVector;
            }
        } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags &ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            const ConstitutiveLaw::StressMeasure this_stress_measure = rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ? ConstitutiveLaw::StressMeasure_PK2 : ConstitutiveLaw::StressMeasure_Kirchhoff;

            //reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this_stress_measure);

                if ( rOutput[point_number].size() != strain_size)
                    rOutput[point_number].resize( strain_size, false );

                rOutput[point_number] = this_constitutive_variables.StrainVector;
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
            std::vector<Vector> stress_vector;

            if( rVariable == CAUCHY_STRESS_TENSOR )
                this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
            else
                this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != dimension )
                    rOutput[point_number].resize( dimension, dimension, false );

                rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
            }
        }
        else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR) {
            std::vector<Vector> strain_vector;
            if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
                CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );
            else
                CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != dimension )
                    rOutput[point_number].resize( dimension, dimension, false );

                rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(strain_vector[point_number]);
            }
        } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);
            Values.SetConstitutiveMatrix(this_constitutive_variables.D); //this is the output parameter

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                if( rOutput[point_number].size2() != this_constitutive_variables.D.size2() )
                    rOutput[point_number].resize( this_constitutive_variables.D.size1() , this_constitutive_variables.D.size2() , false );

                rOutput[point_number] = this_constitutive_variables.D;
            }
        } else if ( rVariable == DEFORMATION_GRADIENT ) { // VARIABLE SET FOR TRANSFER PURPOUSES
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                if( rOutput[point_number].size2() != this_kinematic_variables.F.size2() )
                    rOutput[point_number].resize( this_kinematic_variables.F.size1() , this_kinematic_variables.F.size2() , false );

                rOutput[point_number] = this_kinematic_variables.F;
            }
        }  else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<bool>& rVariable,
    const std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<int>& rVariable,
    const std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    const std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    const std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        for ( IndexType point_number = 0; point_number < integration_points_number; ++point_number ) {
            mConstitutiveLawVector[point_number] = rValues[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    const std::vector<array_1d<double, 6>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    const std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

int  BaseSolidElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    // Basic check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You have called to the CalculateAll from the base class for solid elements" << std::endl;
}

//***********************************************************************
//***********************************************************************

double BaseSolidElement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const
{
    return rThisIntegrationPoints[PointNumber].Weight() * detJ;
}

//***********************************************************************
//***********************************************************************

void BaseSolidElement::CalculateShapeGradientOfMassMatrix(MatrixType& rMassMatrix, ShapeParameter Deriv) const
{
    KRATOS_TRY;

    // Properties
    const auto& r_prop = GetProperties();

    // Geometry information
    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix(mat_size, mat_size);

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // Getting density
    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    const IntegrationMethod integration_method =
        IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);
    Matrix J0(dimension, dimension), DN_DX0_deriv;
    const auto& integration_points = r_geom.IntegrationPoints(integration_method);
    for (unsigned point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[point_number], J0);
        const Matrix& rDN_De = r_geom.ShapeFunctionsLocalGradients(integration_method)[point_number];
        GeometricalSensitivityUtility geometrical_sensitivity(J0, rDN_De);
        double detJ0_deriv;
        geometrical_sensitivity.CalculateSensitivity(Deriv, detJ0_deriv, DN_DX0_deriv);
        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ0_deriv) * thickness;
        const Vector& rN = row(Ncontainer, point_number);

        for (unsigned i = 0; i < r_geom.size(); ++i)
        {
            const unsigned index_i = i * dimension;

            for (unsigned j = 0; j < r_geom.size(); ++j)
            {
                const unsigned index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                for (unsigned k = 0; k < dimension; ++k)
                    rMassMatrix(index_i + k, index_j + k) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    KRATOS_ERROR << "You have called to the CalculateKinematicVariables from the base class for solid elements" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Setting the variables for the CL
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // rotate to local axes strain/F
    RotateToLocalAxes(rValues, rThisKinematicVariables);

    // Actually do the computations in the ConstitutiveLaw in local axes
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done

    // We undo the rotation of strain/F, C, stress
    RotateToGlobalAxes(rValues, rThisKinematicVariables);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::BuildRotationSystem(
    BoundedMatrix<double, 3, 3>& rRotationMatrix,
    const SizeType StrainSize
    )
{
    const array_1d<double, 3>& r_local_axis_1 = this->GetValue(LOCAL_AXIS_1);
    array_1d<double, 3> local_axis_2;
    array_1d<double, 3> local_axis_3;

    if (StrainSize == 6) {
        noalias(local_axis_2) = this->GetValue(LOCAL_AXIS_2);
        noalias(local_axis_3) = MathUtils<double>::CrossProduct(r_local_axis_1, local_axis_2);
    } else if (StrainSize == 3) { // we assume xy plane
        local_axis_2[0] = r_local_axis_1[1];
        local_axis_2[1] = -r_local_axis_1[0];
        local_axis_2[2] = 0.0;
        local_axis_3[0] = 0.0;
        local_axis_3[1] = 0.0;
        local_axis_3[2] = 1.0;
    }
    StructuralMechanicsElementUtilities::InitialCheckLocalAxes(r_local_axis_1, local_axis_2, local_axis_3);
    StructuralMechanicsElementUtilities::BuildRotationMatrix(rRotationMatrix, r_local_axis_1, local_axis_2, local_axis_3);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::RotateToLocalAxes(
    ConstitutiveLaw::Parameters& rValues,
    KinematicVariables& rThisKinematicVariables
    )
{
    if (this->IsElementRotated()) {
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
        BoundedMatrix<double, 3, 3> rotation_matrix;

        BuildRotationSystem(rotation_matrix, strain_size);

        if (UseElementProvidedStrain()) { // we rotate strain
            if (strain_size == 6) {
                BoundedMatrix<double, 6, 6> voigt_rotation_matrix;
                ConstitutiveLawUtilities<6>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
                rValues.GetStrainVector() = prod(voigt_rotation_matrix, rValues.GetStrainVector());
            } else if (strain_size == 3) {
                BoundedMatrix<double, 3, 3> voigt_rotation_matrix;
                ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
                rValues.GetStrainVector() = prod(voigt_rotation_matrix, rValues.GetStrainVector());
            }
        } else { // rotate F
            BoundedMatrix<double, 3, 3> inv_rotation_matrix;
            double aux_det;
            MathUtils<double>::InvertMatrix3(rotation_matrix, inv_rotation_matrix, aux_det);
            rThisKinematicVariables.F = prod(rotation_matrix, rThisKinematicVariables.F);
            rThisKinematicVariables.F = prod(rThisKinematicVariables.F, inv_rotation_matrix);
            rValues.SetDeformationGradientF(rThisKinematicVariables.F);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::RotateToGlobalAxes(
    ConstitutiveLaw::Parameters& rValues,
    KinematicVariables& rThisKinematicVariables
    )
{
    if (this->IsElementRotated()) {
        const auto& r_options = rValues.GetOptions();
        const bool stress_option = r_options.Is(ConstitutiveLaw::COMPUTE_STRESS);
        const bool constitutive_matrix_option = r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
        BoundedMatrix<double, 3, 3> rotation_matrix;

        BuildRotationSystem(rotation_matrix, strain_size);

        // Undo the rotation in strain, stress and C
        if (strain_size == 6) {
            BoundedMatrix<double, 6, 6> voigt_rotation_matrix;
            ConstitutiveLawUtilities<6>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
            rValues.GetStrainVector() = prod(trans(voigt_rotation_matrix), rValues.GetStrainVector());
            if (stress_option)
                rValues.GetStressVector() = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            if (constitutive_matrix_option) {
                BoundedMatrix<double, 6, 6> aux;
                noalias(aux) = prod(trans(voigt_rotation_matrix), rValues.GetConstitutiveMatrix());
                noalias(rValues.GetConstitutiveMatrix()) = prod(aux, voigt_rotation_matrix);
            }
        } else if (strain_size == 3) {
            BoundedMatrix<double, 3, 3> voigt_rotation_matrix;
            ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
            rValues.GetStrainVector() = prod(trans(voigt_rotation_matrix), rValues.GetStrainVector());
            if (stress_option)
                rValues.GetStressVector() = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            if (constitutive_matrix_option) {
                BoundedMatrix<double, 3, 3> aux;
                noalias(aux) = prod(trans(voigt_rotation_matrix), rValues.GetConstitutiveMatrix());
                noalias(rValues.GetConstitutiveMatrix()) = prod(aux, voigt_rotation_matrix);
            }
        }
        // Now undo the rotation in F if required
        if (!UseElementProvidedStrain()) {
            BoundedMatrix<double, 3, 3> inv_rotation_matrix;
            double aux_det;
            MathUtils<double>::InvertMatrix3(rotation_matrix, inv_rotation_matrix, aux_det);
            rThisKinematicVariables.F = prod(inv_rotation_matrix, rThisKinematicVariables.F);
            rThisKinematicVariables.F = prod(rThisKinematicVariables.F, rotation_matrix);
            rValues.SetDeformationGradientF(rThisKinematicVariables.F);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // Assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& BaseSolidElement::CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    DeltaDisplacement.resize(number_of_nodes , dimension, false);

    for ( IndexType i_node = 0; i_node < number_of_nodes; i_node++ ) {
        const array_1d<double, 3 >& current_displacement  = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 >& previous_displacement = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType j_dim = 0; j_dim < dimension; ++j_dim )
            DeltaDisplacement(i_node, j_dim) = current_displacement[j_dim] - previous_displacement[j_dim];
    }

    return DeltaDisplacement;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

double BaseSolidElement::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    const GeometryType& r_geom = GetGeometry();
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geom,
        r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    const Matrix& rDN_De =
        GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    return detJ0;
}

/***********************************************************************************/
/***********************************************************************************/

double BaseSolidElement::CalculateDerivativesOnCurrentConfiguration(
    Matrix& rJ,
    Matrix& rInvJ,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    double detJ;
    rJ = GetGeometry().Jacobian( rJ, PointNumber, ThisIntegrationMethod );
    const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );
    GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ, rDN_DX);
    return detJ;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> BaseSolidElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod( trans( B ), Matrix(prod(D, B)));

const double cA0 = B(0,0)*D(0,0) + B(1,0)*D(1,0) + B(2,0)*D(2,0) + B(3,0)*D(3,0) + B(4,0)*D(4,0) + B(5,0)*D(5,0);
const double cA1 = B(0,0)*D(0,1) + B(1,0)*D(1,1) + B(2,0)*D(2,1) + B(3,0)*D(3,1) + B(4,0)*D(4,1) + B(5,0)*D(5,1);
const double cA2 = B(0,0)*D(0,2) + B(1,0)*D(1,2) + B(2,0)*D(2,2) + B(3,0)*D(3,2) + B(4,0)*D(4,2) + B(5,0)*D(5,2);
const double cA3 = B(0,0)*D(0,3) + B(1,0)*D(1,3) + B(2,0)*D(2,3) + B(3,0)*D(3,3) + B(4,0)*D(4,3) + B(5,0)*D(5,3);
const double cA4 = B(0,0)*D(0,4) + B(1,0)*D(1,4) + B(2,0)*D(2,4) + B(3,0)*D(3,4) + B(4,0)*D(4,4) + B(5,0)*D(5,4);
const double cA5 = B(0,0)*D(0,5) + B(1,0)*D(1,5) + B(2,0)*D(2,5) + B(3,0)*D(3,5) + B(4,0)*D(4,5) + B(5,0)*D(5,5);
const double cA6 = B(0,1)*D(0,0) + B(1,1)*D(1,0) + B(2,1)*D(2,0) + B(3,1)*D(3,0) + B(4,1)*D(4,0) + B(5,1)*D(5,0);
const double cA7 = B(0,1)*D(0,1) + B(1,1)*D(1,1) + B(2,1)*D(2,1) + B(3,1)*D(3,1) + B(4,1)*D(4,1) + B(5,1)*D(5,1);
const double cA8 = B(0,1)*D(0,2) + B(1,1)*D(1,2) + B(2,1)*D(2,2) + B(3,1)*D(3,2) + B(4,1)*D(4,2) + B(5,1)*D(5,2);
const double cA9 = B(0,1)*D(0,3) + B(1,1)*D(1,3) + B(2,1)*D(2,3) + B(3,1)*D(3,3) + B(4,1)*D(4,3) + B(5,1)*D(5,3);
const double cA10 = B(0,1)*D(0,4) + B(1,1)*D(1,4) + B(2,1)*D(2,4) + B(3,1)*D(3,4) + B(4,1)*D(4,4) + B(5,1)*D(5,4);
const double cA11 = B(0,1)*D(0,5) + B(1,1)*D(1,5) + B(2,1)*D(2,5) + B(3,1)*D(3,5) + B(4,1)*D(4,5) + B(5,1)*D(5,5);
const double cA12 = B(0,2)*D(0,0) + B(1,2)*D(1,0) + B(2,2)*D(2,0) + B(3,2)*D(3,0) + B(4,2)*D(4,0) + B(5,2)*D(5,0);
const double cA13 = B(0,2)*D(0,1) + B(1,2)*D(1,1) + B(2,2)*D(2,1) + B(3,2)*D(3,1) + B(4,2)*D(4,1) + B(5,2)*D(5,1);
const double cA14 = B(0,2)*D(0,2) + B(1,2)*D(1,2) + B(2,2)*D(2,2) + B(3,2)*D(3,2) + B(4,2)*D(4,2) + B(5,2)*D(5,2);
const double cA15 = B(0,2)*D(0,3) + B(1,2)*D(1,3) + B(2,2)*D(2,3) + B(3,2)*D(3,3) + B(4,2)*D(4,3) + B(5,2)*D(5,3);
const double cA16 = B(0,2)*D(0,4) + B(1,2)*D(1,4) + B(2,2)*D(2,4) + B(3,2)*D(3,4) + B(4,2)*D(4,4) + B(5,2)*D(5,4);
const double cA17 = B(0,2)*D(0,5) + B(1,2)*D(1,5) + B(2,2)*D(2,5) + B(3,2)*D(3,5) + B(4,2)*D(4,5) + B(5,2)*D(5,5);
const double cA18 = B(0,3)*D(0,0) + B(1,3)*D(1,0) + B(2,3)*D(2,0) + B(3,3)*D(3,0) + B(4,3)*D(4,0) + B(5,3)*D(5,0);
const double cA19 = B(0,3)*D(0,1) + B(1,3)*D(1,1) + B(2,3)*D(2,1) + B(3,3)*D(3,1) + B(4,3)*D(4,1) + B(5,3)*D(5,1);
const double cA20 = B(0,3)*D(0,2) + B(1,3)*D(1,2) + B(2,3)*D(2,2) + B(3,3)*D(3,2) + B(4,3)*D(4,2) + B(5,3)*D(5,2);
const double cA21 = B(0,3)*D(0,3) + B(1,3)*D(1,3) + B(2,3)*D(2,3) + B(3,3)*D(3,3) + B(4,3)*D(4,3) + B(5,3)*D(5,3);
const double cA22 = B(0,3)*D(0,4) + B(1,3)*D(1,4) + B(2,3)*D(2,4) + B(3,3)*D(3,4) + B(4,3)*D(4,4) + B(5,3)*D(5,4);
const double cA23 = B(0,3)*D(0,5) + B(1,3)*D(1,5) + B(2,3)*D(2,5) + B(3,3)*D(3,5) + B(4,3)*D(4,5) + B(5,3)*D(5,5);
const double cA24 = B(0,4)*D(0,0) + B(1,4)*D(1,0) + B(2,4)*D(2,0) + B(3,4)*D(3,0) + B(4,4)*D(4,0) + B(5,4)*D(5,0);
const double cA25 = B(0,4)*D(0,1) + B(1,4)*D(1,1) + B(2,4)*D(2,1) + B(3,4)*D(3,1) + B(4,4)*D(4,1) + B(5,4)*D(5,1);
const double cA26 = B(0,4)*D(0,2) + B(1,4)*D(1,2) + B(2,4)*D(2,2) + B(3,4)*D(3,2) + B(4,4)*D(4,2) + B(5,4)*D(5,2);
const double cA27 = B(0,4)*D(0,3) + B(1,4)*D(1,3) + B(2,4)*D(2,3) + B(3,4)*D(3,3) + B(4,4)*D(4,3) + B(5,4)*D(5,3);
const double cA28 = B(0,4)*D(0,4) + B(1,4)*D(1,4) + B(2,4)*D(2,4) + B(3,4)*D(3,4) + B(4,4)*D(4,4) + B(5,4)*D(5,4);
const double cA29 = B(0,4)*D(0,5) + B(1,4)*D(1,5) + B(2,4)*D(2,5) + B(3,4)*D(3,5) + B(4,4)*D(4,5) + B(5,4)*D(5,5);
const double cA30 = B(0,5)*D(0,0) + B(1,5)*D(1,0) + B(2,5)*D(2,0) + B(3,5)*D(3,0) + B(4,5)*D(4,0) + B(5,5)*D(5,0);
const double cA31 = B(0,5)*D(0,1) + B(1,5)*D(1,1) + B(2,5)*D(2,1) + B(3,5)*D(3,1) + B(4,5)*D(4,1) + B(5,5)*D(5,1);
const double cA32 = B(0,5)*D(0,2) + B(1,5)*D(1,2) + B(2,5)*D(2,2) + B(3,5)*D(3,2) + B(4,5)*D(4,2) + B(5,5)*D(5,2);
const double cA33 = B(0,5)*D(0,3) + B(1,5)*D(1,3) + B(2,5)*D(2,3) + B(3,5)*D(3,3) + B(4,5)*D(4,3) + B(5,5)*D(5,3);
const double cA34 = B(0,5)*D(0,4) + B(1,5)*D(1,4) + B(2,5)*D(2,4) + B(3,5)*D(3,4) + B(4,5)*D(4,4) + B(5,5)*D(5,4);
const double cA35 = B(0,5)*D(0,5) + B(1,5)*D(1,5) + B(2,5)*D(2,5) + B(3,5)*D(3,5) + B(4,5)*D(4,5) + B(5,5)*D(5,5);
const double cA36 = B(0,6)*D(0,0) + B(1,6)*D(1,0) + B(2,6)*D(2,0) + B(3,6)*D(3,0) + B(4,6)*D(4,0) + B(5,6)*D(5,0);
const double cA37 = B(0,6)*D(0,1) + B(1,6)*D(1,1) + B(2,6)*D(2,1) + B(3,6)*D(3,1) + B(4,6)*D(4,1) + B(5,6)*D(5,1);
const double cA38 = B(0,6)*D(0,2) + B(1,6)*D(1,2) + B(2,6)*D(2,2) + B(3,6)*D(3,2) + B(4,6)*D(4,2) + B(5,6)*D(5,2);
const double cA39 = B(0,6)*D(0,3) + B(1,6)*D(1,3) + B(2,6)*D(2,3) + B(3,6)*D(3,3) + B(4,6)*D(4,3) + B(5,6)*D(5,3);
const double cA40 = B(0,6)*D(0,4) + B(1,6)*D(1,4) + B(2,6)*D(2,4) + B(3,6)*D(3,4) + B(4,6)*D(4,4) + B(5,6)*D(5,4);
const double cA41 = B(0,6)*D(0,5) + B(1,6)*D(1,5) + B(2,6)*D(2,5) + B(3,6)*D(3,5) + B(4,6)*D(4,5) + B(5,6)*D(5,5);
const double cA42 = B(0,7)*D(0,0) + B(1,7)*D(1,0) + B(2,7)*D(2,0) + B(3,7)*D(3,0) + B(4,7)*D(4,0) + B(5,7)*D(5,0);
const double cA43 = B(0,7)*D(0,1) + B(1,7)*D(1,1) + B(2,7)*D(2,1) + B(3,7)*D(3,1) + B(4,7)*D(4,1) + B(5,7)*D(5,1);
const double cA44 = B(0,7)*D(0,2) + B(1,7)*D(1,2) + B(2,7)*D(2,2) + B(3,7)*D(3,2) + B(4,7)*D(4,2) + B(5,7)*D(5,2);
const double cA45 = B(0,7)*D(0,3) + B(1,7)*D(1,3) + B(2,7)*D(2,3) + B(3,7)*D(3,3) + B(4,7)*D(4,3) + B(5,7)*D(5,3);
const double cA46 = B(0,7)*D(0,4) + B(1,7)*D(1,4) + B(2,7)*D(2,4) + B(3,7)*D(3,4) + B(4,7)*D(4,4) + B(5,7)*D(5,4);
const double cA47 = B(0,7)*D(0,5) + B(1,7)*D(1,5) + B(2,7)*D(2,5) + B(3,7)*D(3,5) + B(4,7)*D(4,5) + B(5,7)*D(5,5);
const double cA48 = B(0,8)*D(0,0) + B(1,8)*D(1,0) + B(2,8)*D(2,0) + B(3,8)*D(3,0) + B(4,8)*D(4,0) + B(5,8)*D(5,0);
const double cA49 = B(0,8)*D(0,1) + B(1,8)*D(1,1) + B(2,8)*D(2,1) + B(3,8)*D(3,1) + B(4,8)*D(4,1) + B(5,8)*D(5,1);
const double cA50 = B(0,8)*D(0,2) + B(1,8)*D(1,2) + B(2,8)*D(2,2) + B(3,8)*D(3,2) + B(4,8)*D(4,2) + B(5,8)*D(5,2);
const double cA51 = B(0,8)*D(0,3) + B(1,8)*D(1,3) + B(2,8)*D(2,3) + B(3,8)*D(3,3) + B(4,8)*D(4,3) + B(5,8)*D(5,3);
const double cA52 = B(0,8)*D(0,4) + B(1,8)*D(1,4) + B(2,8)*D(2,4) + B(3,8)*D(3,4) + B(4,8)*D(4,4) + B(5,8)*D(5,4);
const double cA53 = B(0,8)*D(0,5) + B(1,8)*D(1,5) + B(2,8)*D(2,5) + B(3,8)*D(3,5) + B(4,8)*D(4,5) + B(5,8)*D(5,5);
const double cA54 = B(0,9)*D(0,0) + B(1,9)*D(1,0) + B(2,9)*D(2,0) + B(3,9)*D(3,0) + B(4,9)*D(4,0) + B(5,9)*D(5,0);
const double cA55 = B(0,9)*D(0,1) + B(1,9)*D(1,1) + B(2,9)*D(2,1) + B(3,9)*D(3,1) + B(4,9)*D(4,1) + B(5,9)*D(5,1);
const double cA56 = B(0,9)*D(0,2) + B(1,9)*D(1,2) + B(2,9)*D(2,2) + B(3,9)*D(3,2) + B(4,9)*D(4,2) + B(5,9)*D(5,2);
const double cA57 = B(0,9)*D(0,3) + B(1,9)*D(1,3) + B(2,9)*D(2,3) + B(3,9)*D(3,3) + B(4,9)*D(4,3) + B(5,9)*D(5,3);
const double cA58 = B(0,9)*D(0,4) + B(1,9)*D(1,4) + B(2,9)*D(2,4) + B(3,9)*D(3,4) + B(4,9)*D(4,4) + B(5,9)*D(5,4);
const double cA59 = B(0,9)*D(0,5) + B(1,9)*D(1,5) + B(2,9)*D(2,5) + B(3,9)*D(3,5) + B(4,9)*D(4,5) + B(5,9)*D(5,5);
const double cA60 = B(0,10)*D(0,0) + B(1,10)*D(1,0) + B(2,10)*D(2,0) + B(3,10)*D(3,0) + B(4,10)*D(4,0) + B(5,10)*D(5,0);
const double cA61 = B(0,10)*D(0,1) + B(1,10)*D(1,1) + B(2,10)*D(2,1) + B(3,10)*D(3,1) + B(4,10)*D(4,1) + B(5,10)*D(5,1);
const double cA62 = B(0,10)*D(0,2) + B(1,10)*D(1,2) + B(2,10)*D(2,2) + B(3,10)*D(3,2) + B(4,10)*D(4,2) + B(5,10)*D(5,2);
const double cA63 = B(0,10)*D(0,3) + B(1,10)*D(1,3) + B(2,10)*D(2,3) + B(3,10)*D(3,3) + B(4,10)*D(4,3) + B(5,10)*D(5,3);
const double cA64 = B(0,10)*D(0,4) + B(1,10)*D(1,4) + B(2,10)*D(2,4) + B(3,10)*D(3,4) + B(4,10)*D(4,4) + B(5,10)*D(5,4);
const double cA65 = B(0,10)*D(0,5) + B(1,10)*D(1,5) + B(2,10)*D(2,5) + B(3,10)*D(3,5) + B(4,10)*D(4,5) + B(5,10)*D(5,5);
const double cA66 = B(0,11)*D(0,0) + B(1,11)*D(1,0) + B(2,11)*D(2,0) + B(3,11)*D(3,0) + B(4,11)*D(4,0) + B(5,11)*D(5,0);
const double cA67 = B(0,11)*D(0,1) + B(1,11)*D(1,1) + B(2,11)*D(2,1) + B(3,11)*D(3,1) + B(4,11)*D(4,1) + B(5,11)*D(5,1);
const double cA68 = B(0,11)*D(0,2) + B(1,11)*D(1,2) + B(2,11)*D(2,2) + B(3,11)*D(3,2) + B(4,11)*D(4,2) + B(5,11)*D(5,2);
const double cA69 = B(0,11)*D(0,3) + B(1,11)*D(1,3) + B(2,11)*D(2,3) + B(3,11)*D(3,3) + B(4,11)*D(4,3) + B(5,11)*D(5,3);
const double cA70 = B(0,11)*D(0,4) + B(1,11)*D(1,4) + B(2,11)*D(2,4) + B(3,11)*D(3,4) + B(4,11)*D(4,4) + B(5,11)*D(5,4);
const double cA71 = B(0,11)*D(0,5) + B(1,11)*D(1,5) + B(2,11)*D(2,5) + B(3,11)*D(3,5) + B(4,11)*D(4,5) + B(5,11)*D(5,5);
const double cA72 = B(0,12)*D(0,0) + B(1,12)*D(1,0) + B(2,12)*D(2,0) + B(3,12)*D(3,0) + B(4,12)*D(4,0) + B(5,12)*D(5,0);
const double cA73 = B(0,12)*D(0,1) + B(1,12)*D(1,1) + B(2,12)*D(2,1) + B(3,12)*D(3,1) + B(4,12)*D(4,1) + B(5,12)*D(5,1);
const double cA74 = B(0,12)*D(0,2) + B(1,12)*D(1,2) + B(2,12)*D(2,2) + B(3,12)*D(3,2) + B(4,12)*D(4,2) + B(5,12)*D(5,2);
const double cA75 = B(0,12)*D(0,3) + B(1,12)*D(1,3) + B(2,12)*D(2,3) + B(3,12)*D(3,3) + B(4,12)*D(4,3) + B(5,12)*D(5,3);
const double cA76 = B(0,12)*D(0,4) + B(1,12)*D(1,4) + B(2,12)*D(2,4) + B(3,12)*D(3,4) + B(4,12)*D(4,4) + B(5,12)*D(5,4);
const double cA77 = B(0,12)*D(0,5) + B(1,12)*D(1,5) + B(2,12)*D(2,5) + B(3,12)*D(3,5) + B(4,12)*D(4,5) + B(5,12)*D(5,5);
const double cA78 = B(0,13)*D(0,0) + B(1,13)*D(1,0) + B(2,13)*D(2,0) + B(3,13)*D(3,0) + B(4,13)*D(4,0) + B(5,13)*D(5,0);
const double cA79 = B(0,13)*D(0,1) + B(1,13)*D(1,1) + B(2,13)*D(2,1) + B(3,13)*D(3,1) + B(4,13)*D(4,1) + B(5,13)*D(5,1);
const double cA80 = B(0,13)*D(0,2) + B(1,13)*D(1,2) + B(2,13)*D(2,2) + B(3,13)*D(3,2) + B(4,13)*D(4,2) + B(5,13)*D(5,2);
const double cA81 = B(0,13)*D(0,3) + B(1,13)*D(1,3) + B(2,13)*D(2,3) + B(3,13)*D(3,3) + B(4,13)*D(4,3) + B(5,13)*D(5,3);
const double cA82 = B(0,13)*D(0,4) + B(1,13)*D(1,4) + B(2,13)*D(2,4) + B(3,13)*D(3,4) + B(4,13)*D(4,4) + B(5,13)*D(5,4);
const double cA83 = B(0,13)*D(0,5) + B(1,13)*D(1,5) + B(2,13)*D(2,5) + B(3,13)*D(3,5) + B(4,13)*D(4,5) + B(5,13)*D(5,5);
const double cA84 = B(0,14)*D(0,0) + B(1,14)*D(1,0) + B(2,14)*D(2,0) + B(3,14)*D(3,0) + B(4,14)*D(4,0) + B(5,14)*D(5,0);
const double cA85 = B(0,14)*D(0,1) + B(1,14)*D(1,1) + B(2,14)*D(2,1) + B(3,14)*D(3,1) + B(4,14)*D(4,1) + B(5,14)*D(5,1);
const double cA86 = B(0,14)*D(0,2) + B(1,14)*D(1,2) + B(2,14)*D(2,2) + B(3,14)*D(3,2) + B(4,14)*D(4,2) + B(5,14)*D(5,2);
const double cA87 = B(0,14)*D(0,3) + B(1,14)*D(1,3) + B(2,14)*D(2,3) + B(3,14)*D(3,3) + B(4,14)*D(4,3) + B(5,14)*D(5,3);
const double cA88 = B(0,14)*D(0,4) + B(1,14)*D(1,4) + B(2,14)*D(2,4) + B(3,14)*D(3,4) + B(4,14)*D(4,4) + B(5,14)*D(5,4);
const double cA89 = B(0,14)*D(0,5) + B(1,14)*D(1,5) + B(2,14)*D(2,5) + B(3,14)*D(3,5) + B(4,14)*D(4,5) + B(5,14)*D(5,5);
const double cA90 = B(0,15)*D(0,0) + B(1,15)*D(1,0) + B(2,15)*D(2,0) + B(3,15)*D(3,0) + B(4,15)*D(4,0) + B(5,15)*D(5,0);
const double cA91 = B(0,15)*D(0,1) + B(1,15)*D(1,1) + B(2,15)*D(2,1) + B(3,15)*D(3,1) + B(4,15)*D(4,1) + B(5,15)*D(5,1);
const double cA92 = B(0,15)*D(0,2) + B(1,15)*D(1,2) + B(2,15)*D(2,2) + B(3,15)*D(3,2) + B(4,15)*D(4,2) + B(5,15)*D(5,2);
const double cA93 = B(0,15)*D(0,3) + B(1,15)*D(1,3) + B(2,15)*D(2,3) + B(3,15)*D(3,3) + B(4,15)*D(4,3) + B(5,15)*D(5,3);
const double cA94 = B(0,15)*D(0,4) + B(1,15)*D(1,4) + B(2,15)*D(2,4) + B(3,15)*D(3,4) + B(4,15)*D(4,4) + B(5,15)*D(5,4);
const double cA95 = B(0,15)*D(0,5) + B(1,15)*D(1,5) + B(2,15)*D(2,5) + B(3,15)*D(3,5) + B(4,15)*D(4,5) + B(5,15)*D(5,5);
const double cA96 = B(0,16)*D(0,0) + B(1,16)*D(1,0) + B(2,16)*D(2,0) + B(3,16)*D(3,0) + B(4,16)*D(4,0) + B(5,16)*D(5,0);
const double cA97 = B(0,16)*D(0,1) + B(1,16)*D(1,1) + B(2,16)*D(2,1) + B(3,16)*D(3,1) + B(4,16)*D(4,1) + B(5,16)*D(5,1);
const double cA98 = B(0,16)*D(0,2) + B(1,16)*D(1,2) + B(2,16)*D(2,2) + B(3,16)*D(3,2) + B(4,16)*D(4,2) + B(5,16)*D(5,2);
const double cA99 = B(0,16)*D(0,3) + B(1,16)*D(1,3) + B(2,16)*D(2,3) + B(3,16)*D(3,3) + B(4,16)*D(4,3) + B(5,16)*D(5,3);
const double cA100 = B(0,16)*D(0,4) + B(1,16)*D(1,4) + B(2,16)*D(2,4) + B(3,16)*D(3,4) + B(4,16)*D(4,4) + B(5,16)*D(5,4);
const double cA101 = B(0,16)*D(0,5) + B(1,16)*D(1,5) + B(2,16)*D(2,5) + B(3,16)*D(3,5) + B(4,16)*D(4,5) + B(5,16)*D(5,5);
const double cA102 = B(0,17)*D(0,0) + B(1,17)*D(1,0) + B(2,17)*D(2,0) + B(3,17)*D(3,0) + B(4,17)*D(4,0) + B(5,17)*D(5,0);
const double cA103 = B(0,17)*D(0,1) + B(1,17)*D(1,1) + B(2,17)*D(2,1) + B(3,17)*D(3,1) + B(4,17)*D(4,1) + B(5,17)*D(5,1);
const double cA104 = B(0,17)*D(0,2) + B(1,17)*D(1,2) + B(2,17)*D(2,2) + B(3,17)*D(3,2) + B(4,17)*D(4,2) + B(5,17)*D(5,2);
const double cA105 = B(0,17)*D(0,3) + B(1,17)*D(1,3) + B(2,17)*D(2,3) + B(3,17)*D(3,3) + B(4,17)*D(4,3) + B(5,17)*D(5,3);
const double cA106 = B(0,17)*D(0,4) + B(1,17)*D(1,4) + B(2,17)*D(2,4) + B(3,17)*D(3,4) + B(4,17)*D(4,4) + B(5,17)*D(5,4);
const double cA107 = B(0,17)*D(0,5) + B(1,17)*D(1,5) + B(2,17)*D(2,5) + B(3,17)*D(3,5) + B(4,17)*D(4,5) + B(5,17)*D(5,5);
const double cA108 = B(0,18)*D(0,0) + B(1,18)*D(1,0) + B(2,18)*D(2,0) + B(3,18)*D(3,0) + B(4,18)*D(4,0) + B(5,18)*D(5,0);
const double cA109 = B(0,18)*D(0,1) + B(1,18)*D(1,1) + B(2,18)*D(2,1) + B(3,18)*D(3,1) + B(4,18)*D(4,1) + B(5,18)*D(5,1);
const double cA110 = B(0,18)*D(0,2) + B(1,18)*D(1,2) + B(2,18)*D(2,2) + B(3,18)*D(3,2) + B(4,18)*D(4,2) + B(5,18)*D(5,2);
const double cA111 = B(0,18)*D(0,3) + B(1,18)*D(1,3) + B(2,18)*D(2,3) + B(3,18)*D(3,3) + B(4,18)*D(4,3) + B(5,18)*D(5,3);
const double cA112 = B(0,18)*D(0,4) + B(1,18)*D(1,4) + B(2,18)*D(2,4) + B(3,18)*D(3,4) + B(4,18)*D(4,4) + B(5,18)*D(5,4);
const double cA113 = B(0,18)*D(0,5) + B(1,18)*D(1,5) + B(2,18)*D(2,5) + B(3,18)*D(3,5) + B(4,18)*D(4,5) + B(5,18)*D(5,5);
const double cA114 = B(0,19)*D(0,0) + B(1,19)*D(1,0) + B(2,19)*D(2,0) + B(3,19)*D(3,0) + B(4,19)*D(4,0) + B(5,19)*D(5,0);
const double cA115 = B(0,19)*D(0,1) + B(1,19)*D(1,1) + B(2,19)*D(2,1) + B(3,19)*D(3,1) + B(4,19)*D(4,1) + B(5,19)*D(5,1);
const double cA116 = B(0,19)*D(0,2) + B(1,19)*D(1,2) + B(2,19)*D(2,2) + B(3,19)*D(3,2) + B(4,19)*D(4,2) + B(5,19)*D(5,2);
const double cA117 = B(0,19)*D(0,3) + B(1,19)*D(1,3) + B(2,19)*D(2,3) + B(3,19)*D(3,3) + B(4,19)*D(4,3) + B(5,19)*D(5,3);
const double cA118 = B(0,19)*D(0,4) + B(1,19)*D(1,4) + B(2,19)*D(2,4) + B(3,19)*D(3,4) + B(4,19)*D(4,4) + B(5,19)*D(5,4);
const double cA119 = B(0,19)*D(0,5) + B(1,19)*D(1,5) + B(2,19)*D(2,5) + B(3,19)*D(3,5) + B(4,19)*D(4,5) + B(5,19)*D(5,5);
const double cA120 = B(0,20)*D(0,0) + B(1,20)*D(1,0) + B(2,20)*D(2,0) + B(3,20)*D(3,0) + B(4,20)*D(4,0) + B(5,20)*D(5,0);
const double cA121 = B(0,20)*D(0,1) + B(1,20)*D(1,1) + B(2,20)*D(2,1) + B(3,20)*D(3,1) + B(4,20)*D(4,1) + B(5,20)*D(5,1);
const double cA122 = B(0,20)*D(0,2) + B(1,20)*D(1,2) + B(2,20)*D(2,2) + B(3,20)*D(3,2) + B(4,20)*D(4,2) + B(5,20)*D(5,2);
const double cA123 = B(0,20)*D(0,3) + B(1,20)*D(1,3) + B(2,20)*D(2,3) + B(3,20)*D(3,3) + B(4,20)*D(4,3) + B(5,20)*D(5,3);
const double cA124 = B(0,20)*D(0,4) + B(1,20)*D(1,4) + B(2,20)*D(2,4) + B(3,20)*D(3,4) + B(4,20)*D(4,4) + B(5,20)*D(5,4);
const double cA125 = B(0,20)*D(0,5) + B(1,20)*D(1,5) + B(2,20)*D(2,5) + B(3,20)*D(3,5) + B(4,20)*D(4,5) + B(5,20)*D(5,5);
const double cA126 = B(0,21)*D(0,0) + B(1,21)*D(1,0) + B(2,21)*D(2,0) + B(3,21)*D(3,0) + B(4,21)*D(4,0) + B(5,21)*D(5,0);
const double cA127 = B(0,21)*D(0,1) + B(1,21)*D(1,1) + B(2,21)*D(2,1) + B(3,21)*D(3,1) + B(4,21)*D(4,1) + B(5,21)*D(5,1);
const double cA128 = B(0,21)*D(0,2) + B(1,21)*D(1,2) + B(2,21)*D(2,2) + B(3,21)*D(3,2) + B(4,21)*D(4,2) + B(5,21)*D(5,2);
const double cA129 = B(0,21)*D(0,3) + B(1,21)*D(1,3) + B(2,21)*D(2,3) + B(3,21)*D(3,3) + B(4,21)*D(4,3) + B(5,21)*D(5,3);
const double cA130 = B(0,21)*D(0,4) + B(1,21)*D(1,4) + B(2,21)*D(2,4) + B(3,21)*D(3,4) + B(4,21)*D(4,4) + B(5,21)*D(5,4);
const double cA131 = B(0,21)*D(0,5) + B(1,21)*D(1,5) + B(2,21)*D(2,5) + B(3,21)*D(3,5) + B(4,21)*D(4,5) + B(5,21)*D(5,5);
const double cA132 = B(0,22)*D(0,0) + B(1,22)*D(1,0) + B(2,22)*D(2,0) + B(3,22)*D(3,0) + B(4,22)*D(4,0) + B(5,22)*D(5,0);
const double cA133 = B(0,22)*D(0,1) + B(1,22)*D(1,1) + B(2,22)*D(2,1) + B(3,22)*D(3,1) + B(4,22)*D(4,1) + B(5,22)*D(5,1);
const double cA134 = B(0,22)*D(0,2) + B(1,22)*D(1,2) + B(2,22)*D(2,2) + B(3,22)*D(3,2) + B(4,22)*D(4,2) + B(5,22)*D(5,2);
const double cA135 = B(0,22)*D(0,3) + B(1,22)*D(1,3) + B(2,22)*D(2,3) + B(3,22)*D(3,3) + B(4,22)*D(4,3) + B(5,22)*D(5,3);
const double cA136 = B(0,22)*D(0,4) + B(1,22)*D(1,4) + B(2,22)*D(2,4) + B(3,22)*D(3,4) + B(4,22)*D(4,4) + B(5,22)*D(5,4);
const double cA137 = B(0,22)*D(0,5) + B(1,22)*D(1,5) + B(2,22)*D(2,5) + B(3,22)*D(3,5) + B(4,22)*D(4,5) + B(5,22)*D(5,5);
const double cA138 = B(0,23)*D(0,0) + B(1,23)*D(1,0) + B(2,23)*D(2,0) + B(3,23)*D(3,0) + B(4,23)*D(4,0) + B(5,23)*D(5,0);
const double cA139 = B(0,23)*D(0,1) + B(1,23)*D(1,1) + B(2,23)*D(2,1) + B(3,23)*D(3,1) + B(4,23)*D(4,1) + B(5,23)*D(5,1);
const double cA140 = B(0,23)*D(0,2) + B(1,23)*D(1,2) + B(2,23)*D(2,2) + B(3,23)*D(3,2) + B(4,23)*D(4,2) + B(5,23)*D(5,2);
const double cA141 = B(0,23)*D(0,3) + B(1,23)*D(1,3) + B(2,23)*D(2,3) + B(3,23)*D(3,3) + B(4,23)*D(4,3) + B(5,23)*D(5,3);
const double cA142 = B(0,23)*D(0,4) + B(1,23)*D(1,4) + B(2,23)*D(2,4) + B(3,23)*D(3,4) + B(4,23)*D(4,4) + B(5,23)*D(5,4);
const double cA143 = B(0,23)*D(0,5) + B(1,23)*D(1,5) + B(2,23)*D(2,5) + B(3,23)*D(3,5) + B(4,23)*D(4,5) + B(5,23)*D(5,5);


rLeftHandSideMatrix(0,0)+=B(0,0)*cA0 + B(1,0)*cA1 + B(2,0)*cA2 + B(3,0)*cA3 + B(4,0)*cA4 + B(5,0)*cA5;
rLeftHandSideMatrix(0,1)+=B(0,1)*cA0 + B(1,1)*cA1 + B(2,1)*cA2 + B(3,1)*cA3 + B(4,1)*cA4 + B(5,1)*cA5;
rLeftHandSideMatrix(0,2)+=B(0,2)*cA0 + B(1,2)*cA1 + B(2,2)*cA2 + B(3,2)*cA3 + B(4,2)*cA4 + B(5,2)*cA5;
rLeftHandSideMatrix(0,3)+=B(0,3)*cA0 + B(1,3)*cA1 + B(2,3)*cA2 + B(3,3)*cA3 + B(4,3)*cA4 + B(5,3)*cA5;
rLeftHandSideMatrix(0,4)+=B(0,4)*cA0 + B(1,4)*cA1 + B(2,4)*cA2 + B(3,4)*cA3 + B(4,4)*cA4 + B(5,4)*cA5;
rLeftHandSideMatrix(0,5)+=B(0,5)*cA0 + B(1,5)*cA1 + B(2,5)*cA2 + B(3,5)*cA3 + B(4,5)*cA4 + B(5,5)*cA5;
rLeftHandSideMatrix(0,6)+=B(0,6)*cA0 + B(1,6)*cA1 + B(2,6)*cA2 + B(3,6)*cA3 + B(4,6)*cA4 + B(5,6)*cA5;
rLeftHandSideMatrix(0,7)+=B(0,7)*cA0 + B(1,7)*cA1 + B(2,7)*cA2 + B(3,7)*cA3 + B(4,7)*cA4 + B(5,7)*cA5;
rLeftHandSideMatrix(0,8)+=B(0,8)*cA0 + B(1,8)*cA1 + B(2,8)*cA2 + B(3,8)*cA3 + B(4,8)*cA4 + B(5,8)*cA5;
rLeftHandSideMatrix(0,9)+=B(0,9)*cA0 + B(1,9)*cA1 + B(2,9)*cA2 + B(3,9)*cA3 + B(4,9)*cA4 + B(5,9)*cA5;
rLeftHandSideMatrix(0,10)+=B(0,10)*cA0 + B(1,10)*cA1 + B(2,10)*cA2 + B(3,10)*cA3 + B(4,10)*cA4 + B(5,10)*cA5;
rLeftHandSideMatrix(0,11)+=B(0,11)*cA0 + B(1,11)*cA1 + B(2,11)*cA2 + B(3,11)*cA3 + B(4,11)*cA4 + B(5,11)*cA5;
rLeftHandSideMatrix(0,12)+=B(0,12)*cA0 + B(1,12)*cA1 + B(2,12)*cA2 + B(3,12)*cA3 + B(4,12)*cA4 + B(5,12)*cA5;
rLeftHandSideMatrix(0,13)+=B(0,13)*cA0 + B(1,13)*cA1 + B(2,13)*cA2 + B(3,13)*cA3 + B(4,13)*cA4 + B(5,13)*cA5;
rLeftHandSideMatrix(0,14)+=B(0,14)*cA0 + B(1,14)*cA1 + B(2,14)*cA2 + B(3,14)*cA3 + B(4,14)*cA4 + B(5,14)*cA5;
rLeftHandSideMatrix(0,15)+=B(0,15)*cA0 + B(1,15)*cA1 + B(2,15)*cA2 + B(3,15)*cA3 + B(4,15)*cA4 + B(5,15)*cA5;
rLeftHandSideMatrix(0,16)+=B(0,16)*cA0 + B(1,16)*cA1 + B(2,16)*cA2 + B(3,16)*cA3 + B(4,16)*cA4 + B(5,16)*cA5;
rLeftHandSideMatrix(0,17)+=B(0,17)*cA0 + B(1,17)*cA1 + B(2,17)*cA2 + B(3,17)*cA3 + B(4,17)*cA4 + B(5,17)*cA5;
rLeftHandSideMatrix(0,18)+=B(0,18)*cA0 + B(1,18)*cA1 + B(2,18)*cA2 + B(3,18)*cA3 + B(4,18)*cA4 + B(5,18)*cA5;
rLeftHandSideMatrix(0,19)+=B(0,19)*cA0 + B(1,19)*cA1 + B(2,19)*cA2 + B(3,19)*cA3 + B(4,19)*cA4 + B(5,19)*cA5;
rLeftHandSideMatrix(0,20)+=B(0,20)*cA0 + B(1,20)*cA1 + B(2,20)*cA2 + B(3,20)*cA3 + B(4,20)*cA4 + B(5,20)*cA5;
rLeftHandSideMatrix(0,21)+=B(0,21)*cA0 + B(1,21)*cA1 + B(2,21)*cA2 + B(3,21)*cA3 + B(4,21)*cA4 + B(5,21)*cA5;
rLeftHandSideMatrix(0,22)+=B(0,22)*cA0 + B(1,22)*cA1 + B(2,22)*cA2 + B(3,22)*cA3 + B(4,22)*cA4 + B(5,22)*cA5;
rLeftHandSideMatrix(0,23)+=B(0,23)*cA0 + B(1,23)*cA1 + B(2,23)*cA2 + B(3,23)*cA3 + B(4,23)*cA4 + B(5,23)*cA5;
rLeftHandSideMatrix(1,0)+=B(0,0)*cA6 + B(1,0)*cA7 + B(2,0)*cA8 + B(3,0)*cA9 + B(4,0)*cA10 + B(5,0)*cA11;
rLeftHandSideMatrix(1,1)+=B(0,1)*cA6 + B(1,1)*cA7 + B(2,1)*cA8 + B(3,1)*cA9 + B(4,1)*cA10 + B(5,1)*cA11;
rLeftHandSideMatrix(1,2)+=B(0,2)*cA6 + B(1,2)*cA7 + B(2,2)*cA8 + B(3,2)*cA9 + B(4,2)*cA10 + B(5,2)*cA11;
rLeftHandSideMatrix(1,3)+=B(0,3)*cA6 + B(1,3)*cA7 + B(2,3)*cA8 + B(3,3)*cA9 + B(4,3)*cA10 + B(5,3)*cA11;
rLeftHandSideMatrix(1,4)+=B(0,4)*cA6 + B(1,4)*cA7 + B(2,4)*cA8 + B(3,4)*cA9 + B(4,4)*cA10 + B(5,4)*cA11;
rLeftHandSideMatrix(1,5)+=B(0,5)*cA6 + B(1,5)*cA7 + B(2,5)*cA8 + B(3,5)*cA9 + B(4,5)*cA10 + B(5,5)*cA11;
rLeftHandSideMatrix(1,6)+=B(0,6)*cA6 + B(1,6)*cA7 + B(2,6)*cA8 + B(3,6)*cA9 + B(4,6)*cA10 + B(5,6)*cA11;
rLeftHandSideMatrix(1,7)+=B(0,7)*cA6 + B(1,7)*cA7 + B(2,7)*cA8 + B(3,7)*cA9 + B(4,7)*cA10 + B(5,7)*cA11;
rLeftHandSideMatrix(1,8)+=B(0,8)*cA6 + B(1,8)*cA7 + B(2,8)*cA8 + B(3,8)*cA9 + B(4,8)*cA10 + B(5,8)*cA11;
rLeftHandSideMatrix(1,9)+=B(0,9)*cA6 + B(1,9)*cA7 + B(2,9)*cA8 + B(3,9)*cA9 + B(4,9)*cA10 + B(5,9)*cA11;
rLeftHandSideMatrix(1,10)+=B(0,10)*cA6 + B(1,10)*cA7 + B(2,10)*cA8 + B(3,10)*cA9 + B(4,10)*cA10 + B(5,10)*cA11;
rLeftHandSideMatrix(1,11)+=B(0,11)*cA6 + B(1,11)*cA7 + B(2,11)*cA8 + B(3,11)*cA9 + B(4,11)*cA10 + B(5,11)*cA11;
rLeftHandSideMatrix(1,12)+=B(0,12)*cA6 + B(1,12)*cA7 + B(2,12)*cA8 + B(3,12)*cA9 + B(4,12)*cA10 + B(5,12)*cA11;
rLeftHandSideMatrix(1,13)+=B(0,13)*cA6 + B(1,13)*cA7 + B(2,13)*cA8 + B(3,13)*cA9 + B(4,13)*cA10 + B(5,13)*cA11;
rLeftHandSideMatrix(1,14)+=B(0,14)*cA6 + B(1,14)*cA7 + B(2,14)*cA8 + B(3,14)*cA9 + B(4,14)*cA10 + B(5,14)*cA11;
rLeftHandSideMatrix(1,15)+=B(0,15)*cA6 + B(1,15)*cA7 + B(2,15)*cA8 + B(3,15)*cA9 + B(4,15)*cA10 + B(5,15)*cA11;
rLeftHandSideMatrix(1,16)+=B(0,16)*cA6 + B(1,16)*cA7 + B(2,16)*cA8 + B(3,16)*cA9 + B(4,16)*cA10 + B(5,16)*cA11;
rLeftHandSideMatrix(1,17)+=B(0,17)*cA6 + B(1,17)*cA7 + B(2,17)*cA8 + B(3,17)*cA9 + B(4,17)*cA10 + B(5,17)*cA11;
rLeftHandSideMatrix(1,18)+=B(0,18)*cA6 + B(1,18)*cA7 + B(2,18)*cA8 + B(3,18)*cA9 + B(4,18)*cA10 + B(5,18)*cA11;
rLeftHandSideMatrix(1,19)+=B(0,19)*cA6 + B(1,19)*cA7 + B(2,19)*cA8 + B(3,19)*cA9 + B(4,19)*cA10 + B(5,19)*cA11;
rLeftHandSideMatrix(1,20)+=B(0,20)*cA6 + B(1,20)*cA7 + B(2,20)*cA8 + B(3,20)*cA9 + B(4,20)*cA10 + B(5,20)*cA11;
rLeftHandSideMatrix(1,21)+=B(0,21)*cA6 + B(1,21)*cA7 + B(2,21)*cA8 + B(3,21)*cA9 + B(4,21)*cA10 + B(5,21)*cA11;
rLeftHandSideMatrix(1,22)+=B(0,22)*cA6 + B(1,22)*cA7 + B(2,22)*cA8 + B(3,22)*cA9 + B(4,22)*cA10 + B(5,22)*cA11;
rLeftHandSideMatrix(1,23)+=B(0,23)*cA6 + B(1,23)*cA7 + B(2,23)*cA8 + B(3,23)*cA9 + B(4,23)*cA10 + B(5,23)*cA11;
rLeftHandSideMatrix(2,0)+=B(0,0)*cA12 + B(1,0)*cA13 + B(2,0)*cA14 + B(3,0)*cA15 + B(4,0)*cA16 + B(5,0)*cA17;
rLeftHandSideMatrix(2,1)+=B(0,1)*cA12 + B(1,1)*cA13 + B(2,1)*cA14 + B(3,1)*cA15 + B(4,1)*cA16 + B(5,1)*cA17;
rLeftHandSideMatrix(2,2)+=B(0,2)*cA12 + B(1,2)*cA13 + B(2,2)*cA14 + B(3,2)*cA15 + B(4,2)*cA16 + B(5,2)*cA17;
rLeftHandSideMatrix(2,3)+=B(0,3)*cA12 + B(1,3)*cA13 + B(2,3)*cA14 + B(3,3)*cA15 + B(4,3)*cA16 + B(5,3)*cA17;
rLeftHandSideMatrix(2,4)+=B(0,4)*cA12 + B(1,4)*cA13 + B(2,4)*cA14 + B(3,4)*cA15 + B(4,4)*cA16 + B(5,4)*cA17;
rLeftHandSideMatrix(2,5)+=B(0,5)*cA12 + B(1,5)*cA13 + B(2,5)*cA14 + B(3,5)*cA15 + B(4,5)*cA16 + B(5,5)*cA17;
rLeftHandSideMatrix(2,6)+=B(0,6)*cA12 + B(1,6)*cA13 + B(2,6)*cA14 + B(3,6)*cA15 + B(4,6)*cA16 + B(5,6)*cA17;
rLeftHandSideMatrix(2,7)+=B(0,7)*cA12 + B(1,7)*cA13 + B(2,7)*cA14 + B(3,7)*cA15 + B(4,7)*cA16 + B(5,7)*cA17;
rLeftHandSideMatrix(2,8)+=B(0,8)*cA12 + B(1,8)*cA13 + B(2,8)*cA14 + B(3,8)*cA15 + B(4,8)*cA16 + B(5,8)*cA17;
rLeftHandSideMatrix(2,9)+=B(0,9)*cA12 + B(1,9)*cA13 + B(2,9)*cA14 + B(3,9)*cA15 + B(4,9)*cA16 + B(5,9)*cA17;
rLeftHandSideMatrix(2,10)+=B(0,10)*cA12 + B(1,10)*cA13 + B(2,10)*cA14 + B(3,10)*cA15 + B(4,10)*cA16 + B(5,10)*cA17;
rLeftHandSideMatrix(2,11)+=B(0,11)*cA12 + B(1,11)*cA13 + B(2,11)*cA14 + B(3,11)*cA15 + B(4,11)*cA16 + B(5,11)*cA17;
rLeftHandSideMatrix(2,12)+=B(0,12)*cA12 + B(1,12)*cA13 + B(2,12)*cA14 + B(3,12)*cA15 + B(4,12)*cA16 + B(5,12)*cA17;
rLeftHandSideMatrix(2,13)+=B(0,13)*cA12 + B(1,13)*cA13 + B(2,13)*cA14 + B(3,13)*cA15 + B(4,13)*cA16 + B(5,13)*cA17;
rLeftHandSideMatrix(2,14)+=B(0,14)*cA12 + B(1,14)*cA13 + B(2,14)*cA14 + B(3,14)*cA15 + B(4,14)*cA16 + B(5,14)*cA17;
rLeftHandSideMatrix(2,15)+=B(0,15)*cA12 + B(1,15)*cA13 + B(2,15)*cA14 + B(3,15)*cA15 + B(4,15)*cA16 + B(5,15)*cA17;
rLeftHandSideMatrix(2,16)+=B(0,16)*cA12 + B(1,16)*cA13 + B(2,16)*cA14 + B(3,16)*cA15 + B(4,16)*cA16 + B(5,16)*cA17;
rLeftHandSideMatrix(2,17)+=B(0,17)*cA12 + B(1,17)*cA13 + B(2,17)*cA14 + B(3,17)*cA15 + B(4,17)*cA16 + B(5,17)*cA17;
rLeftHandSideMatrix(2,18)+=B(0,18)*cA12 + B(1,18)*cA13 + B(2,18)*cA14 + B(3,18)*cA15 + B(4,18)*cA16 + B(5,18)*cA17;
rLeftHandSideMatrix(2,19)+=B(0,19)*cA12 + B(1,19)*cA13 + B(2,19)*cA14 + B(3,19)*cA15 + B(4,19)*cA16 + B(5,19)*cA17;
rLeftHandSideMatrix(2,20)+=B(0,20)*cA12 + B(1,20)*cA13 + B(2,20)*cA14 + B(3,20)*cA15 + B(4,20)*cA16 + B(5,20)*cA17;
rLeftHandSideMatrix(2,21)+=B(0,21)*cA12 + B(1,21)*cA13 + B(2,21)*cA14 + B(3,21)*cA15 + B(4,21)*cA16 + B(5,21)*cA17;
rLeftHandSideMatrix(2,22)+=B(0,22)*cA12 + B(1,22)*cA13 + B(2,22)*cA14 + B(3,22)*cA15 + B(4,22)*cA16 + B(5,22)*cA17;
rLeftHandSideMatrix(2,23)+=B(0,23)*cA12 + B(1,23)*cA13 + B(2,23)*cA14 + B(3,23)*cA15 + B(4,23)*cA16 + B(5,23)*cA17;
rLeftHandSideMatrix(3,0)+=B(0,0)*cA18 + B(1,0)*cA19 + B(2,0)*cA20 + B(3,0)*cA21 + B(4,0)*cA22 + B(5,0)*cA23;
rLeftHandSideMatrix(3,1)+=B(0,1)*cA18 + B(1,1)*cA19 + B(2,1)*cA20 + B(3,1)*cA21 + B(4,1)*cA22 + B(5,1)*cA23;
rLeftHandSideMatrix(3,2)+=B(0,2)*cA18 + B(1,2)*cA19 + B(2,2)*cA20 + B(3,2)*cA21 + B(4,2)*cA22 + B(5,2)*cA23;
rLeftHandSideMatrix(3,3)+=B(0,3)*cA18 + B(1,3)*cA19 + B(2,3)*cA20 + B(3,3)*cA21 + B(4,3)*cA22 + B(5,3)*cA23;
rLeftHandSideMatrix(3,4)+=B(0,4)*cA18 + B(1,4)*cA19 + B(2,4)*cA20 + B(3,4)*cA21 + B(4,4)*cA22 + B(5,4)*cA23;
rLeftHandSideMatrix(3,5)+=B(0,5)*cA18 + B(1,5)*cA19 + B(2,5)*cA20 + B(3,5)*cA21 + B(4,5)*cA22 + B(5,5)*cA23;
rLeftHandSideMatrix(3,6)+=B(0,6)*cA18 + B(1,6)*cA19 + B(2,6)*cA20 + B(3,6)*cA21 + B(4,6)*cA22 + B(5,6)*cA23;
rLeftHandSideMatrix(3,7)+=B(0,7)*cA18 + B(1,7)*cA19 + B(2,7)*cA20 + B(3,7)*cA21 + B(4,7)*cA22 + B(5,7)*cA23;
rLeftHandSideMatrix(3,8)+=B(0,8)*cA18 + B(1,8)*cA19 + B(2,8)*cA20 + B(3,8)*cA21 + B(4,8)*cA22 + B(5,8)*cA23;
rLeftHandSideMatrix(3,9)+=B(0,9)*cA18 + B(1,9)*cA19 + B(2,9)*cA20 + B(3,9)*cA21 + B(4,9)*cA22 + B(5,9)*cA23;
rLeftHandSideMatrix(3,10)+=B(0,10)*cA18 + B(1,10)*cA19 + B(2,10)*cA20 + B(3,10)*cA21 + B(4,10)*cA22 + B(5,10)*cA23;
rLeftHandSideMatrix(3,11)+=B(0,11)*cA18 + B(1,11)*cA19 + B(2,11)*cA20 + B(3,11)*cA21 + B(4,11)*cA22 + B(5,11)*cA23;
rLeftHandSideMatrix(3,12)+=B(0,12)*cA18 + B(1,12)*cA19 + B(2,12)*cA20 + B(3,12)*cA21 + B(4,12)*cA22 + B(5,12)*cA23;
rLeftHandSideMatrix(3,13)+=B(0,13)*cA18 + B(1,13)*cA19 + B(2,13)*cA20 + B(3,13)*cA21 + B(4,13)*cA22 + B(5,13)*cA23;
rLeftHandSideMatrix(3,14)+=B(0,14)*cA18 + B(1,14)*cA19 + B(2,14)*cA20 + B(3,14)*cA21 + B(4,14)*cA22 + B(5,14)*cA23;
rLeftHandSideMatrix(3,15)+=B(0,15)*cA18 + B(1,15)*cA19 + B(2,15)*cA20 + B(3,15)*cA21 + B(4,15)*cA22 + B(5,15)*cA23;
rLeftHandSideMatrix(3,16)+=B(0,16)*cA18 + B(1,16)*cA19 + B(2,16)*cA20 + B(3,16)*cA21 + B(4,16)*cA22 + B(5,16)*cA23;
rLeftHandSideMatrix(3,17)+=B(0,17)*cA18 + B(1,17)*cA19 + B(2,17)*cA20 + B(3,17)*cA21 + B(4,17)*cA22 + B(5,17)*cA23;
rLeftHandSideMatrix(3,18)+=B(0,18)*cA18 + B(1,18)*cA19 + B(2,18)*cA20 + B(3,18)*cA21 + B(4,18)*cA22 + B(5,18)*cA23;
rLeftHandSideMatrix(3,19)+=B(0,19)*cA18 + B(1,19)*cA19 + B(2,19)*cA20 + B(3,19)*cA21 + B(4,19)*cA22 + B(5,19)*cA23;
rLeftHandSideMatrix(3,20)+=B(0,20)*cA18 + B(1,20)*cA19 + B(2,20)*cA20 + B(3,20)*cA21 + B(4,20)*cA22 + B(5,20)*cA23;
rLeftHandSideMatrix(3,21)+=B(0,21)*cA18 + B(1,21)*cA19 + B(2,21)*cA20 + B(3,21)*cA21 + B(4,21)*cA22 + B(5,21)*cA23;
rLeftHandSideMatrix(3,22)+=B(0,22)*cA18 + B(1,22)*cA19 + B(2,22)*cA20 + B(3,22)*cA21 + B(4,22)*cA22 + B(5,22)*cA23;
rLeftHandSideMatrix(3,23)+=B(0,23)*cA18 + B(1,23)*cA19 + B(2,23)*cA20 + B(3,23)*cA21 + B(4,23)*cA22 + B(5,23)*cA23;
rLeftHandSideMatrix(4,0)+=B(0,0)*cA24 + B(1,0)*cA25 + B(2,0)*cA26 + B(3,0)*cA27 + B(4,0)*cA28 + B(5,0)*cA29;
rLeftHandSideMatrix(4,1)+=B(0,1)*cA24 + B(1,1)*cA25 + B(2,1)*cA26 + B(3,1)*cA27 + B(4,1)*cA28 + B(5,1)*cA29;
rLeftHandSideMatrix(4,2)+=B(0,2)*cA24 + B(1,2)*cA25 + B(2,2)*cA26 + B(3,2)*cA27 + B(4,2)*cA28 + B(5,2)*cA29;
rLeftHandSideMatrix(4,3)+=B(0,3)*cA24 + B(1,3)*cA25 + B(2,3)*cA26 + B(3,3)*cA27 + B(4,3)*cA28 + B(5,3)*cA29;
rLeftHandSideMatrix(4,4)+=B(0,4)*cA24 + B(1,4)*cA25 + B(2,4)*cA26 + B(3,4)*cA27 + B(4,4)*cA28 + B(5,4)*cA29;
rLeftHandSideMatrix(4,5)+=B(0,5)*cA24 + B(1,5)*cA25 + B(2,5)*cA26 + B(3,5)*cA27 + B(4,5)*cA28 + B(5,5)*cA29;
rLeftHandSideMatrix(4,6)+=B(0,6)*cA24 + B(1,6)*cA25 + B(2,6)*cA26 + B(3,6)*cA27 + B(4,6)*cA28 + B(5,6)*cA29;
rLeftHandSideMatrix(4,7)+=B(0,7)*cA24 + B(1,7)*cA25 + B(2,7)*cA26 + B(3,7)*cA27 + B(4,7)*cA28 + B(5,7)*cA29;
rLeftHandSideMatrix(4,8)+=B(0,8)*cA24 + B(1,8)*cA25 + B(2,8)*cA26 + B(3,8)*cA27 + B(4,8)*cA28 + B(5,8)*cA29;
rLeftHandSideMatrix(4,9)+=B(0,9)*cA24 + B(1,9)*cA25 + B(2,9)*cA26 + B(3,9)*cA27 + B(4,9)*cA28 + B(5,9)*cA29;
rLeftHandSideMatrix(4,10)+=B(0,10)*cA24 + B(1,10)*cA25 + B(2,10)*cA26 + B(3,10)*cA27 + B(4,10)*cA28 + B(5,10)*cA29;
rLeftHandSideMatrix(4,11)+=B(0,11)*cA24 + B(1,11)*cA25 + B(2,11)*cA26 + B(3,11)*cA27 + B(4,11)*cA28 + B(5,11)*cA29;
rLeftHandSideMatrix(4,12)+=B(0,12)*cA24 + B(1,12)*cA25 + B(2,12)*cA26 + B(3,12)*cA27 + B(4,12)*cA28 + B(5,12)*cA29;
rLeftHandSideMatrix(4,13)+=B(0,13)*cA24 + B(1,13)*cA25 + B(2,13)*cA26 + B(3,13)*cA27 + B(4,13)*cA28 + B(5,13)*cA29;
rLeftHandSideMatrix(4,14)+=B(0,14)*cA24 + B(1,14)*cA25 + B(2,14)*cA26 + B(3,14)*cA27 + B(4,14)*cA28 + B(5,14)*cA29;
rLeftHandSideMatrix(4,15)+=B(0,15)*cA24 + B(1,15)*cA25 + B(2,15)*cA26 + B(3,15)*cA27 + B(4,15)*cA28 + B(5,15)*cA29;
rLeftHandSideMatrix(4,16)+=B(0,16)*cA24 + B(1,16)*cA25 + B(2,16)*cA26 + B(3,16)*cA27 + B(4,16)*cA28 + B(5,16)*cA29;
rLeftHandSideMatrix(4,17)+=B(0,17)*cA24 + B(1,17)*cA25 + B(2,17)*cA26 + B(3,17)*cA27 + B(4,17)*cA28 + B(5,17)*cA29;
rLeftHandSideMatrix(4,18)+=B(0,18)*cA24 + B(1,18)*cA25 + B(2,18)*cA26 + B(3,18)*cA27 + B(4,18)*cA28 + B(5,18)*cA29;
rLeftHandSideMatrix(4,19)+=B(0,19)*cA24 + B(1,19)*cA25 + B(2,19)*cA26 + B(3,19)*cA27 + B(4,19)*cA28 + B(5,19)*cA29;
rLeftHandSideMatrix(4,20)+=B(0,20)*cA24 + B(1,20)*cA25 + B(2,20)*cA26 + B(3,20)*cA27 + B(4,20)*cA28 + B(5,20)*cA29;
rLeftHandSideMatrix(4,21)+=B(0,21)*cA24 + B(1,21)*cA25 + B(2,21)*cA26 + B(3,21)*cA27 + B(4,21)*cA28 + B(5,21)*cA29;
rLeftHandSideMatrix(4,22)+=B(0,22)*cA24 + B(1,22)*cA25 + B(2,22)*cA26 + B(3,22)*cA27 + B(4,22)*cA28 + B(5,22)*cA29;
rLeftHandSideMatrix(4,23)+=B(0,23)*cA24 + B(1,23)*cA25 + B(2,23)*cA26 + B(3,23)*cA27 + B(4,23)*cA28 + B(5,23)*cA29;
rLeftHandSideMatrix(5,0)+=B(0,0)*cA30 + B(1,0)*cA31 + B(2,0)*cA32 + B(3,0)*cA33 + B(4,0)*cA34 + B(5,0)*cA35;
rLeftHandSideMatrix(5,1)+=B(0,1)*cA30 + B(1,1)*cA31 + B(2,1)*cA32 + B(3,1)*cA33 + B(4,1)*cA34 + B(5,1)*cA35;
rLeftHandSideMatrix(5,2)+=B(0,2)*cA30 + B(1,2)*cA31 + B(2,2)*cA32 + B(3,2)*cA33 + B(4,2)*cA34 + B(5,2)*cA35;
rLeftHandSideMatrix(5,3)+=B(0,3)*cA30 + B(1,3)*cA31 + B(2,3)*cA32 + B(3,3)*cA33 + B(4,3)*cA34 + B(5,3)*cA35;
rLeftHandSideMatrix(5,4)+=B(0,4)*cA30 + B(1,4)*cA31 + B(2,4)*cA32 + B(3,4)*cA33 + B(4,4)*cA34 + B(5,4)*cA35;
rLeftHandSideMatrix(5,5)+=B(0,5)*cA30 + B(1,5)*cA31 + B(2,5)*cA32 + B(3,5)*cA33 + B(4,5)*cA34 + B(5,5)*cA35;
rLeftHandSideMatrix(5,6)+=B(0,6)*cA30 + B(1,6)*cA31 + B(2,6)*cA32 + B(3,6)*cA33 + B(4,6)*cA34 + B(5,6)*cA35;
rLeftHandSideMatrix(5,7)+=B(0,7)*cA30 + B(1,7)*cA31 + B(2,7)*cA32 + B(3,7)*cA33 + B(4,7)*cA34 + B(5,7)*cA35;
rLeftHandSideMatrix(5,8)+=B(0,8)*cA30 + B(1,8)*cA31 + B(2,8)*cA32 + B(3,8)*cA33 + B(4,8)*cA34 + B(5,8)*cA35;
rLeftHandSideMatrix(5,9)+=B(0,9)*cA30 + B(1,9)*cA31 + B(2,9)*cA32 + B(3,9)*cA33 + B(4,9)*cA34 + B(5,9)*cA35;
rLeftHandSideMatrix(5,10)+=B(0,10)*cA30 + B(1,10)*cA31 + B(2,10)*cA32 + B(3,10)*cA33 + B(4,10)*cA34 + B(5,10)*cA35;
rLeftHandSideMatrix(5,11)+=B(0,11)*cA30 + B(1,11)*cA31 + B(2,11)*cA32 + B(3,11)*cA33 + B(4,11)*cA34 + B(5,11)*cA35;
rLeftHandSideMatrix(5,12)+=B(0,12)*cA30 + B(1,12)*cA31 + B(2,12)*cA32 + B(3,12)*cA33 + B(4,12)*cA34 + B(5,12)*cA35;
rLeftHandSideMatrix(5,13)+=B(0,13)*cA30 + B(1,13)*cA31 + B(2,13)*cA32 + B(3,13)*cA33 + B(4,13)*cA34 + B(5,13)*cA35;
rLeftHandSideMatrix(5,14)+=B(0,14)*cA30 + B(1,14)*cA31 + B(2,14)*cA32 + B(3,14)*cA33 + B(4,14)*cA34 + B(5,14)*cA35;
rLeftHandSideMatrix(5,15)+=B(0,15)*cA30 + B(1,15)*cA31 + B(2,15)*cA32 + B(3,15)*cA33 + B(4,15)*cA34 + B(5,15)*cA35;
rLeftHandSideMatrix(5,16)+=B(0,16)*cA30 + B(1,16)*cA31 + B(2,16)*cA32 + B(3,16)*cA33 + B(4,16)*cA34 + B(5,16)*cA35;
rLeftHandSideMatrix(5,17)+=B(0,17)*cA30 + B(1,17)*cA31 + B(2,17)*cA32 + B(3,17)*cA33 + B(4,17)*cA34 + B(5,17)*cA35;
rLeftHandSideMatrix(5,18)+=B(0,18)*cA30 + B(1,18)*cA31 + B(2,18)*cA32 + B(3,18)*cA33 + B(4,18)*cA34 + B(5,18)*cA35;
rLeftHandSideMatrix(5,19)+=B(0,19)*cA30 + B(1,19)*cA31 + B(2,19)*cA32 + B(3,19)*cA33 + B(4,19)*cA34 + B(5,19)*cA35;
rLeftHandSideMatrix(5,20)+=B(0,20)*cA30 + B(1,20)*cA31 + B(2,20)*cA32 + B(3,20)*cA33 + B(4,20)*cA34 + B(5,20)*cA35;
rLeftHandSideMatrix(5,21)+=B(0,21)*cA30 + B(1,21)*cA31 + B(2,21)*cA32 + B(3,21)*cA33 + B(4,21)*cA34 + B(5,21)*cA35;
rLeftHandSideMatrix(5,22)+=B(0,22)*cA30 + B(1,22)*cA31 + B(2,22)*cA32 + B(3,22)*cA33 + B(4,22)*cA34 + B(5,22)*cA35;
rLeftHandSideMatrix(5,23)+=B(0,23)*cA30 + B(1,23)*cA31 + B(2,23)*cA32 + B(3,23)*cA33 + B(4,23)*cA34 + B(5,23)*cA35;
rLeftHandSideMatrix(6,0)+=B(0,0)*cA36 + B(1,0)*cA37 + B(2,0)*cA38 + B(3,0)*cA39 + B(4,0)*cA40 + B(5,0)*cA41;
rLeftHandSideMatrix(6,1)+=B(0,1)*cA36 + B(1,1)*cA37 + B(2,1)*cA38 + B(3,1)*cA39 + B(4,1)*cA40 + B(5,1)*cA41;
rLeftHandSideMatrix(6,2)+=B(0,2)*cA36 + B(1,2)*cA37 + B(2,2)*cA38 + B(3,2)*cA39 + B(4,2)*cA40 + B(5,2)*cA41;
rLeftHandSideMatrix(6,3)+=B(0,3)*cA36 + B(1,3)*cA37 + B(2,3)*cA38 + B(3,3)*cA39 + B(4,3)*cA40 + B(5,3)*cA41;
rLeftHandSideMatrix(6,4)+=B(0,4)*cA36 + B(1,4)*cA37 + B(2,4)*cA38 + B(3,4)*cA39 + B(4,4)*cA40 + B(5,4)*cA41;
rLeftHandSideMatrix(6,5)+=B(0,5)*cA36 + B(1,5)*cA37 + B(2,5)*cA38 + B(3,5)*cA39 + B(4,5)*cA40 + B(5,5)*cA41;
rLeftHandSideMatrix(6,6)+=B(0,6)*cA36 + B(1,6)*cA37 + B(2,6)*cA38 + B(3,6)*cA39 + B(4,6)*cA40 + B(5,6)*cA41;
rLeftHandSideMatrix(6,7)+=B(0,7)*cA36 + B(1,7)*cA37 + B(2,7)*cA38 + B(3,7)*cA39 + B(4,7)*cA40 + B(5,7)*cA41;
rLeftHandSideMatrix(6,8)+=B(0,8)*cA36 + B(1,8)*cA37 + B(2,8)*cA38 + B(3,8)*cA39 + B(4,8)*cA40 + B(5,8)*cA41;
rLeftHandSideMatrix(6,9)+=B(0,9)*cA36 + B(1,9)*cA37 + B(2,9)*cA38 + B(3,9)*cA39 + B(4,9)*cA40 + B(5,9)*cA41;
rLeftHandSideMatrix(6,10)+=B(0,10)*cA36 + B(1,10)*cA37 + B(2,10)*cA38 + B(3,10)*cA39 + B(4,10)*cA40 + B(5,10)*cA41;
rLeftHandSideMatrix(6,11)+=B(0,11)*cA36 + B(1,11)*cA37 + B(2,11)*cA38 + B(3,11)*cA39 + B(4,11)*cA40 + B(5,11)*cA41;
rLeftHandSideMatrix(6,12)+=B(0,12)*cA36 + B(1,12)*cA37 + B(2,12)*cA38 + B(3,12)*cA39 + B(4,12)*cA40 + B(5,12)*cA41;
rLeftHandSideMatrix(6,13)+=B(0,13)*cA36 + B(1,13)*cA37 + B(2,13)*cA38 + B(3,13)*cA39 + B(4,13)*cA40 + B(5,13)*cA41;
rLeftHandSideMatrix(6,14)+=B(0,14)*cA36 + B(1,14)*cA37 + B(2,14)*cA38 + B(3,14)*cA39 + B(4,14)*cA40 + B(5,14)*cA41;
rLeftHandSideMatrix(6,15)+=B(0,15)*cA36 + B(1,15)*cA37 + B(2,15)*cA38 + B(3,15)*cA39 + B(4,15)*cA40 + B(5,15)*cA41;
rLeftHandSideMatrix(6,16)+=B(0,16)*cA36 + B(1,16)*cA37 + B(2,16)*cA38 + B(3,16)*cA39 + B(4,16)*cA40 + B(5,16)*cA41;
rLeftHandSideMatrix(6,17)+=B(0,17)*cA36 + B(1,17)*cA37 + B(2,17)*cA38 + B(3,17)*cA39 + B(4,17)*cA40 + B(5,17)*cA41;
rLeftHandSideMatrix(6,18)+=B(0,18)*cA36 + B(1,18)*cA37 + B(2,18)*cA38 + B(3,18)*cA39 + B(4,18)*cA40 + B(5,18)*cA41;
rLeftHandSideMatrix(6,19)+=B(0,19)*cA36 + B(1,19)*cA37 + B(2,19)*cA38 + B(3,19)*cA39 + B(4,19)*cA40 + B(5,19)*cA41;
rLeftHandSideMatrix(6,20)+=B(0,20)*cA36 + B(1,20)*cA37 + B(2,20)*cA38 + B(3,20)*cA39 + B(4,20)*cA40 + B(5,20)*cA41;
rLeftHandSideMatrix(6,21)+=B(0,21)*cA36 + B(1,21)*cA37 + B(2,21)*cA38 + B(3,21)*cA39 + B(4,21)*cA40 + B(5,21)*cA41;
rLeftHandSideMatrix(6,22)+=B(0,22)*cA36 + B(1,22)*cA37 + B(2,22)*cA38 + B(3,22)*cA39 + B(4,22)*cA40 + B(5,22)*cA41;
rLeftHandSideMatrix(6,23)+=B(0,23)*cA36 + B(1,23)*cA37 + B(2,23)*cA38 + B(3,23)*cA39 + B(4,23)*cA40 + B(5,23)*cA41;
rLeftHandSideMatrix(7,0)+=B(0,0)*cA42 + B(1,0)*cA43 + B(2,0)*cA44 + B(3,0)*cA45 + B(4,0)*cA46 + B(5,0)*cA47;
rLeftHandSideMatrix(7,1)+=B(0,1)*cA42 + B(1,1)*cA43 + B(2,1)*cA44 + B(3,1)*cA45 + B(4,1)*cA46 + B(5,1)*cA47;
rLeftHandSideMatrix(7,2)+=B(0,2)*cA42 + B(1,2)*cA43 + B(2,2)*cA44 + B(3,2)*cA45 + B(4,2)*cA46 + B(5,2)*cA47;
rLeftHandSideMatrix(7,3)+=B(0,3)*cA42 + B(1,3)*cA43 + B(2,3)*cA44 + B(3,3)*cA45 + B(4,3)*cA46 + B(5,3)*cA47;
rLeftHandSideMatrix(7,4)+=B(0,4)*cA42 + B(1,4)*cA43 + B(2,4)*cA44 + B(3,4)*cA45 + B(4,4)*cA46 + B(5,4)*cA47;
rLeftHandSideMatrix(7,5)+=B(0,5)*cA42 + B(1,5)*cA43 + B(2,5)*cA44 + B(3,5)*cA45 + B(4,5)*cA46 + B(5,5)*cA47;
rLeftHandSideMatrix(7,6)+=B(0,6)*cA42 + B(1,6)*cA43 + B(2,6)*cA44 + B(3,6)*cA45 + B(4,6)*cA46 + B(5,6)*cA47;
rLeftHandSideMatrix(7,7)+=B(0,7)*cA42 + B(1,7)*cA43 + B(2,7)*cA44 + B(3,7)*cA45 + B(4,7)*cA46 + B(5,7)*cA47;
rLeftHandSideMatrix(7,8)+=B(0,8)*cA42 + B(1,8)*cA43 + B(2,8)*cA44 + B(3,8)*cA45 + B(4,8)*cA46 + B(5,8)*cA47;
rLeftHandSideMatrix(7,9)+=B(0,9)*cA42 + B(1,9)*cA43 + B(2,9)*cA44 + B(3,9)*cA45 + B(4,9)*cA46 + B(5,9)*cA47;
rLeftHandSideMatrix(7,10)+=B(0,10)*cA42 + B(1,10)*cA43 + B(2,10)*cA44 + B(3,10)*cA45 + B(4,10)*cA46 + B(5,10)*cA47;
rLeftHandSideMatrix(7,11)+=B(0,11)*cA42 + B(1,11)*cA43 + B(2,11)*cA44 + B(3,11)*cA45 + B(4,11)*cA46 + B(5,11)*cA47;
rLeftHandSideMatrix(7,12)+=B(0,12)*cA42 + B(1,12)*cA43 + B(2,12)*cA44 + B(3,12)*cA45 + B(4,12)*cA46 + B(5,12)*cA47;
rLeftHandSideMatrix(7,13)+=B(0,13)*cA42 + B(1,13)*cA43 + B(2,13)*cA44 + B(3,13)*cA45 + B(4,13)*cA46 + B(5,13)*cA47;
rLeftHandSideMatrix(7,14)+=B(0,14)*cA42 + B(1,14)*cA43 + B(2,14)*cA44 + B(3,14)*cA45 + B(4,14)*cA46 + B(5,14)*cA47;
rLeftHandSideMatrix(7,15)+=B(0,15)*cA42 + B(1,15)*cA43 + B(2,15)*cA44 + B(3,15)*cA45 + B(4,15)*cA46 + B(5,15)*cA47;
rLeftHandSideMatrix(7,16)+=B(0,16)*cA42 + B(1,16)*cA43 + B(2,16)*cA44 + B(3,16)*cA45 + B(4,16)*cA46 + B(5,16)*cA47;
rLeftHandSideMatrix(7,17)+=B(0,17)*cA42 + B(1,17)*cA43 + B(2,17)*cA44 + B(3,17)*cA45 + B(4,17)*cA46 + B(5,17)*cA47;
rLeftHandSideMatrix(7,18)+=B(0,18)*cA42 + B(1,18)*cA43 + B(2,18)*cA44 + B(3,18)*cA45 + B(4,18)*cA46 + B(5,18)*cA47;
rLeftHandSideMatrix(7,19)+=B(0,19)*cA42 + B(1,19)*cA43 + B(2,19)*cA44 + B(3,19)*cA45 + B(4,19)*cA46 + B(5,19)*cA47;
rLeftHandSideMatrix(7,20)+=B(0,20)*cA42 + B(1,20)*cA43 + B(2,20)*cA44 + B(3,20)*cA45 + B(4,20)*cA46 + B(5,20)*cA47;
rLeftHandSideMatrix(7,21)+=B(0,21)*cA42 + B(1,21)*cA43 + B(2,21)*cA44 + B(3,21)*cA45 + B(4,21)*cA46 + B(5,21)*cA47;
rLeftHandSideMatrix(7,22)+=B(0,22)*cA42 + B(1,22)*cA43 + B(2,22)*cA44 + B(3,22)*cA45 + B(4,22)*cA46 + B(5,22)*cA47;
rLeftHandSideMatrix(7,23)+=B(0,23)*cA42 + B(1,23)*cA43 + B(2,23)*cA44 + B(3,23)*cA45 + B(4,23)*cA46 + B(5,23)*cA47;
rLeftHandSideMatrix(8,0)+=B(0,0)*cA48 + B(1,0)*cA49 + B(2,0)*cA50 + B(3,0)*cA51 + B(4,0)*cA52 + B(5,0)*cA53;
rLeftHandSideMatrix(8,1)+=B(0,1)*cA48 + B(1,1)*cA49 + B(2,1)*cA50 + B(3,1)*cA51 + B(4,1)*cA52 + B(5,1)*cA53;
rLeftHandSideMatrix(8,2)+=B(0,2)*cA48 + B(1,2)*cA49 + B(2,2)*cA50 + B(3,2)*cA51 + B(4,2)*cA52 + B(5,2)*cA53;
rLeftHandSideMatrix(8,3)+=B(0,3)*cA48 + B(1,3)*cA49 + B(2,3)*cA50 + B(3,3)*cA51 + B(4,3)*cA52 + B(5,3)*cA53;
rLeftHandSideMatrix(8,4)+=B(0,4)*cA48 + B(1,4)*cA49 + B(2,4)*cA50 + B(3,4)*cA51 + B(4,4)*cA52 + B(5,4)*cA53;
rLeftHandSideMatrix(8,5)+=B(0,5)*cA48 + B(1,5)*cA49 + B(2,5)*cA50 + B(3,5)*cA51 + B(4,5)*cA52 + B(5,5)*cA53;
rLeftHandSideMatrix(8,6)+=B(0,6)*cA48 + B(1,6)*cA49 + B(2,6)*cA50 + B(3,6)*cA51 + B(4,6)*cA52 + B(5,6)*cA53;
rLeftHandSideMatrix(8,7)+=B(0,7)*cA48 + B(1,7)*cA49 + B(2,7)*cA50 + B(3,7)*cA51 + B(4,7)*cA52 + B(5,7)*cA53;
rLeftHandSideMatrix(8,8)+=B(0,8)*cA48 + B(1,8)*cA49 + B(2,8)*cA50 + B(3,8)*cA51 + B(4,8)*cA52 + B(5,8)*cA53;
rLeftHandSideMatrix(8,9)+=B(0,9)*cA48 + B(1,9)*cA49 + B(2,9)*cA50 + B(3,9)*cA51 + B(4,9)*cA52 + B(5,9)*cA53;
rLeftHandSideMatrix(8,10)+=B(0,10)*cA48 + B(1,10)*cA49 + B(2,10)*cA50 + B(3,10)*cA51 + B(4,10)*cA52 + B(5,10)*cA53;
rLeftHandSideMatrix(8,11)+=B(0,11)*cA48 + B(1,11)*cA49 + B(2,11)*cA50 + B(3,11)*cA51 + B(4,11)*cA52 + B(5,11)*cA53;
rLeftHandSideMatrix(8,12)+=B(0,12)*cA48 + B(1,12)*cA49 + B(2,12)*cA50 + B(3,12)*cA51 + B(4,12)*cA52 + B(5,12)*cA53;
rLeftHandSideMatrix(8,13)+=B(0,13)*cA48 + B(1,13)*cA49 + B(2,13)*cA50 + B(3,13)*cA51 + B(4,13)*cA52 + B(5,13)*cA53;
rLeftHandSideMatrix(8,14)+=B(0,14)*cA48 + B(1,14)*cA49 + B(2,14)*cA50 + B(3,14)*cA51 + B(4,14)*cA52 + B(5,14)*cA53;
rLeftHandSideMatrix(8,15)+=B(0,15)*cA48 + B(1,15)*cA49 + B(2,15)*cA50 + B(3,15)*cA51 + B(4,15)*cA52 + B(5,15)*cA53;
rLeftHandSideMatrix(8,16)+=B(0,16)*cA48 + B(1,16)*cA49 + B(2,16)*cA50 + B(3,16)*cA51 + B(4,16)*cA52 + B(5,16)*cA53;
rLeftHandSideMatrix(8,17)+=B(0,17)*cA48 + B(1,17)*cA49 + B(2,17)*cA50 + B(3,17)*cA51 + B(4,17)*cA52 + B(5,17)*cA53;
rLeftHandSideMatrix(8,18)+=B(0,18)*cA48 + B(1,18)*cA49 + B(2,18)*cA50 + B(3,18)*cA51 + B(4,18)*cA52 + B(5,18)*cA53;
rLeftHandSideMatrix(8,19)+=B(0,19)*cA48 + B(1,19)*cA49 + B(2,19)*cA50 + B(3,19)*cA51 + B(4,19)*cA52 + B(5,19)*cA53;
rLeftHandSideMatrix(8,20)+=B(0,20)*cA48 + B(1,20)*cA49 + B(2,20)*cA50 + B(3,20)*cA51 + B(4,20)*cA52 + B(5,20)*cA53;
rLeftHandSideMatrix(8,21)+=B(0,21)*cA48 + B(1,21)*cA49 + B(2,21)*cA50 + B(3,21)*cA51 + B(4,21)*cA52 + B(5,21)*cA53;
rLeftHandSideMatrix(8,22)+=B(0,22)*cA48 + B(1,22)*cA49 + B(2,22)*cA50 + B(3,22)*cA51 + B(4,22)*cA52 + B(5,22)*cA53;
rLeftHandSideMatrix(8,23)+=B(0,23)*cA48 + B(1,23)*cA49 + B(2,23)*cA50 + B(3,23)*cA51 + B(4,23)*cA52 + B(5,23)*cA53;
rLeftHandSideMatrix(9,0)+=B(0,0)*cA54 + B(1,0)*cA55 + B(2,0)*cA56 + B(3,0)*cA57 + B(4,0)*cA58 + B(5,0)*cA59;
rLeftHandSideMatrix(9,1)+=B(0,1)*cA54 + B(1,1)*cA55 + B(2,1)*cA56 + B(3,1)*cA57 + B(4,1)*cA58 + B(5,1)*cA59;
rLeftHandSideMatrix(9,2)+=B(0,2)*cA54 + B(1,2)*cA55 + B(2,2)*cA56 + B(3,2)*cA57 + B(4,2)*cA58 + B(5,2)*cA59;
rLeftHandSideMatrix(9,3)+=B(0,3)*cA54 + B(1,3)*cA55 + B(2,3)*cA56 + B(3,3)*cA57 + B(4,3)*cA58 + B(5,3)*cA59;
rLeftHandSideMatrix(9,4)+=B(0,4)*cA54 + B(1,4)*cA55 + B(2,4)*cA56 + B(3,4)*cA57 + B(4,4)*cA58 + B(5,4)*cA59;
rLeftHandSideMatrix(9,5)+=B(0,5)*cA54 + B(1,5)*cA55 + B(2,5)*cA56 + B(3,5)*cA57 + B(4,5)*cA58 + B(5,5)*cA59;
rLeftHandSideMatrix(9,6)+=B(0,6)*cA54 + B(1,6)*cA55 + B(2,6)*cA56 + B(3,6)*cA57 + B(4,6)*cA58 + B(5,6)*cA59;
rLeftHandSideMatrix(9,7)+=B(0,7)*cA54 + B(1,7)*cA55 + B(2,7)*cA56 + B(3,7)*cA57 + B(4,7)*cA58 + B(5,7)*cA59;
rLeftHandSideMatrix(9,8)+=B(0,8)*cA54 + B(1,8)*cA55 + B(2,8)*cA56 + B(3,8)*cA57 + B(4,8)*cA58 + B(5,8)*cA59;
rLeftHandSideMatrix(9,9)+=B(0,9)*cA54 + B(1,9)*cA55 + B(2,9)*cA56 + B(3,9)*cA57 + B(4,9)*cA58 + B(5,9)*cA59;
rLeftHandSideMatrix(9,10)+=B(0,10)*cA54 + B(1,10)*cA55 + B(2,10)*cA56 + B(3,10)*cA57 + B(4,10)*cA58 + B(5,10)*cA59;
rLeftHandSideMatrix(9,11)+=B(0,11)*cA54 + B(1,11)*cA55 + B(2,11)*cA56 + B(3,11)*cA57 + B(4,11)*cA58 + B(5,11)*cA59;
rLeftHandSideMatrix(9,12)+=B(0,12)*cA54 + B(1,12)*cA55 + B(2,12)*cA56 + B(3,12)*cA57 + B(4,12)*cA58 + B(5,12)*cA59;
rLeftHandSideMatrix(9,13)+=B(0,13)*cA54 + B(1,13)*cA55 + B(2,13)*cA56 + B(3,13)*cA57 + B(4,13)*cA58 + B(5,13)*cA59;
rLeftHandSideMatrix(9,14)+=B(0,14)*cA54 + B(1,14)*cA55 + B(2,14)*cA56 + B(3,14)*cA57 + B(4,14)*cA58 + B(5,14)*cA59;
rLeftHandSideMatrix(9,15)+=B(0,15)*cA54 + B(1,15)*cA55 + B(2,15)*cA56 + B(3,15)*cA57 + B(4,15)*cA58 + B(5,15)*cA59;
rLeftHandSideMatrix(9,16)+=B(0,16)*cA54 + B(1,16)*cA55 + B(2,16)*cA56 + B(3,16)*cA57 + B(4,16)*cA58 + B(5,16)*cA59;
rLeftHandSideMatrix(9,17)+=B(0,17)*cA54 + B(1,17)*cA55 + B(2,17)*cA56 + B(3,17)*cA57 + B(4,17)*cA58 + B(5,17)*cA59;
rLeftHandSideMatrix(9,18)+=B(0,18)*cA54 + B(1,18)*cA55 + B(2,18)*cA56 + B(3,18)*cA57 + B(4,18)*cA58 + B(5,18)*cA59;
rLeftHandSideMatrix(9,19)+=B(0,19)*cA54 + B(1,19)*cA55 + B(2,19)*cA56 + B(3,19)*cA57 + B(4,19)*cA58 + B(5,19)*cA59;
rLeftHandSideMatrix(9,20)+=B(0,20)*cA54 + B(1,20)*cA55 + B(2,20)*cA56 + B(3,20)*cA57 + B(4,20)*cA58 + B(5,20)*cA59;
rLeftHandSideMatrix(9,21)+=B(0,21)*cA54 + B(1,21)*cA55 + B(2,21)*cA56 + B(3,21)*cA57 + B(4,21)*cA58 + B(5,21)*cA59;
rLeftHandSideMatrix(9,22)+=B(0,22)*cA54 + B(1,22)*cA55 + B(2,22)*cA56 + B(3,22)*cA57 + B(4,22)*cA58 + B(5,22)*cA59;
rLeftHandSideMatrix(9,23)+=B(0,23)*cA54 + B(1,23)*cA55 + B(2,23)*cA56 + B(3,23)*cA57 + B(4,23)*cA58 + B(5,23)*cA59;
rLeftHandSideMatrix(10,0)+=B(0,0)*cA60 + B(1,0)*cA61 + B(2,0)*cA62 + B(3,0)*cA63 + B(4,0)*cA64 + B(5,0)*cA65;
rLeftHandSideMatrix(10,1)+=B(0,1)*cA60 + B(1,1)*cA61 + B(2,1)*cA62 + B(3,1)*cA63 + B(4,1)*cA64 + B(5,1)*cA65;
rLeftHandSideMatrix(10,2)+=B(0,2)*cA60 + B(1,2)*cA61 + B(2,2)*cA62 + B(3,2)*cA63 + B(4,2)*cA64 + B(5,2)*cA65;
rLeftHandSideMatrix(10,3)+=B(0,3)*cA60 + B(1,3)*cA61 + B(2,3)*cA62 + B(3,3)*cA63 + B(4,3)*cA64 + B(5,3)*cA65;
rLeftHandSideMatrix(10,4)+=B(0,4)*cA60 + B(1,4)*cA61 + B(2,4)*cA62 + B(3,4)*cA63 + B(4,4)*cA64 + B(5,4)*cA65;
rLeftHandSideMatrix(10,5)+=B(0,5)*cA60 + B(1,5)*cA61 + B(2,5)*cA62 + B(3,5)*cA63 + B(4,5)*cA64 + B(5,5)*cA65;
rLeftHandSideMatrix(10,6)+=B(0,6)*cA60 + B(1,6)*cA61 + B(2,6)*cA62 + B(3,6)*cA63 + B(4,6)*cA64 + B(5,6)*cA65;
rLeftHandSideMatrix(10,7)+=B(0,7)*cA60 + B(1,7)*cA61 + B(2,7)*cA62 + B(3,7)*cA63 + B(4,7)*cA64 + B(5,7)*cA65;
rLeftHandSideMatrix(10,8)+=B(0,8)*cA60 + B(1,8)*cA61 + B(2,8)*cA62 + B(3,8)*cA63 + B(4,8)*cA64 + B(5,8)*cA65;
rLeftHandSideMatrix(10,9)+=B(0,9)*cA60 + B(1,9)*cA61 + B(2,9)*cA62 + B(3,9)*cA63 + B(4,9)*cA64 + B(5,9)*cA65;
rLeftHandSideMatrix(10,10)+=B(0,10)*cA60 + B(1,10)*cA61 + B(2,10)*cA62 + B(3,10)*cA63 + B(4,10)*cA64 + B(5,10)*cA65;
rLeftHandSideMatrix(10,11)+=B(0,11)*cA60 + B(1,11)*cA61 + B(2,11)*cA62 + B(3,11)*cA63 + B(4,11)*cA64 + B(5,11)*cA65;
rLeftHandSideMatrix(10,12)+=B(0,12)*cA60 + B(1,12)*cA61 + B(2,12)*cA62 + B(3,12)*cA63 + B(4,12)*cA64 + B(5,12)*cA65;
rLeftHandSideMatrix(10,13)+=B(0,13)*cA60 + B(1,13)*cA61 + B(2,13)*cA62 + B(3,13)*cA63 + B(4,13)*cA64 + B(5,13)*cA65;
rLeftHandSideMatrix(10,14)+=B(0,14)*cA60 + B(1,14)*cA61 + B(2,14)*cA62 + B(3,14)*cA63 + B(4,14)*cA64 + B(5,14)*cA65;
rLeftHandSideMatrix(10,15)+=B(0,15)*cA60 + B(1,15)*cA61 + B(2,15)*cA62 + B(3,15)*cA63 + B(4,15)*cA64 + B(5,15)*cA65;
rLeftHandSideMatrix(10,16)+=B(0,16)*cA60 + B(1,16)*cA61 + B(2,16)*cA62 + B(3,16)*cA63 + B(4,16)*cA64 + B(5,16)*cA65;
rLeftHandSideMatrix(10,17)+=B(0,17)*cA60 + B(1,17)*cA61 + B(2,17)*cA62 + B(3,17)*cA63 + B(4,17)*cA64 + B(5,17)*cA65;
rLeftHandSideMatrix(10,18)+=B(0,18)*cA60 + B(1,18)*cA61 + B(2,18)*cA62 + B(3,18)*cA63 + B(4,18)*cA64 + B(5,18)*cA65;
rLeftHandSideMatrix(10,19)+=B(0,19)*cA60 + B(1,19)*cA61 + B(2,19)*cA62 + B(3,19)*cA63 + B(4,19)*cA64 + B(5,19)*cA65;
rLeftHandSideMatrix(10,20)+=B(0,20)*cA60 + B(1,20)*cA61 + B(2,20)*cA62 + B(3,20)*cA63 + B(4,20)*cA64 + B(5,20)*cA65;
rLeftHandSideMatrix(10,21)+=B(0,21)*cA60 + B(1,21)*cA61 + B(2,21)*cA62 + B(3,21)*cA63 + B(4,21)*cA64 + B(5,21)*cA65;
rLeftHandSideMatrix(10,22)+=B(0,22)*cA60 + B(1,22)*cA61 + B(2,22)*cA62 + B(3,22)*cA63 + B(4,22)*cA64 + B(5,22)*cA65;
rLeftHandSideMatrix(10,23)+=B(0,23)*cA60 + B(1,23)*cA61 + B(2,23)*cA62 + B(3,23)*cA63 + B(4,23)*cA64 + B(5,23)*cA65;
rLeftHandSideMatrix(11,0)+=B(0,0)*cA66 + B(1,0)*cA67 + B(2,0)*cA68 + B(3,0)*cA69 + B(4,0)*cA70 + B(5,0)*cA71;
rLeftHandSideMatrix(11,1)+=B(0,1)*cA66 + B(1,1)*cA67 + B(2,1)*cA68 + B(3,1)*cA69 + B(4,1)*cA70 + B(5,1)*cA71;
rLeftHandSideMatrix(11,2)+=B(0,2)*cA66 + B(1,2)*cA67 + B(2,2)*cA68 + B(3,2)*cA69 + B(4,2)*cA70 + B(5,2)*cA71;
rLeftHandSideMatrix(11,3)+=B(0,3)*cA66 + B(1,3)*cA67 + B(2,3)*cA68 + B(3,3)*cA69 + B(4,3)*cA70 + B(5,3)*cA71;
rLeftHandSideMatrix(11,4)+=B(0,4)*cA66 + B(1,4)*cA67 + B(2,4)*cA68 + B(3,4)*cA69 + B(4,4)*cA70 + B(5,4)*cA71;
rLeftHandSideMatrix(11,5)+=B(0,5)*cA66 + B(1,5)*cA67 + B(2,5)*cA68 + B(3,5)*cA69 + B(4,5)*cA70 + B(5,5)*cA71;
rLeftHandSideMatrix(11,6)+=B(0,6)*cA66 + B(1,6)*cA67 + B(2,6)*cA68 + B(3,6)*cA69 + B(4,6)*cA70 + B(5,6)*cA71;
rLeftHandSideMatrix(11,7)+=B(0,7)*cA66 + B(1,7)*cA67 + B(2,7)*cA68 + B(3,7)*cA69 + B(4,7)*cA70 + B(5,7)*cA71;
rLeftHandSideMatrix(11,8)+=B(0,8)*cA66 + B(1,8)*cA67 + B(2,8)*cA68 + B(3,8)*cA69 + B(4,8)*cA70 + B(5,8)*cA71;
rLeftHandSideMatrix(11,9)+=B(0,9)*cA66 + B(1,9)*cA67 + B(2,9)*cA68 + B(3,9)*cA69 + B(4,9)*cA70 + B(5,9)*cA71;
rLeftHandSideMatrix(11,10)+=B(0,10)*cA66 + B(1,10)*cA67 + B(2,10)*cA68 + B(3,10)*cA69 + B(4,10)*cA70 + B(5,10)*cA71;
rLeftHandSideMatrix(11,11)+=B(0,11)*cA66 + B(1,11)*cA67 + B(2,11)*cA68 + B(3,11)*cA69 + B(4,11)*cA70 + B(5,11)*cA71;
rLeftHandSideMatrix(11,12)+=B(0,12)*cA66 + B(1,12)*cA67 + B(2,12)*cA68 + B(3,12)*cA69 + B(4,12)*cA70 + B(5,12)*cA71;
rLeftHandSideMatrix(11,13)+=B(0,13)*cA66 + B(1,13)*cA67 + B(2,13)*cA68 + B(3,13)*cA69 + B(4,13)*cA70 + B(5,13)*cA71;
rLeftHandSideMatrix(11,14)+=B(0,14)*cA66 + B(1,14)*cA67 + B(2,14)*cA68 + B(3,14)*cA69 + B(4,14)*cA70 + B(5,14)*cA71;
rLeftHandSideMatrix(11,15)+=B(0,15)*cA66 + B(1,15)*cA67 + B(2,15)*cA68 + B(3,15)*cA69 + B(4,15)*cA70 + B(5,15)*cA71;
rLeftHandSideMatrix(11,16)+=B(0,16)*cA66 + B(1,16)*cA67 + B(2,16)*cA68 + B(3,16)*cA69 + B(4,16)*cA70 + B(5,16)*cA71;
rLeftHandSideMatrix(11,17)+=B(0,17)*cA66 + B(1,17)*cA67 + B(2,17)*cA68 + B(3,17)*cA69 + B(4,17)*cA70 + B(5,17)*cA71;
rLeftHandSideMatrix(11,18)+=B(0,18)*cA66 + B(1,18)*cA67 + B(2,18)*cA68 + B(3,18)*cA69 + B(4,18)*cA70 + B(5,18)*cA71;
rLeftHandSideMatrix(11,19)+=B(0,19)*cA66 + B(1,19)*cA67 + B(2,19)*cA68 + B(3,19)*cA69 + B(4,19)*cA70 + B(5,19)*cA71;
rLeftHandSideMatrix(11,20)+=B(0,20)*cA66 + B(1,20)*cA67 + B(2,20)*cA68 + B(3,20)*cA69 + B(4,20)*cA70 + B(5,20)*cA71;
rLeftHandSideMatrix(11,21)+=B(0,21)*cA66 + B(1,21)*cA67 + B(2,21)*cA68 + B(3,21)*cA69 + B(4,21)*cA70 + B(5,21)*cA71;
rLeftHandSideMatrix(11,22)+=B(0,22)*cA66 + B(1,22)*cA67 + B(2,22)*cA68 + B(3,22)*cA69 + B(4,22)*cA70 + B(5,22)*cA71;
rLeftHandSideMatrix(11,23)+=B(0,23)*cA66 + B(1,23)*cA67 + B(2,23)*cA68 + B(3,23)*cA69 + B(4,23)*cA70 + B(5,23)*cA71;
rLeftHandSideMatrix(12,0)+=B(0,0)*cA72 + B(1,0)*cA73 + B(2,0)*cA74 + B(3,0)*cA75 + B(4,0)*cA76 + B(5,0)*cA77;
rLeftHandSideMatrix(12,1)+=B(0,1)*cA72 + B(1,1)*cA73 + B(2,1)*cA74 + B(3,1)*cA75 + B(4,1)*cA76 + B(5,1)*cA77;
rLeftHandSideMatrix(12,2)+=B(0,2)*cA72 + B(1,2)*cA73 + B(2,2)*cA74 + B(3,2)*cA75 + B(4,2)*cA76 + B(5,2)*cA77;
rLeftHandSideMatrix(12,3)+=B(0,3)*cA72 + B(1,3)*cA73 + B(2,3)*cA74 + B(3,3)*cA75 + B(4,3)*cA76 + B(5,3)*cA77;
rLeftHandSideMatrix(12,4)+=B(0,4)*cA72 + B(1,4)*cA73 + B(2,4)*cA74 + B(3,4)*cA75 + B(4,4)*cA76 + B(5,4)*cA77;
rLeftHandSideMatrix(12,5)+=B(0,5)*cA72 + B(1,5)*cA73 + B(2,5)*cA74 + B(3,5)*cA75 + B(4,5)*cA76 + B(5,5)*cA77;
rLeftHandSideMatrix(12,6)+=B(0,6)*cA72 + B(1,6)*cA73 + B(2,6)*cA74 + B(3,6)*cA75 + B(4,6)*cA76 + B(5,6)*cA77;
rLeftHandSideMatrix(12,7)+=B(0,7)*cA72 + B(1,7)*cA73 + B(2,7)*cA74 + B(3,7)*cA75 + B(4,7)*cA76 + B(5,7)*cA77;
rLeftHandSideMatrix(12,8)+=B(0,8)*cA72 + B(1,8)*cA73 + B(2,8)*cA74 + B(3,8)*cA75 + B(4,8)*cA76 + B(5,8)*cA77;
rLeftHandSideMatrix(12,9)+=B(0,9)*cA72 + B(1,9)*cA73 + B(2,9)*cA74 + B(3,9)*cA75 + B(4,9)*cA76 + B(5,9)*cA77;
rLeftHandSideMatrix(12,10)+=B(0,10)*cA72 + B(1,10)*cA73 + B(2,10)*cA74 + B(3,10)*cA75 + B(4,10)*cA76 + B(5,10)*cA77;
rLeftHandSideMatrix(12,11)+=B(0,11)*cA72 + B(1,11)*cA73 + B(2,11)*cA74 + B(3,11)*cA75 + B(4,11)*cA76 + B(5,11)*cA77;
rLeftHandSideMatrix(12,12)+=B(0,12)*cA72 + B(1,12)*cA73 + B(2,12)*cA74 + B(3,12)*cA75 + B(4,12)*cA76 + B(5,12)*cA77;
rLeftHandSideMatrix(12,13)+=B(0,13)*cA72 + B(1,13)*cA73 + B(2,13)*cA74 + B(3,13)*cA75 + B(4,13)*cA76 + B(5,13)*cA77;
rLeftHandSideMatrix(12,14)+=B(0,14)*cA72 + B(1,14)*cA73 + B(2,14)*cA74 + B(3,14)*cA75 + B(4,14)*cA76 + B(5,14)*cA77;
rLeftHandSideMatrix(12,15)+=B(0,15)*cA72 + B(1,15)*cA73 + B(2,15)*cA74 + B(3,15)*cA75 + B(4,15)*cA76 + B(5,15)*cA77;
rLeftHandSideMatrix(12,16)+=B(0,16)*cA72 + B(1,16)*cA73 + B(2,16)*cA74 + B(3,16)*cA75 + B(4,16)*cA76 + B(5,16)*cA77;
rLeftHandSideMatrix(12,17)+=B(0,17)*cA72 + B(1,17)*cA73 + B(2,17)*cA74 + B(3,17)*cA75 + B(4,17)*cA76 + B(5,17)*cA77;
rLeftHandSideMatrix(12,18)+=B(0,18)*cA72 + B(1,18)*cA73 + B(2,18)*cA74 + B(3,18)*cA75 + B(4,18)*cA76 + B(5,18)*cA77;
rLeftHandSideMatrix(12,19)+=B(0,19)*cA72 + B(1,19)*cA73 + B(2,19)*cA74 + B(3,19)*cA75 + B(4,19)*cA76 + B(5,19)*cA77;
rLeftHandSideMatrix(12,20)+=B(0,20)*cA72 + B(1,20)*cA73 + B(2,20)*cA74 + B(3,20)*cA75 + B(4,20)*cA76 + B(5,20)*cA77;
rLeftHandSideMatrix(12,21)+=B(0,21)*cA72 + B(1,21)*cA73 + B(2,21)*cA74 + B(3,21)*cA75 + B(4,21)*cA76 + B(5,21)*cA77;
rLeftHandSideMatrix(12,22)+=B(0,22)*cA72 + B(1,22)*cA73 + B(2,22)*cA74 + B(3,22)*cA75 + B(4,22)*cA76 + B(5,22)*cA77;
rLeftHandSideMatrix(12,23)+=B(0,23)*cA72 + B(1,23)*cA73 + B(2,23)*cA74 + B(3,23)*cA75 + B(4,23)*cA76 + B(5,23)*cA77;
rLeftHandSideMatrix(13,0)+=B(0,0)*cA78 + B(1,0)*cA79 + B(2,0)*cA80 + B(3,0)*cA81 + B(4,0)*cA82 + B(5,0)*cA83;
rLeftHandSideMatrix(13,1)+=B(0,1)*cA78 + B(1,1)*cA79 + B(2,1)*cA80 + B(3,1)*cA81 + B(4,1)*cA82 + B(5,1)*cA83;
rLeftHandSideMatrix(13,2)+=B(0,2)*cA78 + B(1,2)*cA79 + B(2,2)*cA80 + B(3,2)*cA81 + B(4,2)*cA82 + B(5,2)*cA83;
rLeftHandSideMatrix(13,3)+=B(0,3)*cA78 + B(1,3)*cA79 + B(2,3)*cA80 + B(3,3)*cA81 + B(4,3)*cA82 + B(5,3)*cA83;
rLeftHandSideMatrix(13,4)+=B(0,4)*cA78 + B(1,4)*cA79 + B(2,4)*cA80 + B(3,4)*cA81 + B(4,4)*cA82 + B(5,4)*cA83;
rLeftHandSideMatrix(13,5)+=B(0,5)*cA78 + B(1,5)*cA79 + B(2,5)*cA80 + B(3,5)*cA81 + B(4,5)*cA82 + B(5,5)*cA83;
rLeftHandSideMatrix(13,6)+=B(0,6)*cA78 + B(1,6)*cA79 + B(2,6)*cA80 + B(3,6)*cA81 + B(4,6)*cA82 + B(5,6)*cA83;
rLeftHandSideMatrix(13,7)+=B(0,7)*cA78 + B(1,7)*cA79 + B(2,7)*cA80 + B(3,7)*cA81 + B(4,7)*cA82 + B(5,7)*cA83;
rLeftHandSideMatrix(13,8)+=B(0,8)*cA78 + B(1,8)*cA79 + B(2,8)*cA80 + B(3,8)*cA81 + B(4,8)*cA82 + B(5,8)*cA83;
rLeftHandSideMatrix(13,9)+=B(0,9)*cA78 + B(1,9)*cA79 + B(2,9)*cA80 + B(3,9)*cA81 + B(4,9)*cA82 + B(5,9)*cA83;
rLeftHandSideMatrix(13,10)+=B(0,10)*cA78 + B(1,10)*cA79 + B(2,10)*cA80 + B(3,10)*cA81 + B(4,10)*cA82 + B(5,10)*cA83;
rLeftHandSideMatrix(13,11)+=B(0,11)*cA78 + B(1,11)*cA79 + B(2,11)*cA80 + B(3,11)*cA81 + B(4,11)*cA82 + B(5,11)*cA83;
rLeftHandSideMatrix(13,12)+=B(0,12)*cA78 + B(1,12)*cA79 + B(2,12)*cA80 + B(3,12)*cA81 + B(4,12)*cA82 + B(5,12)*cA83;
rLeftHandSideMatrix(13,13)+=B(0,13)*cA78 + B(1,13)*cA79 + B(2,13)*cA80 + B(3,13)*cA81 + B(4,13)*cA82 + B(5,13)*cA83;
rLeftHandSideMatrix(13,14)+=B(0,14)*cA78 + B(1,14)*cA79 + B(2,14)*cA80 + B(3,14)*cA81 + B(4,14)*cA82 + B(5,14)*cA83;
rLeftHandSideMatrix(13,15)+=B(0,15)*cA78 + B(1,15)*cA79 + B(2,15)*cA80 + B(3,15)*cA81 + B(4,15)*cA82 + B(5,15)*cA83;
rLeftHandSideMatrix(13,16)+=B(0,16)*cA78 + B(1,16)*cA79 + B(2,16)*cA80 + B(3,16)*cA81 + B(4,16)*cA82 + B(5,16)*cA83;
rLeftHandSideMatrix(13,17)+=B(0,17)*cA78 + B(1,17)*cA79 + B(2,17)*cA80 + B(3,17)*cA81 + B(4,17)*cA82 + B(5,17)*cA83;
rLeftHandSideMatrix(13,18)+=B(0,18)*cA78 + B(1,18)*cA79 + B(2,18)*cA80 + B(3,18)*cA81 + B(4,18)*cA82 + B(5,18)*cA83;
rLeftHandSideMatrix(13,19)+=B(0,19)*cA78 + B(1,19)*cA79 + B(2,19)*cA80 + B(3,19)*cA81 + B(4,19)*cA82 + B(5,19)*cA83;
rLeftHandSideMatrix(13,20)+=B(0,20)*cA78 + B(1,20)*cA79 + B(2,20)*cA80 + B(3,20)*cA81 + B(4,20)*cA82 + B(5,20)*cA83;
rLeftHandSideMatrix(13,21)+=B(0,21)*cA78 + B(1,21)*cA79 + B(2,21)*cA80 + B(3,21)*cA81 + B(4,21)*cA82 + B(5,21)*cA83;
rLeftHandSideMatrix(13,22)+=B(0,22)*cA78 + B(1,22)*cA79 + B(2,22)*cA80 + B(3,22)*cA81 + B(4,22)*cA82 + B(5,22)*cA83;
rLeftHandSideMatrix(13,23)+=B(0,23)*cA78 + B(1,23)*cA79 + B(2,23)*cA80 + B(3,23)*cA81 + B(4,23)*cA82 + B(5,23)*cA83;
rLeftHandSideMatrix(14,0)+=B(0,0)*cA84 + B(1,0)*cA85 + B(2,0)*cA86 + B(3,0)*cA87 + B(4,0)*cA88 + B(5,0)*cA89;
rLeftHandSideMatrix(14,1)+=B(0,1)*cA84 + B(1,1)*cA85 + B(2,1)*cA86 + B(3,1)*cA87 + B(4,1)*cA88 + B(5,1)*cA89;
rLeftHandSideMatrix(14,2)+=B(0,2)*cA84 + B(1,2)*cA85 + B(2,2)*cA86 + B(3,2)*cA87 + B(4,2)*cA88 + B(5,2)*cA89;
rLeftHandSideMatrix(14,3)+=B(0,3)*cA84 + B(1,3)*cA85 + B(2,3)*cA86 + B(3,3)*cA87 + B(4,3)*cA88 + B(5,3)*cA89;
rLeftHandSideMatrix(14,4)+=B(0,4)*cA84 + B(1,4)*cA85 + B(2,4)*cA86 + B(3,4)*cA87 + B(4,4)*cA88 + B(5,4)*cA89;
rLeftHandSideMatrix(14,5)+=B(0,5)*cA84 + B(1,5)*cA85 + B(2,5)*cA86 + B(3,5)*cA87 + B(4,5)*cA88 + B(5,5)*cA89;
rLeftHandSideMatrix(14,6)+=B(0,6)*cA84 + B(1,6)*cA85 + B(2,6)*cA86 + B(3,6)*cA87 + B(4,6)*cA88 + B(5,6)*cA89;
rLeftHandSideMatrix(14,7)+=B(0,7)*cA84 + B(1,7)*cA85 + B(2,7)*cA86 + B(3,7)*cA87 + B(4,7)*cA88 + B(5,7)*cA89;
rLeftHandSideMatrix(14,8)+=B(0,8)*cA84 + B(1,8)*cA85 + B(2,8)*cA86 + B(3,8)*cA87 + B(4,8)*cA88 + B(5,8)*cA89;
rLeftHandSideMatrix(14,9)+=B(0,9)*cA84 + B(1,9)*cA85 + B(2,9)*cA86 + B(3,9)*cA87 + B(4,9)*cA88 + B(5,9)*cA89;
rLeftHandSideMatrix(14,10)+=B(0,10)*cA84 + B(1,10)*cA85 + B(2,10)*cA86 + B(3,10)*cA87 + B(4,10)*cA88 + B(5,10)*cA89;
rLeftHandSideMatrix(14,11)+=B(0,11)*cA84 + B(1,11)*cA85 + B(2,11)*cA86 + B(3,11)*cA87 + B(4,11)*cA88 + B(5,11)*cA89;
rLeftHandSideMatrix(14,12)+=B(0,12)*cA84 + B(1,12)*cA85 + B(2,12)*cA86 + B(3,12)*cA87 + B(4,12)*cA88 + B(5,12)*cA89;
rLeftHandSideMatrix(14,13)+=B(0,13)*cA84 + B(1,13)*cA85 + B(2,13)*cA86 + B(3,13)*cA87 + B(4,13)*cA88 + B(5,13)*cA89;
rLeftHandSideMatrix(14,14)+=B(0,14)*cA84 + B(1,14)*cA85 + B(2,14)*cA86 + B(3,14)*cA87 + B(4,14)*cA88 + B(5,14)*cA89;
rLeftHandSideMatrix(14,15)+=B(0,15)*cA84 + B(1,15)*cA85 + B(2,15)*cA86 + B(3,15)*cA87 + B(4,15)*cA88 + B(5,15)*cA89;
rLeftHandSideMatrix(14,16)+=B(0,16)*cA84 + B(1,16)*cA85 + B(2,16)*cA86 + B(3,16)*cA87 + B(4,16)*cA88 + B(5,16)*cA89;
rLeftHandSideMatrix(14,17)+=B(0,17)*cA84 + B(1,17)*cA85 + B(2,17)*cA86 + B(3,17)*cA87 + B(4,17)*cA88 + B(5,17)*cA89;
rLeftHandSideMatrix(14,18)+=B(0,18)*cA84 + B(1,18)*cA85 + B(2,18)*cA86 + B(3,18)*cA87 + B(4,18)*cA88 + B(5,18)*cA89;
rLeftHandSideMatrix(14,19)+=B(0,19)*cA84 + B(1,19)*cA85 + B(2,19)*cA86 + B(3,19)*cA87 + B(4,19)*cA88 + B(5,19)*cA89;
rLeftHandSideMatrix(14,20)+=B(0,20)*cA84 + B(1,20)*cA85 + B(2,20)*cA86 + B(3,20)*cA87 + B(4,20)*cA88 + B(5,20)*cA89;
rLeftHandSideMatrix(14,21)+=B(0,21)*cA84 + B(1,21)*cA85 + B(2,21)*cA86 + B(3,21)*cA87 + B(4,21)*cA88 + B(5,21)*cA89;
rLeftHandSideMatrix(14,22)+=B(0,22)*cA84 + B(1,22)*cA85 + B(2,22)*cA86 + B(3,22)*cA87 + B(4,22)*cA88 + B(5,22)*cA89;
rLeftHandSideMatrix(14,23)+=B(0,23)*cA84 + B(1,23)*cA85 + B(2,23)*cA86 + B(3,23)*cA87 + B(4,23)*cA88 + B(5,23)*cA89;
rLeftHandSideMatrix(15,0)+=B(0,0)*cA90 + B(1,0)*cA91 + B(2,0)*cA92 + B(3,0)*cA93 + B(4,0)*cA94 + B(5,0)*cA95;
rLeftHandSideMatrix(15,1)+=B(0,1)*cA90 + B(1,1)*cA91 + B(2,1)*cA92 + B(3,1)*cA93 + B(4,1)*cA94 + B(5,1)*cA95;
rLeftHandSideMatrix(15,2)+=B(0,2)*cA90 + B(1,2)*cA91 + B(2,2)*cA92 + B(3,2)*cA93 + B(4,2)*cA94 + B(5,2)*cA95;
rLeftHandSideMatrix(15,3)+=B(0,3)*cA90 + B(1,3)*cA91 + B(2,3)*cA92 + B(3,3)*cA93 + B(4,3)*cA94 + B(5,3)*cA95;
rLeftHandSideMatrix(15,4)+=B(0,4)*cA90 + B(1,4)*cA91 + B(2,4)*cA92 + B(3,4)*cA93 + B(4,4)*cA94 + B(5,4)*cA95;
rLeftHandSideMatrix(15,5)+=B(0,5)*cA90 + B(1,5)*cA91 + B(2,5)*cA92 + B(3,5)*cA93 + B(4,5)*cA94 + B(5,5)*cA95;
rLeftHandSideMatrix(15,6)+=B(0,6)*cA90 + B(1,6)*cA91 + B(2,6)*cA92 + B(3,6)*cA93 + B(4,6)*cA94 + B(5,6)*cA95;
rLeftHandSideMatrix(15,7)+=B(0,7)*cA90 + B(1,7)*cA91 + B(2,7)*cA92 + B(3,7)*cA93 + B(4,7)*cA94 + B(5,7)*cA95;
rLeftHandSideMatrix(15,8)+=B(0,8)*cA90 + B(1,8)*cA91 + B(2,8)*cA92 + B(3,8)*cA93 + B(4,8)*cA94 + B(5,8)*cA95;
rLeftHandSideMatrix(15,9)+=B(0,9)*cA90 + B(1,9)*cA91 + B(2,9)*cA92 + B(3,9)*cA93 + B(4,9)*cA94 + B(5,9)*cA95;
rLeftHandSideMatrix(15,10)+=B(0,10)*cA90 + B(1,10)*cA91 + B(2,10)*cA92 + B(3,10)*cA93 + B(4,10)*cA94 + B(5,10)*cA95;
rLeftHandSideMatrix(15,11)+=B(0,11)*cA90 + B(1,11)*cA91 + B(2,11)*cA92 + B(3,11)*cA93 + B(4,11)*cA94 + B(5,11)*cA95;
rLeftHandSideMatrix(15,12)+=B(0,12)*cA90 + B(1,12)*cA91 + B(2,12)*cA92 + B(3,12)*cA93 + B(4,12)*cA94 + B(5,12)*cA95;
rLeftHandSideMatrix(15,13)+=B(0,13)*cA90 + B(1,13)*cA91 + B(2,13)*cA92 + B(3,13)*cA93 + B(4,13)*cA94 + B(5,13)*cA95;
rLeftHandSideMatrix(15,14)+=B(0,14)*cA90 + B(1,14)*cA91 + B(2,14)*cA92 + B(3,14)*cA93 + B(4,14)*cA94 + B(5,14)*cA95;
rLeftHandSideMatrix(15,15)+=B(0,15)*cA90 + B(1,15)*cA91 + B(2,15)*cA92 + B(3,15)*cA93 + B(4,15)*cA94 + B(5,15)*cA95;
rLeftHandSideMatrix(15,16)+=B(0,16)*cA90 + B(1,16)*cA91 + B(2,16)*cA92 + B(3,16)*cA93 + B(4,16)*cA94 + B(5,16)*cA95;
rLeftHandSideMatrix(15,17)+=B(0,17)*cA90 + B(1,17)*cA91 + B(2,17)*cA92 + B(3,17)*cA93 + B(4,17)*cA94 + B(5,17)*cA95;
rLeftHandSideMatrix(15,18)+=B(0,18)*cA90 + B(1,18)*cA91 + B(2,18)*cA92 + B(3,18)*cA93 + B(4,18)*cA94 + B(5,18)*cA95;
rLeftHandSideMatrix(15,19)+=B(0,19)*cA90 + B(1,19)*cA91 + B(2,19)*cA92 + B(3,19)*cA93 + B(4,19)*cA94 + B(5,19)*cA95;
rLeftHandSideMatrix(15,20)+=B(0,20)*cA90 + B(1,20)*cA91 + B(2,20)*cA92 + B(3,20)*cA93 + B(4,20)*cA94 + B(5,20)*cA95;
rLeftHandSideMatrix(15,21)+=B(0,21)*cA90 + B(1,21)*cA91 + B(2,21)*cA92 + B(3,21)*cA93 + B(4,21)*cA94 + B(5,21)*cA95;
rLeftHandSideMatrix(15,22)+=B(0,22)*cA90 + B(1,22)*cA91 + B(2,22)*cA92 + B(3,22)*cA93 + B(4,22)*cA94 + B(5,22)*cA95;
rLeftHandSideMatrix(15,23)+=B(0,23)*cA90 + B(1,23)*cA91 + B(2,23)*cA92 + B(3,23)*cA93 + B(4,23)*cA94 + B(5,23)*cA95;
rLeftHandSideMatrix(16,0)+=B(0,0)*cA96 + B(1,0)*cA97 + B(2,0)*cA98 + B(3,0)*cA99 + B(4,0)*cA100 + B(5,0)*cA101;
rLeftHandSideMatrix(16,1)+=B(0,1)*cA96 + B(1,1)*cA97 + B(2,1)*cA98 + B(3,1)*cA99 + B(4,1)*cA100 + B(5,1)*cA101;
rLeftHandSideMatrix(16,2)+=B(0,2)*cA96 + B(1,2)*cA97 + B(2,2)*cA98 + B(3,2)*cA99 + B(4,2)*cA100 + B(5,2)*cA101;
rLeftHandSideMatrix(16,3)+=B(0,3)*cA96 + B(1,3)*cA97 + B(2,3)*cA98 + B(3,3)*cA99 + B(4,3)*cA100 + B(5,3)*cA101;
rLeftHandSideMatrix(16,4)+=B(0,4)*cA96 + B(1,4)*cA97 + B(2,4)*cA98 + B(3,4)*cA99 + B(4,4)*cA100 + B(5,4)*cA101;
rLeftHandSideMatrix(16,5)+=B(0,5)*cA96 + B(1,5)*cA97 + B(2,5)*cA98 + B(3,5)*cA99 + B(4,5)*cA100 + B(5,5)*cA101;
rLeftHandSideMatrix(16,6)+=B(0,6)*cA96 + B(1,6)*cA97 + B(2,6)*cA98 + B(3,6)*cA99 + B(4,6)*cA100 + B(5,6)*cA101;
rLeftHandSideMatrix(16,7)+=B(0,7)*cA96 + B(1,7)*cA97 + B(2,7)*cA98 + B(3,7)*cA99 + B(4,7)*cA100 + B(5,7)*cA101;
rLeftHandSideMatrix(16,8)+=B(0,8)*cA96 + B(1,8)*cA97 + B(2,8)*cA98 + B(3,8)*cA99 + B(4,8)*cA100 + B(5,8)*cA101;
rLeftHandSideMatrix(16,9)+=B(0,9)*cA96 + B(1,9)*cA97 + B(2,9)*cA98 + B(3,9)*cA99 + B(4,9)*cA100 + B(5,9)*cA101;
rLeftHandSideMatrix(16,10)+=B(0,10)*cA96 + B(1,10)*cA97 + B(2,10)*cA98 + B(3,10)*cA99 + B(4,10)*cA100 + B(5,10)*cA101;
rLeftHandSideMatrix(16,11)+=B(0,11)*cA96 + B(1,11)*cA97 + B(2,11)*cA98 + B(3,11)*cA99 + B(4,11)*cA100 + B(5,11)*cA101;
rLeftHandSideMatrix(16,12)+=B(0,12)*cA96 + B(1,12)*cA97 + B(2,12)*cA98 + B(3,12)*cA99 + B(4,12)*cA100 + B(5,12)*cA101;
rLeftHandSideMatrix(16,13)+=B(0,13)*cA96 + B(1,13)*cA97 + B(2,13)*cA98 + B(3,13)*cA99 + B(4,13)*cA100 + B(5,13)*cA101;
rLeftHandSideMatrix(16,14)+=B(0,14)*cA96 + B(1,14)*cA97 + B(2,14)*cA98 + B(3,14)*cA99 + B(4,14)*cA100 + B(5,14)*cA101;
rLeftHandSideMatrix(16,15)+=B(0,15)*cA96 + B(1,15)*cA97 + B(2,15)*cA98 + B(3,15)*cA99 + B(4,15)*cA100 + B(5,15)*cA101;
rLeftHandSideMatrix(16,16)+=B(0,16)*cA96 + B(1,16)*cA97 + B(2,16)*cA98 + B(3,16)*cA99 + B(4,16)*cA100 + B(5,16)*cA101;
rLeftHandSideMatrix(16,17)+=B(0,17)*cA96 + B(1,17)*cA97 + B(2,17)*cA98 + B(3,17)*cA99 + B(4,17)*cA100 + B(5,17)*cA101;
rLeftHandSideMatrix(16,18)+=B(0,18)*cA96 + B(1,18)*cA97 + B(2,18)*cA98 + B(3,18)*cA99 + B(4,18)*cA100 + B(5,18)*cA101;
rLeftHandSideMatrix(16,19)+=B(0,19)*cA96 + B(1,19)*cA97 + B(2,19)*cA98 + B(3,19)*cA99 + B(4,19)*cA100 + B(5,19)*cA101;
rLeftHandSideMatrix(16,20)+=B(0,20)*cA96 + B(1,20)*cA97 + B(2,20)*cA98 + B(3,20)*cA99 + B(4,20)*cA100 + B(5,20)*cA101;
rLeftHandSideMatrix(16,21)+=B(0,21)*cA96 + B(1,21)*cA97 + B(2,21)*cA98 + B(3,21)*cA99 + B(4,21)*cA100 + B(5,21)*cA101;
rLeftHandSideMatrix(16,22)+=B(0,22)*cA96 + B(1,22)*cA97 + B(2,22)*cA98 + B(3,22)*cA99 + B(4,22)*cA100 + B(5,22)*cA101;
rLeftHandSideMatrix(16,23)+=B(0,23)*cA96 + B(1,23)*cA97 + B(2,23)*cA98 + B(3,23)*cA99 + B(4,23)*cA100 + B(5,23)*cA101;
rLeftHandSideMatrix(17,0)+=B(0,0)*cA102 + B(1,0)*cA103 + B(2,0)*cA104 + B(3,0)*cA105 + B(4,0)*cA106 + B(5,0)*cA107;
rLeftHandSideMatrix(17,1)+=B(0,1)*cA102 + B(1,1)*cA103 + B(2,1)*cA104 + B(3,1)*cA105 + B(4,1)*cA106 + B(5,1)*cA107;
rLeftHandSideMatrix(17,2)+=B(0,2)*cA102 + B(1,2)*cA103 + B(2,2)*cA104 + B(3,2)*cA105 + B(4,2)*cA106 + B(5,2)*cA107;
rLeftHandSideMatrix(17,3)+=B(0,3)*cA102 + B(1,3)*cA103 + B(2,3)*cA104 + B(3,3)*cA105 + B(4,3)*cA106 + B(5,3)*cA107;
rLeftHandSideMatrix(17,4)+=B(0,4)*cA102 + B(1,4)*cA103 + B(2,4)*cA104 + B(3,4)*cA105 + B(4,4)*cA106 + B(5,4)*cA107;
rLeftHandSideMatrix(17,5)+=B(0,5)*cA102 + B(1,5)*cA103 + B(2,5)*cA104 + B(3,5)*cA105 + B(4,5)*cA106 + B(5,5)*cA107;
rLeftHandSideMatrix(17,6)+=B(0,6)*cA102 + B(1,6)*cA103 + B(2,6)*cA104 + B(3,6)*cA105 + B(4,6)*cA106 + B(5,6)*cA107;
rLeftHandSideMatrix(17,7)+=B(0,7)*cA102 + B(1,7)*cA103 + B(2,7)*cA104 + B(3,7)*cA105 + B(4,7)*cA106 + B(5,7)*cA107;
rLeftHandSideMatrix(17,8)+=B(0,8)*cA102 + B(1,8)*cA103 + B(2,8)*cA104 + B(3,8)*cA105 + B(4,8)*cA106 + B(5,8)*cA107;
rLeftHandSideMatrix(17,9)+=B(0,9)*cA102 + B(1,9)*cA103 + B(2,9)*cA104 + B(3,9)*cA105 + B(4,9)*cA106 + B(5,9)*cA107;
rLeftHandSideMatrix(17,10)+=B(0,10)*cA102 + B(1,10)*cA103 + B(2,10)*cA104 + B(3,10)*cA105 + B(4,10)*cA106 + B(5,10)*cA107;
rLeftHandSideMatrix(17,11)+=B(0,11)*cA102 + B(1,11)*cA103 + B(2,11)*cA104 + B(3,11)*cA105 + B(4,11)*cA106 + B(5,11)*cA107;
rLeftHandSideMatrix(17,12)+=B(0,12)*cA102 + B(1,12)*cA103 + B(2,12)*cA104 + B(3,12)*cA105 + B(4,12)*cA106 + B(5,12)*cA107;
rLeftHandSideMatrix(17,13)+=B(0,13)*cA102 + B(1,13)*cA103 + B(2,13)*cA104 + B(3,13)*cA105 + B(4,13)*cA106 + B(5,13)*cA107;
rLeftHandSideMatrix(17,14)+=B(0,14)*cA102 + B(1,14)*cA103 + B(2,14)*cA104 + B(3,14)*cA105 + B(4,14)*cA106 + B(5,14)*cA107;
rLeftHandSideMatrix(17,15)+=B(0,15)*cA102 + B(1,15)*cA103 + B(2,15)*cA104 + B(3,15)*cA105 + B(4,15)*cA106 + B(5,15)*cA107;
rLeftHandSideMatrix(17,16)+=B(0,16)*cA102 + B(1,16)*cA103 + B(2,16)*cA104 + B(3,16)*cA105 + B(4,16)*cA106 + B(5,16)*cA107;
rLeftHandSideMatrix(17,17)+=B(0,17)*cA102 + B(1,17)*cA103 + B(2,17)*cA104 + B(3,17)*cA105 + B(4,17)*cA106 + B(5,17)*cA107;
rLeftHandSideMatrix(17,18)+=B(0,18)*cA102 + B(1,18)*cA103 + B(2,18)*cA104 + B(3,18)*cA105 + B(4,18)*cA106 + B(5,18)*cA107;
rLeftHandSideMatrix(17,19)+=B(0,19)*cA102 + B(1,19)*cA103 + B(2,19)*cA104 + B(3,19)*cA105 + B(4,19)*cA106 + B(5,19)*cA107;
rLeftHandSideMatrix(17,20)+=B(0,20)*cA102 + B(1,20)*cA103 + B(2,20)*cA104 + B(3,20)*cA105 + B(4,20)*cA106 + B(5,20)*cA107;
rLeftHandSideMatrix(17,21)+=B(0,21)*cA102 + B(1,21)*cA103 + B(2,21)*cA104 + B(3,21)*cA105 + B(4,21)*cA106 + B(5,21)*cA107;
rLeftHandSideMatrix(17,22)+=B(0,22)*cA102 + B(1,22)*cA103 + B(2,22)*cA104 + B(3,22)*cA105 + B(4,22)*cA106 + B(5,22)*cA107;
rLeftHandSideMatrix(17,23)+=B(0,23)*cA102 + B(1,23)*cA103 + B(2,23)*cA104 + B(3,23)*cA105 + B(4,23)*cA106 + B(5,23)*cA107;
rLeftHandSideMatrix(18,0)+=B(0,0)*cA108 + B(1,0)*cA109 + B(2,0)*cA110 + B(3,0)*cA111 + B(4,0)*cA112 + B(5,0)*cA113;
rLeftHandSideMatrix(18,1)+=B(0,1)*cA108 + B(1,1)*cA109 + B(2,1)*cA110 + B(3,1)*cA111 + B(4,1)*cA112 + B(5,1)*cA113;
rLeftHandSideMatrix(18,2)+=B(0,2)*cA108 + B(1,2)*cA109 + B(2,2)*cA110 + B(3,2)*cA111 + B(4,2)*cA112 + B(5,2)*cA113;
rLeftHandSideMatrix(18,3)+=B(0,3)*cA108 + B(1,3)*cA109 + B(2,3)*cA110 + B(3,3)*cA111 + B(4,3)*cA112 + B(5,3)*cA113;
rLeftHandSideMatrix(18,4)+=B(0,4)*cA108 + B(1,4)*cA109 + B(2,4)*cA110 + B(3,4)*cA111 + B(4,4)*cA112 + B(5,4)*cA113;
rLeftHandSideMatrix(18,5)+=B(0,5)*cA108 + B(1,5)*cA109 + B(2,5)*cA110 + B(3,5)*cA111 + B(4,5)*cA112 + B(5,5)*cA113;
rLeftHandSideMatrix(18,6)+=B(0,6)*cA108 + B(1,6)*cA109 + B(2,6)*cA110 + B(3,6)*cA111 + B(4,6)*cA112 + B(5,6)*cA113;
rLeftHandSideMatrix(18,7)+=B(0,7)*cA108 + B(1,7)*cA109 + B(2,7)*cA110 + B(3,7)*cA111 + B(4,7)*cA112 + B(5,7)*cA113;
rLeftHandSideMatrix(18,8)+=B(0,8)*cA108 + B(1,8)*cA109 + B(2,8)*cA110 + B(3,8)*cA111 + B(4,8)*cA112 + B(5,8)*cA113;
rLeftHandSideMatrix(18,9)+=B(0,9)*cA108 + B(1,9)*cA109 + B(2,9)*cA110 + B(3,9)*cA111 + B(4,9)*cA112 + B(5,9)*cA113;
rLeftHandSideMatrix(18,10)+=B(0,10)*cA108 + B(1,10)*cA109 + B(2,10)*cA110 + B(3,10)*cA111 + B(4,10)*cA112 + B(5,10)*cA113;
rLeftHandSideMatrix(18,11)+=B(0,11)*cA108 + B(1,11)*cA109 + B(2,11)*cA110 + B(3,11)*cA111 + B(4,11)*cA112 + B(5,11)*cA113;
rLeftHandSideMatrix(18,12)+=B(0,12)*cA108 + B(1,12)*cA109 + B(2,12)*cA110 + B(3,12)*cA111 + B(4,12)*cA112 + B(5,12)*cA113;
rLeftHandSideMatrix(18,13)+=B(0,13)*cA108 + B(1,13)*cA109 + B(2,13)*cA110 + B(3,13)*cA111 + B(4,13)*cA112 + B(5,13)*cA113;
rLeftHandSideMatrix(18,14)+=B(0,14)*cA108 + B(1,14)*cA109 + B(2,14)*cA110 + B(3,14)*cA111 + B(4,14)*cA112 + B(5,14)*cA113;
rLeftHandSideMatrix(18,15)+=B(0,15)*cA108 + B(1,15)*cA109 + B(2,15)*cA110 + B(3,15)*cA111 + B(4,15)*cA112 + B(5,15)*cA113;
rLeftHandSideMatrix(18,16)+=B(0,16)*cA108 + B(1,16)*cA109 + B(2,16)*cA110 + B(3,16)*cA111 + B(4,16)*cA112 + B(5,16)*cA113;
rLeftHandSideMatrix(18,17)+=B(0,17)*cA108 + B(1,17)*cA109 + B(2,17)*cA110 + B(3,17)*cA111 + B(4,17)*cA112 + B(5,17)*cA113;
rLeftHandSideMatrix(18,18)+=B(0,18)*cA108 + B(1,18)*cA109 + B(2,18)*cA110 + B(3,18)*cA111 + B(4,18)*cA112 + B(5,18)*cA113;
rLeftHandSideMatrix(18,19)+=B(0,19)*cA108 + B(1,19)*cA109 + B(2,19)*cA110 + B(3,19)*cA111 + B(4,19)*cA112 + B(5,19)*cA113;
rLeftHandSideMatrix(18,20)+=B(0,20)*cA108 + B(1,20)*cA109 + B(2,20)*cA110 + B(3,20)*cA111 + B(4,20)*cA112 + B(5,20)*cA113;
rLeftHandSideMatrix(18,21)+=B(0,21)*cA108 + B(1,21)*cA109 + B(2,21)*cA110 + B(3,21)*cA111 + B(4,21)*cA112 + B(5,21)*cA113;
rLeftHandSideMatrix(18,22)+=B(0,22)*cA108 + B(1,22)*cA109 + B(2,22)*cA110 + B(3,22)*cA111 + B(4,22)*cA112 + B(5,22)*cA113;
rLeftHandSideMatrix(18,23)+=B(0,23)*cA108 + B(1,23)*cA109 + B(2,23)*cA110 + B(3,23)*cA111 + B(4,23)*cA112 + B(5,23)*cA113;
rLeftHandSideMatrix(19,0)+=B(0,0)*cA114 + B(1,0)*cA115 + B(2,0)*cA116 + B(3,0)*cA117 + B(4,0)*cA118 + B(5,0)*cA119;
rLeftHandSideMatrix(19,1)+=B(0,1)*cA114 + B(1,1)*cA115 + B(2,1)*cA116 + B(3,1)*cA117 + B(4,1)*cA118 + B(5,1)*cA119;
rLeftHandSideMatrix(19,2)+=B(0,2)*cA114 + B(1,2)*cA115 + B(2,2)*cA116 + B(3,2)*cA117 + B(4,2)*cA118 + B(5,2)*cA119;
rLeftHandSideMatrix(19,3)+=B(0,3)*cA114 + B(1,3)*cA115 + B(2,3)*cA116 + B(3,3)*cA117 + B(4,3)*cA118 + B(5,3)*cA119;
rLeftHandSideMatrix(19,4)+=B(0,4)*cA114 + B(1,4)*cA115 + B(2,4)*cA116 + B(3,4)*cA117 + B(4,4)*cA118 + B(5,4)*cA119;
rLeftHandSideMatrix(19,5)+=B(0,5)*cA114 + B(1,5)*cA115 + B(2,5)*cA116 + B(3,5)*cA117 + B(4,5)*cA118 + B(5,5)*cA119;
rLeftHandSideMatrix(19,6)+=B(0,6)*cA114 + B(1,6)*cA115 + B(2,6)*cA116 + B(3,6)*cA117 + B(4,6)*cA118 + B(5,6)*cA119;
rLeftHandSideMatrix(19,7)+=B(0,7)*cA114 + B(1,7)*cA115 + B(2,7)*cA116 + B(3,7)*cA117 + B(4,7)*cA118 + B(5,7)*cA119;
rLeftHandSideMatrix(19,8)+=B(0,8)*cA114 + B(1,8)*cA115 + B(2,8)*cA116 + B(3,8)*cA117 + B(4,8)*cA118 + B(5,8)*cA119;
rLeftHandSideMatrix(19,9)+=B(0,9)*cA114 + B(1,9)*cA115 + B(2,9)*cA116 + B(3,9)*cA117 + B(4,9)*cA118 + B(5,9)*cA119;
rLeftHandSideMatrix(19,10)+=B(0,10)*cA114 + B(1,10)*cA115 + B(2,10)*cA116 + B(3,10)*cA117 + B(4,10)*cA118 + B(5,10)*cA119;
rLeftHandSideMatrix(19,11)+=B(0,11)*cA114 + B(1,11)*cA115 + B(2,11)*cA116 + B(3,11)*cA117 + B(4,11)*cA118 + B(5,11)*cA119;
rLeftHandSideMatrix(19,12)+=B(0,12)*cA114 + B(1,12)*cA115 + B(2,12)*cA116 + B(3,12)*cA117 + B(4,12)*cA118 + B(5,12)*cA119;
rLeftHandSideMatrix(19,13)+=B(0,13)*cA114 + B(1,13)*cA115 + B(2,13)*cA116 + B(3,13)*cA117 + B(4,13)*cA118 + B(5,13)*cA119;
rLeftHandSideMatrix(19,14)+=B(0,14)*cA114 + B(1,14)*cA115 + B(2,14)*cA116 + B(3,14)*cA117 + B(4,14)*cA118 + B(5,14)*cA119;
rLeftHandSideMatrix(19,15)+=B(0,15)*cA114 + B(1,15)*cA115 + B(2,15)*cA116 + B(3,15)*cA117 + B(4,15)*cA118 + B(5,15)*cA119;
rLeftHandSideMatrix(19,16)+=B(0,16)*cA114 + B(1,16)*cA115 + B(2,16)*cA116 + B(3,16)*cA117 + B(4,16)*cA118 + B(5,16)*cA119;
rLeftHandSideMatrix(19,17)+=B(0,17)*cA114 + B(1,17)*cA115 + B(2,17)*cA116 + B(3,17)*cA117 + B(4,17)*cA118 + B(5,17)*cA119;
rLeftHandSideMatrix(19,18)+=B(0,18)*cA114 + B(1,18)*cA115 + B(2,18)*cA116 + B(3,18)*cA117 + B(4,18)*cA118 + B(5,18)*cA119;
rLeftHandSideMatrix(19,19)+=B(0,19)*cA114 + B(1,19)*cA115 + B(2,19)*cA116 + B(3,19)*cA117 + B(4,19)*cA118 + B(5,19)*cA119;
rLeftHandSideMatrix(19,20)+=B(0,20)*cA114 + B(1,20)*cA115 + B(2,20)*cA116 + B(3,20)*cA117 + B(4,20)*cA118 + B(5,20)*cA119;
rLeftHandSideMatrix(19,21)+=B(0,21)*cA114 + B(1,21)*cA115 + B(2,21)*cA116 + B(3,21)*cA117 + B(4,21)*cA118 + B(5,21)*cA119;
rLeftHandSideMatrix(19,22)+=B(0,22)*cA114 + B(1,22)*cA115 + B(2,22)*cA116 + B(3,22)*cA117 + B(4,22)*cA118 + B(5,22)*cA119;
rLeftHandSideMatrix(19,23)+=B(0,23)*cA114 + B(1,23)*cA115 + B(2,23)*cA116 + B(3,23)*cA117 + B(4,23)*cA118 + B(5,23)*cA119;
rLeftHandSideMatrix(20,0)+=B(0,0)*cA120 + B(1,0)*cA121 + B(2,0)*cA122 + B(3,0)*cA123 + B(4,0)*cA124 + B(5,0)*cA125;
rLeftHandSideMatrix(20,1)+=B(0,1)*cA120 + B(1,1)*cA121 + B(2,1)*cA122 + B(3,1)*cA123 + B(4,1)*cA124 + B(5,1)*cA125;
rLeftHandSideMatrix(20,2)+=B(0,2)*cA120 + B(1,2)*cA121 + B(2,2)*cA122 + B(3,2)*cA123 + B(4,2)*cA124 + B(5,2)*cA125;
rLeftHandSideMatrix(20,3)+=B(0,3)*cA120 + B(1,3)*cA121 + B(2,3)*cA122 + B(3,3)*cA123 + B(4,3)*cA124 + B(5,3)*cA125;
rLeftHandSideMatrix(20,4)+=B(0,4)*cA120 + B(1,4)*cA121 + B(2,4)*cA122 + B(3,4)*cA123 + B(4,4)*cA124 + B(5,4)*cA125;
rLeftHandSideMatrix(20,5)+=B(0,5)*cA120 + B(1,5)*cA121 + B(2,5)*cA122 + B(3,5)*cA123 + B(4,5)*cA124 + B(5,5)*cA125;
rLeftHandSideMatrix(20,6)+=B(0,6)*cA120 + B(1,6)*cA121 + B(2,6)*cA122 + B(3,6)*cA123 + B(4,6)*cA124 + B(5,6)*cA125;
rLeftHandSideMatrix(20,7)+=B(0,7)*cA120 + B(1,7)*cA121 + B(2,7)*cA122 + B(3,7)*cA123 + B(4,7)*cA124 + B(5,7)*cA125;
rLeftHandSideMatrix(20,8)+=B(0,8)*cA120 + B(1,8)*cA121 + B(2,8)*cA122 + B(3,8)*cA123 + B(4,8)*cA124 + B(5,8)*cA125;
rLeftHandSideMatrix(20,9)+=B(0,9)*cA120 + B(1,9)*cA121 + B(2,9)*cA122 + B(3,9)*cA123 + B(4,9)*cA124 + B(5,9)*cA125;
rLeftHandSideMatrix(20,10)+=B(0,10)*cA120 + B(1,10)*cA121 + B(2,10)*cA122 + B(3,10)*cA123 + B(4,10)*cA124 + B(5,10)*cA125;
rLeftHandSideMatrix(20,11)+=B(0,11)*cA120 + B(1,11)*cA121 + B(2,11)*cA122 + B(3,11)*cA123 + B(4,11)*cA124 + B(5,11)*cA125;
rLeftHandSideMatrix(20,12)+=B(0,12)*cA120 + B(1,12)*cA121 + B(2,12)*cA122 + B(3,12)*cA123 + B(4,12)*cA124 + B(5,12)*cA125;
rLeftHandSideMatrix(20,13)+=B(0,13)*cA120 + B(1,13)*cA121 + B(2,13)*cA122 + B(3,13)*cA123 + B(4,13)*cA124 + B(5,13)*cA125;
rLeftHandSideMatrix(20,14)+=B(0,14)*cA120 + B(1,14)*cA121 + B(2,14)*cA122 + B(3,14)*cA123 + B(4,14)*cA124 + B(5,14)*cA125;
rLeftHandSideMatrix(20,15)+=B(0,15)*cA120 + B(1,15)*cA121 + B(2,15)*cA122 + B(3,15)*cA123 + B(4,15)*cA124 + B(5,15)*cA125;
rLeftHandSideMatrix(20,16)+=B(0,16)*cA120 + B(1,16)*cA121 + B(2,16)*cA122 + B(3,16)*cA123 + B(4,16)*cA124 + B(5,16)*cA125;
rLeftHandSideMatrix(20,17)+=B(0,17)*cA120 + B(1,17)*cA121 + B(2,17)*cA122 + B(3,17)*cA123 + B(4,17)*cA124 + B(5,17)*cA125;
rLeftHandSideMatrix(20,18)+=B(0,18)*cA120 + B(1,18)*cA121 + B(2,18)*cA122 + B(3,18)*cA123 + B(4,18)*cA124 + B(5,18)*cA125;
rLeftHandSideMatrix(20,19)+=B(0,19)*cA120 + B(1,19)*cA121 + B(2,19)*cA122 + B(3,19)*cA123 + B(4,19)*cA124 + B(5,19)*cA125;
rLeftHandSideMatrix(20,20)+=B(0,20)*cA120 + B(1,20)*cA121 + B(2,20)*cA122 + B(3,20)*cA123 + B(4,20)*cA124 + B(5,20)*cA125;
rLeftHandSideMatrix(20,21)+=B(0,21)*cA120 + B(1,21)*cA121 + B(2,21)*cA122 + B(3,21)*cA123 + B(4,21)*cA124 + B(5,21)*cA125;
rLeftHandSideMatrix(20,22)+=B(0,22)*cA120 + B(1,22)*cA121 + B(2,22)*cA122 + B(3,22)*cA123 + B(4,22)*cA124 + B(5,22)*cA125;
rLeftHandSideMatrix(20,23)+=B(0,23)*cA120 + B(1,23)*cA121 + B(2,23)*cA122 + B(3,23)*cA123 + B(4,23)*cA124 + B(5,23)*cA125;
rLeftHandSideMatrix(21,0)+=B(0,0)*cA126 + B(1,0)*cA127 + B(2,0)*cA128 + B(3,0)*cA129 + B(4,0)*cA130 + B(5,0)*cA131;
rLeftHandSideMatrix(21,1)+=B(0,1)*cA126 + B(1,1)*cA127 + B(2,1)*cA128 + B(3,1)*cA129 + B(4,1)*cA130 + B(5,1)*cA131;
rLeftHandSideMatrix(21,2)+=B(0,2)*cA126 + B(1,2)*cA127 + B(2,2)*cA128 + B(3,2)*cA129 + B(4,2)*cA130 + B(5,2)*cA131;
rLeftHandSideMatrix(21,3)+=B(0,3)*cA126 + B(1,3)*cA127 + B(2,3)*cA128 + B(3,3)*cA129 + B(4,3)*cA130 + B(5,3)*cA131;
rLeftHandSideMatrix(21,4)+=B(0,4)*cA126 + B(1,4)*cA127 + B(2,4)*cA128 + B(3,4)*cA129 + B(4,4)*cA130 + B(5,4)*cA131;
rLeftHandSideMatrix(21,5)+=B(0,5)*cA126 + B(1,5)*cA127 + B(2,5)*cA128 + B(3,5)*cA129 + B(4,5)*cA130 + B(5,5)*cA131;
rLeftHandSideMatrix(21,6)+=B(0,6)*cA126 + B(1,6)*cA127 + B(2,6)*cA128 + B(3,6)*cA129 + B(4,6)*cA130 + B(5,6)*cA131;
rLeftHandSideMatrix(21,7)+=B(0,7)*cA126 + B(1,7)*cA127 + B(2,7)*cA128 + B(3,7)*cA129 + B(4,7)*cA130 + B(5,7)*cA131;
rLeftHandSideMatrix(21,8)+=B(0,8)*cA126 + B(1,8)*cA127 + B(2,8)*cA128 + B(3,8)*cA129 + B(4,8)*cA130 + B(5,8)*cA131;
rLeftHandSideMatrix(21,9)+=B(0,9)*cA126 + B(1,9)*cA127 + B(2,9)*cA128 + B(3,9)*cA129 + B(4,9)*cA130 + B(5,9)*cA131;
rLeftHandSideMatrix(21,10)+=B(0,10)*cA126 + B(1,10)*cA127 + B(2,10)*cA128 + B(3,10)*cA129 + B(4,10)*cA130 + B(5,10)*cA131;
rLeftHandSideMatrix(21,11)+=B(0,11)*cA126 + B(1,11)*cA127 + B(2,11)*cA128 + B(3,11)*cA129 + B(4,11)*cA130 + B(5,11)*cA131;
rLeftHandSideMatrix(21,12)+=B(0,12)*cA126 + B(1,12)*cA127 + B(2,12)*cA128 + B(3,12)*cA129 + B(4,12)*cA130 + B(5,12)*cA131;
rLeftHandSideMatrix(21,13)+=B(0,13)*cA126 + B(1,13)*cA127 + B(2,13)*cA128 + B(3,13)*cA129 + B(4,13)*cA130 + B(5,13)*cA131;
rLeftHandSideMatrix(21,14)+=B(0,14)*cA126 + B(1,14)*cA127 + B(2,14)*cA128 + B(3,14)*cA129 + B(4,14)*cA130 + B(5,14)*cA131;
rLeftHandSideMatrix(21,15)+=B(0,15)*cA126 + B(1,15)*cA127 + B(2,15)*cA128 + B(3,15)*cA129 + B(4,15)*cA130 + B(5,15)*cA131;
rLeftHandSideMatrix(21,16)+=B(0,16)*cA126 + B(1,16)*cA127 + B(2,16)*cA128 + B(3,16)*cA129 + B(4,16)*cA130 + B(5,16)*cA131;
rLeftHandSideMatrix(21,17)+=B(0,17)*cA126 + B(1,17)*cA127 + B(2,17)*cA128 + B(3,17)*cA129 + B(4,17)*cA130 + B(5,17)*cA131;
rLeftHandSideMatrix(21,18)+=B(0,18)*cA126 + B(1,18)*cA127 + B(2,18)*cA128 + B(3,18)*cA129 + B(4,18)*cA130 + B(5,18)*cA131;
rLeftHandSideMatrix(21,19)+=B(0,19)*cA126 + B(1,19)*cA127 + B(2,19)*cA128 + B(3,19)*cA129 + B(4,19)*cA130 + B(5,19)*cA131;
rLeftHandSideMatrix(21,20)+=B(0,20)*cA126 + B(1,20)*cA127 + B(2,20)*cA128 + B(3,20)*cA129 + B(4,20)*cA130 + B(5,20)*cA131;
rLeftHandSideMatrix(21,21)+=B(0,21)*cA126 + B(1,21)*cA127 + B(2,21)*cA128 + B(3,21)*cA129 + B(4,21)*cA130 + B(5,21)*cA131;
rLeftHandSideMatrix(21,22)+=B(0,22)*cA126 + B(1,22)*cA127 + B(2,22)*cA128 + B(3,22)*cA129 + B(4,22)*cA130 + B(5,22)*cA131;
rLeftHandSideMatrix(21,23)+=B(0,23)*cA126 + B(1,23)*cA127 + B(2,23)*cA128 + B(3,23)*cA129 + B(4,23)*cA130 + B(5,23)*cA131;
rLeftHandSideMatrix(22,0)+=B(0,0)*cA132 + B(1,0)*cA133 + B(2,0)*cA134 + B(3,0)*cA135 + B(4,0)*cA136 + B(5,0)*cA137;
rLeftHandSideMatrix(22,1)+=B(0,1)*cA132 + B(1,1)*cA133 + B(2,1)*cA134 + B(3,1)*cA135 + B(4,1)*cA136 + B(5,1)*cA137;
rLeftHandSideMatrix(22,2)+=B(0,2)*cA132 + B(1,2)*cA133 + B(2,2)*cA134 + B(3,2)*cA135 + B(4,2)*cA136 + B(5,2)*cA137;
rLeftHandSideMatrix(22,3)+=B(0,3)*cA132 + B(1,3)*cA133 + B(2,3)*cA134 + B(3,3)*cA135 + B(4,3)*cA136 + B(5,3)*cA137;
rLeftHandSideMatrix(22,4)+=B(0,4)*cA132 + B(1,4)*cA133 + B(2,4)*cA134 + B(3,4)*cA135 + B(4,4)*cA136 + B(5,4)*cA137;
rLeftHandSideMatrix(22,5)+=B(0,5)*cA132 + B(1,5)*cA133 + B(2,5)*cA134 + B(3,5)*cA135 + B(4,5)*cA136 + B(5,5)*cA137;
rLeftHandSideMatrix(22,6)+=B(0,6)*cA132 + B(1,6)*cA133 + B(2,6)*cA134 + B(3,6)*cA135 + B(4,6)*cA136 + B(5,6)*cA137;
rLeftHandSideMatrix(22,7)+=B(0,7)*cA132 + B(1,7)*cA133 + B(2,7)*cA134 + B(3,7)*cA135 + B(4,7)*cA136 + B(5,7)*cA137;
rLeftHandSideMatrix(22,8)+=B(0,8)*cA132 + B(1,8)*cA133 + B(2,8)*cA134 + B(3,8)*cA135 + B(4,8)*cA136 + B(5,8)*cA137;
rLeftHandSideMatrix(22,9)+=B(0,9)*cA132 + B(1,9)*cA133 + B(2,9)*cA134 + B(3,9)*cA135 + B(4,9)*cA136 + B(5,9)*cA137;
rLeftHandSideMatrix(22,10)+=B(0,10)*cA132 + B(1,10)*cA133 + B(2,10)*cA134 + B(3,10)*cA135 + B(4,10)*cA136 + B(5,10)*cA137;
rLeftHandSideMatrix(22,11)+=B(0,11)*cA132 + B(1,11)*cA133 + B(2,11)*cA134 + B(3,11)*cA135 + B(4,11)*cA136 + B(5,11)*cA137;
rLeftHandSideMatrix(22,12)+=B(0,12)*cA132 + B(1,12)*cA133 + B(2,12)*cA134 + B(3,12)*cA135 + B(4,12)*cA136 + B(5,12)*cA137;
rLeftHandSideMatrix(22,13)+=B(0,13)*cA132 + B(1,13)*cA133 + B(2,13)*cA134 + B(3,13)*cA135 + B(4,13)*cA136 + B(5,13)*cA137;
rLeftHandSideMatrix(22,14)+=B(0,14)*cA132 + B(1,14)*cA133 + B(2,14)*cA134 + B(3,14)*cA135 + B(4,14)*cA136 + B(5,14)*cA137;
rLeftHandSideMatrix(22,15)+=B(0,15)*cA132 + B(1,15)*cA133 + B(2,15)*cA134 + B(3,15)*cA135 + B(4,15)*cA136 + B(5,15)*cA137;
rLeftHandSideMatrix(22,16)+=B(0,16)*cA132 + B(1,16)*cA133 + B(2,16)*cA134 + B(3,16)*cA135 + B(4,16)*cA136 + B(5,16)*cA137;
rLeftHandSideMatrix(22,17)+=B(0,17)*cA132 + B(1,17)*cA133 + B(2,17)*cA134 + B(3,17)*cA135 + B(4,17)*cA136 + B(5,17)*cA137;
rLeftHandSideMatrix(22,18)+=B(0,18)*cA132 + B(1,18)*cA133 + B(2,18)*cA134 + B(3,18)*cA135 + B(4,18)*cA136 + B(5,18)*cA137;
rLeftHandSideMatrix(22,19)+=B(0,19)*cA132 + B(1,19)*cA133 + B(2,19)*cA134 + B(3,19)*cA135 + B(4,19)*cA136 + B(5,19)*cA137;
rLeftHandSideMatrix(22,20)+=B(0,20)*cA132 + B(1,20)*cA133 + B(2,20)*cA134 + B(3,20)*cA135 + B(4,20)*cA136 + B(5,20)*cA137;
rLeftHandSideMatrix(22,21)+=B(0,21)*cA132 + B(1,21)*cA133 + B(2,21)*cA134 + B(3,21)*cA135 + B(4,21)*cA136 + B(5,21)*cA137;
rLeftHandSideMatrix(22,22)+=B(0,22)*cA132 + B(1,22)*cA133 + B(2,22)*cA134 + B(3,22)*cA135 + B(4,22)*cA136 + B(5,22)*cA137;
rLeftHandSideMatrix(22,23)+=B(0,23)*cA132 + B(1,23)*cA133 + B(2,23)*cA134 + B(3,23)*cA135 + B(4,23)*cA136 + B(5,23)*cA137;
rLeftHandSideMatrix(23,0)+=B(0,0)*cA138 + B(1,0)*cA139 + B(2,0)*cA140 + B(3,0)*cA141 + B(4,0)*cA142 + B(5,0)*cA143;
rLeftHandSideMatrix(23,1)+=B(0,1)*cA138 + B(1,1)*cA139 + B(2,1)*cA140 + B(3,1)*cA141 + B(4,1)*cA142 + B(5,1)*cA143;
rLeftHandSideMatrix(23,2)+=B(0,2)*cA138 + B(1,2)*cA139 + B(2,2)*cA140 + B(3,2)*cA141 + B(4,2)*cA142 + B(5,2)*cA143;
rLeftHandSideMatrix(23,3)+=B(0,3)*cA138 + B(1,3)*cA139 + B(2,3)*cA140 + B(3,3)*cA141 + B(4,3)*cA142 + B(5,3)*cA143;
rLeftHandSideMatrix(23,4)+=B(0,4)*cA138 + B(1,4)*cA139 + B(2,4)*cA140 + B(3,4)*cA141 + B(4,4)*cA142 + B(5,4)*cA143;
rLeftHandSideMatrix(23,5)+=B(0,5)*cA138 + B(1,5)*cA139 + B(2,5)*cA140 + B(3,5)*cA141 + B(4,5)*cA142 + B(5,5)*cA143;
rLeftHandSideMatrix(23,6)+=B(0,6)*cA138 + B(1,6)*cA139 + B(2,6)*cA140 + B(3,6)*cA141 + B(4,6)*cA142 + B(5,6)*cA143;
rLeftHandSideMatrix(23,7)+=B(0,7)*cA138 + B(1,7)*cA139 + B(2,7)*cA140 + B(3,7)*cA141 + B(4,7)*cA142 + B(5,7)*cA143;
rLeftHandSideMatrix(23,8)+=B(0,8)*cA138 + B(1,8)*cA139 + B(2,8)*cA140 + B(3,8)*cA141 + B(4,8)*cA142 + B(5,8)*cA143;
rLeftHandSideMatrix(23,9)+=B(0,9)*cA138 + B(1,9)*cA139 + B(2,9)*cA140 + B(3,9)*cA141 + B(4,9)*cA142 + B(5,9)*cA143;
rLeftHandSideMatrix(23,10)+=B(0,10)*cA138 + B(1,10)*cA139 + B(2,10)*cA140 + B(3,10)*cA141 + B(4,10)*cA142 + B(5,10)*cA143;
rLeftHandSideMatrix(23,11)+=B(0,11)*cA138 + B(1,11)*cA139 + B(2,11)*cA140 + B(3,11)*cA141 + B(4,11)*cA142 + B(5,11)*cA143;
rLeftHandSideMatrix(23,12)+=B(0,12)*cA138 + B(1,12)*cA139 + B(2,12)*cA140 + B(3,12)*cA141 + B(4,12)*cA142 + B(5,12)*cA143;
rLeftHandSideMatrix(23,13)+=B(0,13)*cA138 + B(1,13)*cA139 + B(2,13)*cA140 + B(3,13)*cA141 + B(4,13)*cA142 + B(5,13)*cA143;
rLeftHandSideMatrix(23,14)+=B(0,14)*cA138 + B(1,14)*cA139 + B(2,14)*cA140 + B(3,14)*cA141 + B(4,14)*cA142 + B(5,14)*cA143;
rLeftHandSideMatrix(23,15)+=B(0,15)*cA138 + B(1,15)*cA139 + B(2,15)*cA140 + B(3,15)*cA141 + B(4,15)*cA142 + B(5,15)*cA143;
rLeftHandSideMatrix(23,16)+=B(0,16)*cA138 + B(1,16)*cA139 + B(2,16)*cA140 + B(3,16)*cA141 + B(4,16)*cA142 + B(5,16)*cA143;
rLeftHandSideMatrix(23,17)+=B(0,17)*cA138 + B(1,17)*cA139 + B(2,17)*cA140 + B(3,17)*cA141 + B(4,17)*cA142 + B(5,17)*cA143;
rLeftHandSideMatrix(23,18)+=B(0,18)*cA138 + B(1,18)*cA139 + B(2,18)*cA140 + B(3,18)*cA141 + B(4,18)*cA142 + B(5,18)*cA143;
rLeftHandSideMatrix(23,19)+=B(0,19)*cA138 + B(1,19)*cA139 + B(2,19)*cA140 + B(3,19)*cA141 + B(4,19)*cA142 + B(5,19)*cA143;
rLeftHandSideMatrix(23,20)+=B(0,20)*cA138 + B(1,20)*cA139 + B(2,20)*cA140 + B(3,20)*cA141 + B(4,20)*cA142 + B(5,20)*cA143;
rLeftHandSideMatrix(23,21)+=B(0,21)*cA138 + B(1,21)*cA139 + B(2,21)*cA140 + B(3,21)*cA141 + B(4,21)*cA142 + B(5,21)*cA143;
rLeftHandSideMatrix(23,22)+=B(0,22)*cA138 + B(1,22)*cA139 + B(2,22)*cA140 + B(3,22)*cA141 + B(4,22)*cA142 + B(5,22)*cA143;
rLeftHandSideMatrix(23,23)+=B(0,23)*cA138 + B(1,23)*cA139 + B(2,23)*cA140 + B(3,23)*cA141 + B(4,23)*cA142 + B(5,23)*cA143;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKg(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& DN_DX,
    const Vector& StressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix stress_tensor_x_weigth = IntegrationWeight * MathUtils<double>::StressVectorToTensor( StressVector );
    Matrix reduced_Kg(DN_DX.size1(), DN_DX.size1());
    MathUtils<double>::BDBtProductOperation(reduced_Kg, stress_tensor_x_weigth, DN_DX);
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, reduced_Kg, dimension );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    this->CalculateAndAddExtForceContribution( rThisKinematicVariables.N, rCurrentProcessInfo, rBodyForce, rRightHandSideVector, IntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rStressVector );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddExtForceContribution(
    const Vector& rN,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    VectorType& rRightHandSideVector,
    const double Weight
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const SizeType index = dimension * i;

        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[index + j] += Weight * rN[i] * rBodyForce[j];
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rLumpedMassVector.size() != mat_size)
        rLumpedMassVector.resize( mat_size, false );

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = GetGeometry().DomainSize() * density * thickness;

    Vector lumping_factors;
    lumping_factors = GetGeometry().LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < dimension; ++j ) {
            IndexType index = i * dimension + j;
            rLumpedMassVector[index] = temp;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateDampingMatrixWithLumpedMass(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
        beta = GetProperties()[RAYLEIGH_BETA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];

    // Compose the Damping Matrix:
    // Rayleigh Damping Matrix: alpha*M + beta*K

    // 2.-Calculate mass matrix:
    if (alpha > std::numeric_limits<double>::epsilon()) {
        VectorType temp_vector(mat_size);
        this->CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);
        for (IndexType i = 0; i < mat_size; ++i)
            rDampingMatrix(i, i) += alpha * temp_vector[i];
    }

    // 3.-Calculate StiffnessMatrix:
    if (beta > std::numeric_limits<double>::epsilon()) {
        MatrixType stiffness_matrix( mat_size, mat_size );
        VectorType residual_vector( mat_size );

        this->CalculateAll(stiffness_matrix, residual_vector, rCurrentProcessInfo, true, false);

        noalias( rDampingMatrix ) += beta  * stiffness_matrix;
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters BaseSolidElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","implicit","explicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["INTEGRATION_WEIGHT","STRAIN_ENERGY","ERROR_INTEGRATION_POINT","VON_MISES_STRESS","INSITU_STRESS","CAUCHY_STRESS_VECTOR","PK2_STRESS_VECTOR","GREEN_LAGRANGE_STRAIN_VECTOR","ALMANSI_STRAIN_VECTOR","CAUCHY_STRESS_TENSOR","PK2_STRESS_TENSOR","GREEN_LAGRANGE_STRAIN_TENSOR","ALMANSI_STRAIN_TENSOR","CONSTITUTIVE_MATRIX","DEFORMATION_GRADIENT","CONSTITUTIVE_LAW"],
            "nodal_historical"       : ["DISPLACEMENT","VELOCITY","ACCELERATION"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Triangle2D6", "Quadrilateral2D4", "Quadrilateral2D8", "Quadrilateral2D9","Tetrahedra3D4", "Prism3D6", "Prism3D15", "Hexahedra3D8", "Hexahedra3D20", "Hexahedra3D27", "Tetrahedra3D10"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : -1,
        "documentation"   : "This is a pure displacement element"
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

bool BaseSolidElement::IsElementRotated() const
{
    if (mConstitutiveLawVector[0]->GetStrainSize() == 6) {
        return (this->Has(LOCAL_AXIS_1) && this->Has(LOCAL_AXIS_2));
    } else if (mConstitutiveLawVector[0]->GetStrainSize() == 3) {
        return (this->Has(LOCAL_AXIS_1));
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}
} // Namespace Kratos
