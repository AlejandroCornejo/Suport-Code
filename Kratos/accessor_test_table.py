import KratosMultiphysics as KM

current_model = KM.Model()
model_part = current_model.CreateModelPart("Main")
model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)

node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

submodel = model_part.CreateSubModelPart("solid")

geom = KM.Triangle2D3(node1,node2,node3)

N = KM.Vector(3)

# properties_1 = model_part.CreateNewProperties(1)
test_settings = KM.Parameters("""
    {
        "properties" : [{
            "model_part_name" : "Structure.solid",
            "properties_id"   : 1,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "LinearElasticPlaneStrain2DLaw"
                },
                "Variables"        : {
                    "DENSITY"       : 7850.0,
                    "YOUNG_MODULUS" : 1.5e11,
                    "POISSON_RATIO" : 0.29
                },
                "Tables"           : {
                    "TEMPERATURE_vs_YOUNG_MODULUS" : {
                        "input_variable"  : "TEMPERATURE",
                        "output_variable" : "YOUNG_MODULUS",
                        "data"            : [[-1.0,   2.1e11],
                                            [1000,    2.1e11],
                                            [3500,    2.1e11],
                                            [1e6 ,    2.1e11]]
                    }
                },
                "Accessors"        : {
                    "accessor_table_T_E" : {
                        "accessor_type"             : "table_accessor",
                        "table_input_variable"      : "TEMPERATURE",
                        "table_output_variable"     : "YOUNG_MODULUS",
                        "table_input_variable_type" : "nodal_historical"
                    }
                }
            }
        }]
    }
""")

KM.ReadMaterialsUtility(test_settings, current_model)

