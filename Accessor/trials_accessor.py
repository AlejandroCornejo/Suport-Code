
import KratosMultiphysics as KM

# '''
{
    "model_part_name" : "Parts_exterior",
    "properties_id"   : 2,
    "Material"        : {
        "constitutive_law" : {
            "name" : "LinearElastic3DLaw"
        },
        "Variables"        : {
            "DENSITY"       : 7850.0,
            "YOUNG_MODULUS" : 2.1e11,
            "POISSON_RATIO" : 0.29
        },
        "Tables"           : {
            "TEMPERATURE_vs_YOUNG_MODULUS" : { // we have a name for each table
                "input_variable"  : "TEMPERATURE",
                "output_variable" : "YOUNG_MODULUS",
                "data"            : [[-1.0,   2.1e11],
                                     [1000,   2.1e11],
                                     [3500, 2.1e8 ],
                                     [1e6 ,  2.1e8 ]]
            }
        },
        "Accessors"        : {
            "accessor_table_T_E" : {
                "accessor_type"  : "table_accessor", // here we indicate the table accessor
                "table_input_variable"  : "TEMPERATURE",
                "table_output_variable" : "YOUNG_MODULUS"
            }
        }
    }
}
#'''

