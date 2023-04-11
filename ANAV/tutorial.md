
-- Prestressed concrete in Kratos --

1. Create geometry in GiD, volumes of concrete and lines for steel tendons.
	Obs: * Volumes must be in one layer and lines in another one.
	     * We identify each Tendon by creating a Group in GiD named Tendon_Whateveryouwant
2. Calculate intersections of lines with hexas by means of the Asco_PT02.gid problemtype. It generats a .intersections file
3. At the beginning of each tendon intersection block we add:       D = 1.121e-2    Ep = 6.858e-3
    which defines the diameter of the tendon and its imposed strain.
	Example:
	    Begin Tendon - Hexahedra intersection: Tendon_H1      D = 1.121e-2    Ep = 6.858e-3
	       227	        0.24066	    -4998.91279	        0.05066		        0.24066	    -4998.81395	        0.05066
	       228	        0.24066	    -4998.81395	        0.05066		        0.24066	    -4998.71512	        0.05066
		   ...
		End Tendon - Hexahedra intersection
4. We define the material props of concrete, boundary conditions, loads, etc with the Kratos problemtype for GiD.
    * This will generate several files:
	        - .mdpa: defines the coordinates of nodes, connectivities of elements, submodelparts...
			- ProjectParameters.json: defines problem general parameters: time step, load values, imposed displacement values, printing...
			- StructuralMaterials.json: defines the name of the constitutive laws to be used and its material properties.
			- MainKratos.py: the main file to be run --> python MainKratos.py
5. Now we add the custom prestressing info in the standard Kratos files:
    * Create a new submodel part in the .mdpa file that includes all the elements intersected for each tendon (we can create one submodel per tendon or one submodel for several tendons)
	    Example:
				Begin SubModelPart Tendon_H1
				  Begin SubModelPartTables
				  End SubModelPartTables
				  Begin SubModelPartNodes
				  End SubModelPartNodes
				  Begin SubModelPartElements
					227     228     229     230     231     232     233     234     235     236
					237     238     239     240     241     242     243     244     245     246
					247     248     249     250     251     252     253     254     255     256
					3043    3045    3047    3049    3051    3053    3055    3057    3059    3061
					3063    3065    3067    3069    3071    3073    3075    3077    3079    3081
					3083    3085    3087    3089    3091    3093    3095    3657    3659    3661
					3663    3665    3667    3669    3671    3021   
				  End SubModelPartElements
				  Begin SubModelPartConditions
				  End SubModelPartConditions
				End SubModelPart
	
	* Include a block in the ProjectParameters.json in the "list_of_other_processes" to compute the % participation and orientation of the steel tendons.
	Example:
	    {
            "python_module" : "set_up_pre_stressed_oriented_composite_materials",
            "kratos_module" : "KratosMultiphysics.ConstitutiveLawsApplication",
            "process_name"  : "SetUpPreStressedOrientedCompositeMaterials",
            "Parameters"    : {
                "model_part_name"        : "Structure",
				"intersection_file_name" : "geometry.intersections",  // The file created and modified in steps 2 and 3
				"echo_level"             : 1
            }
        }
	
	* Add a new material property in the StructuralMaterials.json that defines the behaviour of the composite prestressed concrete. We must use the Serial-Parallel Rule of Mixtures.
	Example:
		{
			"model_part_name" : "Structure.tendon_H", // ok
			"properties_id"   : 8,
			"Material"        : {
				"constitutive_law" : {
					"name" : "SerialParallelRuleOfMixturesLaw",
					"combination_factors"          : [0.99999, 0.00001],  // dummy % participation of concrete and steel, will be automatically computed
					"parallel_behaviour_directions" : [1,0,0,0,0,0]       // do not touch
				},
				"Variables"        : {
					"DENSITY"       : 2000.0,
					"TANGENT_OPERATOR_ESTIMATION" : 2,
					"SYMMETRIZE_TANGENT_OPERATOR" : false
				},
				"Tables"           : {}
			},
			"sub_properties" : [{
				"properties_id"   : 81,
				"Material"        : {
					"constitutive_law" : {
						"name" : "LinearElastic3DLaw"  // concrete
				},
				"Variables"        : {
					"YOUNG_MODULUS" : 4.1049E+10,
					"POISSON_RATIO" : 0.2
					},
					"Tables"           : {}
				}
			},{
				"properties_id"   : 82,
				"Material"        : {
					"constitutive_law" : {
						"name" : "LinearElastic3DLaw"  // prestressing steel
				},
				"Variables"        : {
					"YOUNG_MODULUS" : 195150e6,
					"POISSON_RATIO" : 0.3
					},
					"Tables"           : {}
				}
			}]
		}