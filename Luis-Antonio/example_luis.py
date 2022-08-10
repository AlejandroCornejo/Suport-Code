
import KratosMultiphysics as KM
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility

# These are our settings...
our_settings = KM.Parameters("""{
            "name"             : "csv_table",
            "filename"         : "sample.csv",
            "delimiter"        : ",",
            "skiprows"         : 1,
            "first_column_id"  : 0,
            "second_column_id" : 1,
            "table_id"         : 0,
            "na_replace"       : 0.0
        }""")



print(our_settings["filename"].GetString()) # sample.csv
print(our_settings["filename"].GetString().split(".")[0]) #sample
our_settings["filename"].SetString("sampleLuis.csv")
print(our_settings["filename"].GetString())

a



model = KM.Model()
sub_model = model.CreateModelPart("trial_model")
print(model.HasModelPart("trial_model"))


# Here we are creating a Table by reading a CSV
table_sample = ReadCsvTableUtility(our_settings).Read(sub_model) # By default the model_part is a null one

#print(sub_model)

# Let's print some stuff to be sure everything is ok...
print(table_sample) # The whole table
print("---------------")
print("How much is the value interpolated of 3.77? --> ", table_sample.GetValue(3.77))
