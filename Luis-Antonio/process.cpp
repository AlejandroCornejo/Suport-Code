void SetAutomatedInitialStateProcess::ExecuteInitialize()
{
    KRATOS_TRY
    array_1d<double, 6> strain_vector;
    array_1d<double, 6> initial_stress_vector;
    // const array_1d<double, 3> element_centroid;
    double element_radial_cordinate;

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        // int TableId = mThisParameters["initial_state_table"];
        // KRATOS_WATCH (TableId)
        int ElemId = rElement.Id()
        KRATOS_WATCH (ElemId)
        const array_1d<double, 3> element_centroid = rElement.GetGeometry().Center();
        element_radial_cordinate = sqrt(element_centroid[0] * element_centroid[0] + element_centroid[1] * element_centroid[1]);
        initial_stress_vector[0] = mrThisModelPart.GetTable(0).GetValue(element_radial_cordinate);
        initial_stress_vector[1] = 0;
        initial_stress_vector[2] = 0;
        initial_stress_vector[3] = 0;
        initial_stress_vector[4] = 0;
        initial_stress_vector[5] = 0;
        rElement.SetValue(INITIAL_STRESS_VECTOR, initial_stress_vector);
        KRATOS_WATCH(element_centroid)
        KRATOS_WATCH(initial_stress_vector)
    });
    KRATOS_CATCH("")
}