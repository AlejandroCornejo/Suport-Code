	std::vector<double> result;
	std::vector<double> damages(8);
	r_process_info = mrMopdelPart.GetProcessInfo();
	for elem in modelpart:
		elem.CalculateOnIntegrationPoints(THRESHOLD, result,r_process_info); // threshold vector (8 components for 3d hexas)
		// implement a way to get damage for this element (...)  --> d
		
		// setting damage to all IP
		for comp in damages:
			damages[comp] = d;
		elem.SetValuesOnIntegrationPoints(DAMAGE, damages, r_process_info);
		
		// setting threhsold to all IP
		for comp in damage:
			result[comp] /= (1-d)
		elem.SetValuesOnIntegrationPoints(THRESHOLD, result, r_process_info);
		
		
		
		
		
		