


	Vector v(6); // or v.resize(6, false);
	Matrix A(6,6);

	noalias(A) = ...; // always set the size
	noalias(v) = ...;


	noalias(v) = prod(A, v); // NO!!
	noalias(v) += ...; // OK!!