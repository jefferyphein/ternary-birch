// Magma code for computing Hecke eigenvalues associated to rational
// eigenforms.

function readGenusFile(N, Ts)
	FP := Open("./data/genus/" cat IntegerToString(N) cat ".txt", "r");

	// Determine the genus size, then skip that many lines.
	temp := Gets(FP);
	genusSize := StringToInteger(temp);
	for n in [1..genusSize] do
		temp := Gets(FP);
	end for;

	blank := Gets(FP);
	temp := Gets(FP);
	while temp ne "" do
		array := Split(temp, " ");
		d := StringToInteger(array[1]);
		p := StringToInteger(array[2]);
		lines := StringToInteger(array[3]);
		dim := StringToInteger(array[4]);

		if not IsDefined(Ts, d) then
			Ts[d] := AssociativeArray();
		end if;

		// Does this matrix already exist in our data structure?
		found := IsDefined(Ts[d], p);

		// If not, initialize a new matrix.
		if not found then
			M := SparseMatrix(Integers(), dim, dim);
		end if;

		for n in [1..lines] do
			temp := Gets(FP);

			// If the matrix hasn't been computed yet, populate it.
			if not found then
				array := Split(temp, " ");
				row := StringToInteger(array[1]);
				entries := StringToInteger(array[2]);

				for m in [1..entries] do
					col := StringToInteger(array[2*m+1]);
					value := StringToInteger(array[2*m+2]);
					M[col+1][row+1] := value;
				end for;
			end if;
		end for;

		if not found then
			Ts[d][p] := M;
		end if;

		blank := Gets(FP);
		blank := Gets(FP);
		temp := Gets(FP);
	end while;

	delete FP;

	return Ts;
end function;

function readEigenforms(N)
	filename := "./data/vecs/" cat IntegerToString(N) cat ".txt";
	forms := [];
	lines := Split(Read(filename));
	for n in [1..(#lines - 1) div 2] do
		vecline := lines[2*n];
		eigs := AssociativeArray();
		values := Split(lines[1+2*n], " ");
		num := StringToInteger(values[1]);
		for n in [1..num] do
			p := StringToInteger(values[2*n]);
			ap := StringToInteger(values[1+2*n]);
			eigs[p] := ap;
		end for;
		Append(~forms, < vecline, eigs >);
	end for;
	return forms;
//	basename := IntegerToString(N) cat ".txt";
//	filename := "./data/vecs/" cat basename;
//
//	forms := [];
//
//	FP := Open(filename, "r");
//	condline := Gets(FP);
//	blank := Gets(FP);
//	vecline := Gets(FP);
//	while not IsEof(vecline) do
//		temp := Split(vecline, " ");
//		cond := StringToInteger(temp[1]);
//
//		eigs := AssociativeArray();
//
//		line := Gets(FP);
//		token := Split(line, " ");
//		if #token ne 0 then
//			num := StringToInteger(token[1]);
//			for n in [1..num] do
//				p := StringToInteger(token[2*n]);
//				ap := StringToInteger(token[2*n+1]);
//				eigs[p] := ap;
//			end for;
//
//			Append(~forms, < vecline, eigs >);
//		end if;
//
//		blank := Gets(FP);
//
//		vecline := Gets(FP);
//	end while;
//
//	return forms;
end function;

function getQuadForm(N)
	// Make sure our paperwork is in order.
	assert N gt 0;
	assert IsSquarefree(N);

	K := QNF();
	R := Integers(K);

	// If N is divisible by an even number of primes, ensure that the largest
	// such prime ramifies, this will be useful when we distinguish eigenvectors
	// based on whethey they arise from a lower level or not.
	primes := PrimeDivisors(N);
	max := #primes mod 2 eq 0 select Max(PrimeDivisors(N)) else -1;

	primes := { p*R : p in PrimeDivisors(N) | p mod 4 eq 1 xor p eq max };
	if #primes mod 2 eq 1 then
		primes join:= { 2*R };
	end if;

	F := QuadraticFormWithInvariants(3, K!N, primes, [0]);
	L := MaximalIntegralLattice(F);
	G := 2 * ChangeRing(GramMatrix(L), Rationals());
	L := LatticeWithGram(G);
	G := GramMatrix(LLL(L));

	assert Determinant(G) div 2 eq N;

	return [ G[1][1] div 2, G[2][2] div 2, G[3][3] div 2,
		     G[2][3], G[1][3], G[1][2] ];
end function;

procedure stepOne(N)
	assert N gt 0;
	assert IsSquarefree(N);

	// Determine the appropriate quadratic form for this discriminant.
	coeffs := getQuadForm(N);

	// Determine the first good prime for this quadratic form.
	p := 1;
	repeat p := NextPrime(p);
	until N mod p ne 0;

	// The filename to save/load to/from.
	filename := "./data/genus/" cat IntegerToString(N) cat ".txt";

	// Beginning of the command-line.
	commandLine := "./birch -a -x -o " cat filename cat " -p " cat IntegerToString(p);

	// If this file already exists, we're going to use it as an input file
	// instead of creating it.
	try
		ptr := Open(filename, "r");
		ptr := []; // Close file pointer.
		commandLine cat:= " -i " cat filename;
	catch e
		commandLine cat:= &cat[ " " cat IntegerToString(x) : x in coeffs ];
	end try;
	
	// Compute the Hecke operators.
	commandLine;
	System(commandLine);
end procedure;

procedure stepTwo(N)
	// The file basename.
	basename := IntegerToString(N) cat ".txt";

	// Input/output filenames.
	vecfilename := "./data/vecs/" cat basename;
	genfilename := "./data/genus/" cat basename;

	// Have we already computed eigenvectors? If so, we're done.
	try
		ptr := Open("./data/vecs/" cat basename, "r");
		ptr := []; // Close the file pointer.
		return;
	catch e
		;
	end try;

	// The largest prime for which we've computed Hecke operators.
	largest := 0;


	function search(Ts, d, pos, basis, eigs, vecs, output, Ms)
		if pos gt #Ts then
			ps := [ data[1] : data in Ts ];

			// Get the next good prime.
			p := Max(ps);
			repeat p := NextPrime(p);
			until N mod p ne 0;

			if not IsDefined(Ms[d], p) then
				// Need to compute more Hecke operators.
				commandLine := "./birch -a -x -o " cat genfilename cat
					           " -i " cat genfilename cat
							   " -p " cat IntegerToString(p);
				commandLine;
				System(commandLine);

				// Read genus file to obtain the newly computed Hecke operator.
				Ms := readGenusFile(N, Ms);
			end if;
			Append(~Ts, < p, Ms[d][p] >);
		end if;

		// Determine the prime and the associated Hecke operator.
		p := Ts[pos][1];
		T := Ts[pos][2];

		// The bound on the (absolute) size of the integral eigenvalue.
		x := Floor(2 * Sqrt(p));

		// Project the Hecke matrix onto the image of the basis.
		if pos ne 1 then
			T := Solution(basis, basis * T);
		end if;

		// Loop over all valid eigenvalues.
		for e in [-x..x] do
			// Substract off e*Id from the matrix so we can compute the eigenspace.
			M := T;
			for n in [1..Nrows(T)] do
				M[n][n] -:= e;
			end for;

			// The kernel is the eigenspace.
			ker := Kernel(M);

			// The dimension of the eigenspace.
			dim := Dimension(ker);

			// If the dimension is 1, we have an eigenvector.
			if dim eq 1 then
				// Update the data structure.
				temp := Append(eigs, < p, e >);

				// Obtain the eigenvector in the original global coordinates.
				if pos eq 1 then
					vec := Basis(ker)[1];
				else
					vec := Vector(Matrix(Basis(ker)) * basis);
				end if;

				// Build the output string to be written to file.
				string := IntegerToString(d);
				string cat:= &cat[ " " cat IntegerToString(x) : x in Eltseq(vec) ];

				// Write to file.
				output cat:= string cat "\n";
				output cat:= "0\n";
				output cat:= "\n";
				Append(~vecs, temp);
			elif dim gt 1 then
				// Update the data structure.
				temp := Append(eigs, < p, e >);

				// Convert the basis of the eigenspace into a matrix.
				B := Matrix(Basis(ker));

				// Project the basis into the original coordinates.
				if pos ne 1 then B *:= basis; end if;

				// Continue searching via the next Hecke operator.
				vecs, Ms, output := search(Ts, d, pos+1, B, temp, vecs, output, Ms);
			end if;
		end for;

		return vecs, Ms, output;
	end function;

	// The eigenvector data file does not exist, so we need to create it.
	mats := readGenusFile(N, AssociativeArray());

	// Compute rational eigenvectors and generate output to be written to file.
	output := "";
	found := [];
	keys1 := Keys(mats);
	for d in keys1 do
		keys2 := Keys(mats[d]);
		Ts := [];
		for p in keys2 do
			Append(~Ts, < p, mats[d][p] >);
		end for;
		vecs, mats, output := search(Ts, d, 1, [], [], [], output, mats);
		Append(~found, vecs);
	end for;

	// Determine the distinct conductors which have rational eigenvectors.
	keys := [ x : x in keys1 ];
	conds := [];
	for n in [1..#keys] do
		if #found[n] gt 0 then
			conds cat:= [keys[n]];
		end if;
	end for;
	firstLine := IntegerToString(#conds) cat &cat[ " " cat IntegerToString(x) : x in conds ];

	// Write output eigenvector data file.
	OP := Open(vecfilename, "w");
	Puts(OP, firstLine);
	Puts(OP, "");
	Puts(OP, output);
	delete OP;
end procedure;

procedure stepThree(N : upTo := 1000)
	assert N gt 0;
	assert IsSquarefree(N);
	assert upTo ge 1;

	basename := IntegerToString(N) cat ".txt";
	vecfilename := "./data/vecs/" cat basename;
	genfilename := "./data/genus/" cat basename;

	// Compute eigenvalues.
	commandLine := "./birch -x -u " cat IntegerToString(upTo) cat " -i " cat genfilename cat
				   " -e " cat vecfilename;
	commandLine;
	System(commandLine);

	// If the number of prime divisors is odd, then we're done, all rational
	// eigenvectors of level N have been found.
	if #PrimeDivisors(N) mod 2 eq 1 then return; end if;

	// Otherwise, we need to rule out eigenforms arising from lower levels.
	// We chose the original quadratic form in such a way that its largest prime
	// ramified, and so we need only remove it and check if the eigenvalues
	// appear at a lower level.
	D := N div Max(PrimeDivisors(N));

	// Compute all eigenvectors at the lower level. All primes will be unramified
	// and so all eigenvectors found here will be associated to elliptic curves
	// with conductor D.
	birch(D : StepFour := false);

	// Load the eigenvectors into memory.
	forms_D := readEigenforms(D);
	forms_N := readEigenforms(N);

	// If no forms found at larger level, we're done.
	if #forms_N eq 0 then return; end if;

	// If no forms found at smaller level, we're done.
	if #forms_D eq 0 then return; end if;

	// A list of primes for which to compare Hecke eigenvalues.
	keys := Keys(forms_N[1][2]);
	for n in [2..#forms_N] do keys meet:= Keys(forms_N[n][2]); end for;
	for d in [1..#forms_D] do keys meet:= Keys(forms_D[d][2]); end for;
	keys := { x : x in keys | N mod x ne 0 };

	// TODO: Sturm bounds should probably go here.
	assert #keys gt 50;

	conductors := {};

	n := 1;
	repeat
		// Did we find a match?
		found := false;

		for d in [1..#forms_D] do
			// Initialize match to true to get things started.
			match := true;

			// Compare Hecke eigenvalues to see if we have a match.
			for key in keys do
				match := match and forms_D[d][2][key] eq forms_N[n][2][key];
			end for;

			if match then
				// If a match is found at a lower level, remove it from the list as
				// it has level D.
				Remove(~forms_N, n);

				// Indicate that we found a match at a lower level.
				found := true;

				// And then break out of the for-loop as no more comparisons needed.
				break;
			end if;
		end for;

		// If no match found, the current eigenform has level N.
		if not found then
			conductors join:= { Split(forms_N[n][1], " ")[1] };
			n +:= 1;
		end if;
	until n gt #forms_N;

	// Open a file stream for writing.
	eigfilename := "./data/vecs/" cat IntegerToString(N) cat ".txt";
	OP := Open(eigfilename, "w");

	// All distinct character conductors.
	conductors := Sort([ x : x in conductors ]);

	string := IntegerToString(#conductors) cat &cat[ " " cat x : x in conductors ];
	Puts(OP, string);
	Puts(OP, "");

	for form in forms_N do
		// Organize the keys in increasing order.
		keys := Keys(form[2]);
		keys := Sort([ x : x in keys ]);

		// Build the string containing all Hecke eigenvalue data.
		string := IntegerToString(#keys);
		string cat:= &cat[ " " cat IntegerToString(p) cat " " cat
						   IntegerToString(form[2][p]) : p in keys ];

		Puts(OP, form[1]);
		Puts(OP, string);
		Puts(OP, "");
	end for;

	delete OP;
end procedure;

procedure stepFour(N)
	QQ := RationalsAsNumberField();
	ZZ := Integers(QQ);
	filename := "./data/vecs/" cat IntegerToString(N) cat ".txt";
	lines := Split(Read(filename));
	for n in [1..(#lines-1) div 2] do
		values := [ StringToInteger(x) : x in Split(lines[1+2*n], " ") ];
		primes := [ values[n]*ZZ : n in [1..#values] | n mod 2 eq 0 and N mod values[n] ne 0 ];
		traces := [ values[n] : n in [3..#values] | n mod 2 eq 1 and N mod values[n-1] ne 0 ];
		Es := EllipticCurveSearch(N*ZZ, 2000 : Primes := primes, Traces := traces);
		for E in Es do
			Norm(Conductor(E)), E;
		end for;
	end for;
end procedure;

intrinsic birch(N::RngIntElt : StepFour := true) {}
	// This code (currently) only works for squarefree discriminants.
	if N lt 0 then N := -N; end if;
	require IsSquarefree(N): "Discriminant must be squarefree.";

	stepOne(N);
	stepTwo(N);
	stepThree(N);
	if StepFour then
		stepFour(N);
	end if;
end intrinsic;
