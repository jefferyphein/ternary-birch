//freeze;
//SetDebugOnError(true);


/******************************************************************
*                                                                 *
* Spinor genera of lattices over number fields                    *
* BONGs of lattices over number fields                            *
* Genera of lattices over number fields                           *
*                                                                 *
* See D. Lorch, "Einklassige Geschlechter orthogonaler Gruppen",  * 
*                                                                 *
* David Lorch, 2015                                               *
*                                                                 *
******************************************************************/

// we fix some primes for the genus computations and add them as a global property of the field / lattice:

declare attributes LatMod:
  GenusComp;
  
GenusComputationFormat:=recformat<
  RCG:GrpAb, // a ray class group associated to some lattice L
  hom_RCG:Map, 

  MM:SeqEnum, // list of prime powers defining RCG
  M:RngOrdIdl, // M = &*MM
  RayPrimes, // like MM, but is a list of the prime ideals only (no exponents)

  RCG_Subgroup:GrpAb, // a subgroup of RCG
  RCG_FactorGroup:GrpAb, // = RCG modulo RCG_Subgroup
  hom_RCG_FactorGroup:Map, 

  CriticalPrimes:SeqEnum, // a set of primes needed to list all the proper spinor genera in Genus(L)
  GroupElementsToPrimes: Assoc, // maps the elements of RCG_Factorgroup to primes representing them
  
  Status: RngIntElt // 0 = not initialized; 1 = RCG initialized, 2 = RCG & factor group initialized, 3 = fully initialized
>; 

intrinsic deletedata(L::LatMod)
  { deletes genus data stored in L }
  if assigned L`GenusComp then delete L`GenusComp; end if;
  L`GenusComp := rec< GenusComputationFormat | Status:=0, GroupElementsToPrimes:=AssociativeArray()>;
  return;
end intrinsic;
  
// 0. Helper functions, and some basic things which might later be turned into intrinsics:

function GetExtra(L, Primes, Extra)
  // return more primes, the number of which is specified by Extra, which are coprime to Primes and BadPrimes(L).
  assert Extra ne 0;
  K:= NumberField(BaseRing(L));
  Coprime:= &*(Set(Primes) join BadPrimes(L));

  Bnd:= 30;
  repeat
    Ps:= PrimesUpTo(Bnd, K : coprime_to:= 2*Coprime);
    Bnd +:= 10;
  until #Ps ge Extra;
  return Ps[1..Extra];
end function;

// TODO: replace the Jordan splitting in lat.m by this one.
function DetailedJordanSplitting(M, F, p)
  O:= Order(p);
  even:= Minimum(p) eq 2;
  if even then pi:= UniformizingElement(p); end if;
//  B:= LocalBasis(M, p: Type:= "Submodule");
  B:= LocalBasis(M, p);
  n:= Ncols(F);
  k:= 1;

  oldval:= Infinity();
  Blocks:= []; Exponents:= [];
  Steps_temp := [];
  S:= Matrix(B);
  Steps := []; // list of numbers of base vectors forming the 1x1 or 2x2 matrices that comprise the Jordan decomposition
  while k le n do
    G:= S*F*Transpose(S);

    // find an element <i, j> with minimal p-valuation.
    // Take it from the diagonal, if possible, and from the lowest-possible row number.
    X:= [ Valuation(G[i,i], p): i in [k..n] ];
    m, ii:= Minimum( X ); ii +:= k-1;
    pair:= <ii, ii>;

    for i in [k..n], j in [i+1..n] do
      tmp:= Valuation(G[i,j], p);
      if tmp lt m then m:= tmp; pair:= <i,j>; end if;
    end for;
    if m ne oldval then Append(~Blocks, k); oldval:= m; Append(~Exponents, m); end if;

    if even and pair[1] ne pair[2] then
      SwapRows(~S, pair[1],   k); // swap f_1 and e_i
      SwapRows(~S, pair[2], k+1); // swap f_2 and e_j

      T12:= (S[k] * F * Matrix(1, Eltseq(S[k+1])))[1];
      S[k] *:= pi^Valuation(T12, p)/T12;
      T := func<i,j|(S[i] * F * Matrix(1, Eltseq(S[j])))[1]>;
      T11 := T(k,k); T22 := T(k+1, k+1); T12:= T(k, k+1);

      d := T11*T22 - T12^2;
      for l in [k+2..n] do
        tl := T12*T(k+1,l)-T22*T(k  ,l); // t_k from step 4
        ul := T12*T(k  ,l)-T11*T(k+1,l); // u_k from step 4
        S[l] +:= (tl/d)*S[k] + (ul/d)*S[k+1];
      end for;
      k +:= 2;
      Append(~Steps_temp, 2);
    else
      if pair[1] eq pair[2] then
        SwapRows(~S, pair[1], k);
      else
        S[pair[1]] +:= S[pair[2]];
        SwapRows(~S, pair[1], k);
      end if;
      nrm:= (S[k] * F * Matrix(1, Eltseq(S[k])))[1];
      X:= S[k] * F * Transpose(S); // T(e_k, f_i), 1<=i<=n
      for l in [k+1..n] do S[l] -:= X[l]/nrm * S[k]; end for;
      k+:= 1;
      Append(~Steps_temp, 1);
    end if;
  end while;

  Append(~Blocks, n+1);
  Matrices:= [* RowSubmatrix(S, Blocks[i], Blocks[i+1] - Blocks[i]) : i in [1..#Blocks-1] *];
  
  Steps := [];
  j := 1;
  for i in [1..#Blocks-1] do 
    // return a finer decomposition of the blocks into 1x1 and 2x2 matrices:
    start := Blocks[i];
    len := Blocks[i+1]-Blocks[i];
    while len gt 0 do
      if Steps_temp[j] eq 2 then Append(~Steps, [start, start+1]); len -:= 2; start +:= 2; else Append(~Steps, [start]); len -:= 1; start +:= 1; end if;
      j +:= 1;
    end while;
  end for;
  return Matrices, [* m*F*Transpose(m): m in Matrices *], Exponents, Blocks, Steps;
end function;

intrinsic DetailedJordanDecomposition(L::LatMod, p::RngOrdIdl) -> []
{ A Jordan decomposition of the completion of L at p -- more detailed output than JordanDecomposition }
  return DetailedJordanSplitting(L`Module, L`Form, p); 
end intrinsic;



function ScalesAndNorms(G, p, Uniformizer)
  // G is a splitting of a LatticeModule L into Gram matrices,
  // which can e.g. be calculated by JordanDecomposition(), but need not
  // be a Jordan decomposition.

  // We calculate the invariants of the entries of G, viewed as individual lattices.
  
  e := RamificationIndex(p);
  sL := [Minimum({Valuation(g[i,j], p): i in [1..j], j in [1..Ncols(g)]}): g in G];

  aL:= []; uL:= []; wL:= [];
  sL := [Minimum({Valuation(x, p): x in Eltseq(g)}): g in G];
  for i in [1..#G] do
    // similar, but not equal to, the code for obtaining a genus symbol
    // (for the genus symbol, we consider a scaling L^(s(L_i)) in O'Meara's notation)
    GG := G[i];

    D:= Diagonal(GG);
    // 2*s(L_i) eq n_i?
    if e+sL[i] le Minimum({ Valuation(d, p) : d in D }) then
      Append(~aL, FieldOfFractions(Order(p))!Uniformizer^(e+sL[i]));
    else
      Append(~aL, D[b] where _, b is Minimum([Valuation(x,p): x in D]));
    end if;
    Append(~uL, Valuation(aL[i], p));
    assert uL[i] eq Valuation(Norm(LatticeModule(G[i])), p);
    Append(~wL, Minimum({e+sL[i]} join {uL[i] + QuadraticDefect(d/aL[i], p): d in D}));
  end for;

  return sL, aL, uL, wL; // scales, norm generators, norm valuations, weight valuations of a (Jordan) decomposition of L
end function;


function NormUpscaled(G, p)
  // calculate Norm(L^{scale(L_i)}), for all i.
  // for the definition of L^{a}, cf. § 82 I, p.237, in O'Meara
  
  // G is a splitting of a LatticeModule L into Gram matrices,
  // which can e.g. be calculated by JordanDecomposition(), but need not
  // be a Jordan decomposition.
    t:= #G;
    sL := [Minimum({Valuation(g[i,j], p): i in [1..j], j in [1..Ncols(g)]}): g in G];
    e := RamificationIndex(p);
    Uniformizer := UniformizingElement(p);
    
    aL:= []; uL:= []; wL:= [];
    for i in [1..t] do
      // Die Norm ist 2*Scale + <die Ideale die von den Diagonalelementen erzeugt werden>, cf. § 94 O'Meara.
      GG:= DiagonalJoin(< j lt i select Uniformizer^(2*(sL[i]-sL[j])) * G[j] else G[j] : j in [1..t] >);
      D:= Diagonal(GG);
      // Is 2*s(L_i) eq n_i? (We always have scale \supseteq norm \supseteq 2*scale)
      vprintf Genus, 2: "val(2*scale) = %o, val(diagonal)=%o\n", e+sL[i], Minimum({ Valuation(d, p) : d in D }); 
      // um einen Normerzeuger zu finden, muss man einfach Erzeuger des größeren der Ideale 2*Scale bzw. <Diagonalelemente> nehmen:
      if e+sL[i] le Minimum({ Valuation(d, p) : d in D }) then
        // 20150505: note that "le", not "eq", is correct here (we are not comparing 2*s(L) with norm(L), but 2*s(L) with only the diagonal elements).
        // 2*scale ist das größere Ideal
        vprintf Genus, 2: "NU: case 1\n"; 
        Append(~aL, FieldOfFractions(Order(p))!Uniformizer^(e+sL[i]));
      else
        // die Diagonalelemente bilden das größere Ideal:
        vprintf Genus, 2: "NU: case 2\n"; 
        Append(~aL, D[b] where _, b is Minimum([Valuation(x,p): x in D]));
      end if;
      Append(~uL, Valuation(aL[i], p));
      // Append(~wL, Minimum({e+sL[i]} join {uL[i] + QuadraticDefect(d/aL[i], p): d in D}));
    end for;

 // uL: the valuations of the upscaled norms
 // aL: the generators of the upscaled norms
 return uL, aL;
end function;

intrinsic IsSquarefreeLattice(L::LatMod, p::RngOrdIdl)->BoolElt
{ Check whether L is a squarefree lattice at the prime p. }
  // L any LatticeModule
   _,_,E := JordanDecomposition(L, p);
  return Set(E) subset {0, 1};

// TODO: is there a way to check "being squarefree" without doing expensive calculations (Jordan decomp., dual quotient ...)?
// scale integral testen, und beweisen dass es keine p^2-modulare Komponente gibt: mit p^-1 reskalieren und dualisieren. ist dieses Gitter integral?
end intrinsic;

intrinsic IsSquarefreeLatticeEverywhere(L::LatMod)->BoolElt
{ Check whether L is a squarefree lattice. }
  badP := [];
  for p in BadPrimes(L) do if not IsSquarefreeLattice(L, p) then Append(~badP, p); end if; end for;
  return #badP eq 0, badP;
end intrinsic;

intrinsic IsPrimitiveLattice(L::LatMod, p::RngOrdIdl)->BoolElt
{ Check whether L is primitiveat p, i.e., any Jordan decomposition consists of a 0-valuated component A and a 1-valuated component B, and rank(A) ge rank(B). }
  assert Minimum(p) eq 2;
  dims, scale := GenusSymbol(L, p);
  return (Set(scale) subset {0,1}) and (dims[1] le (Rank(L)/2));
end intrinsic;

intrinsic IsUnverzweigt(L::LatMod, p::RngOrdIdl) -> BoolElt
{ check whether L is unimodular at p and Norm(L) equals 2*Scale(L) at p. }
  // see Pfeuffer, p.11
  return IsUnimodular(L, p) and (Valuation(Norm(L), p) eq Valuation(2*Scale(L), p));
end intrinsic;

intrinsic IsUnverzweigtEverywhere(L::LatMod) -> BoolElt
{ check that IsUnverzweigt(L, p) holds for every p over 2 }
  if #[b: b in BadPrimes(L) | Minimum(b) ne 2] gt 0 then return false; end if;
  return &and[IsUnverzweigt(L, p):  p in [x[1]: x in Factorization(2*Integers(BaseRing(L)))]];
end intrinsic;

  
// 1. BONGs and good BONGs, as defined by C. Beli:

intrinsic RelativeQuadraticDefect(a::RngElt, p::RngOrdIdl) -> RngIntElt
{ Returns the valuation of the relative quadratic defect of a at the prime ideal p, as defined by C. Beli }
  q := QuadraticDefect(a, p);
  return q eq Infinity() select Infinity() else q - Valuation(a, p);
end intrinsic;

function IsInA(a, p) // needed for G_function
  e := Valuation(Order(p)!2, p);
  // cf. Lemma 3.5 in Beli: A is the set of all a \in F^* for which either (-a\not\in (F^*)^2 and D(-a) \subseteq O_p) or (-a\in (F^*)^2 and ord(a)\geq -2e).
  if (Type(QuadraticDefect(-a, p)) ne Infty) and (QuadraticDefect(-a, p) ge 0) then return true; end if;
  if (Type(QuadraticDefect(-a, p)) eq Infty) and (Valuation(a, p) ge -2*e) then return true; end if;
  return false;
end function;


function HasPropertyA(L, p) // needed for SpinorNorm
  error if not Minimum(p) eq 2, "property A only applies to lattices over dyadic fields";
  rL, sL, wL, aL, fL, G := GenusSymbol(L, p); 
  nL := [Valuation(aL[i], p): i in [1..#aL]];

  for i in [1..#G] do if rL[i] gt 2 then
    vprintf Genus : "[at a prime over 2]: property A is violated because there is a %o-dimensional Jordan component\n", Rank(G[i]);
    return false;
  end if; end for;
  
  // GenusSymbol: rL, sL, wL, aL, fL, G, note that aL only contains valuations
  for i in [1..#sL] do for j in [i+1..#sL] do
    if not ((0 lt nL[j] - nL[i]) and (nL[j] - nL[i] lt 2*(sL[j] - sL[i]))) then
    vprintf Genus : "[at a prime over 2]: property A is violated at %o, %o (norm/scale valuations do not fit)\n", i, j;
    return false; 
    end if;
  end for; end for;
  return true;
end function;

function N_function(a, g, p)
  // g is the mapping obtained by LocalMultiplicativeGroupModSquares(p).
  // This is expensive to calculate, so we pass it as an argument.

  b := a@@g;


  // project a into F^*/(O^*)^2:
  V := Parent(b);  // Full F_2-vector space of some dimension over GF(2)

  // treat the squares separately:
  if QuadraticDefect(a, p) eq Infinity() then return V; end if;

  B := Basis(V);
  
  // cf. paragraph before 1.2 in Beli 2003:
  // N(a) := N(F(\sqrt(a))/F), i.e. the values of the norm mapping of the quadratic extension
  // and N(a) is the orthogonal complement of (<a>(F^*)^2). 
  preimages := V!0;
  for i in [1..#B] do
    b := B[i];
    preim := b@g;
    // convert the Hasse symbol (as an element of the multiplicative group C_2)
    // to an element of the additive group F_2:
    // printf "Hilbert symbol with %o-th basis vector: %o\n", Index(B, b), HilbertSymbol(a, preim, p);
    preimages[i] := HilbertSymbol(a, preim, p) eq 1 select 0 else 1;
  end for;
  KernelGenerators := {B[i]: i in [1..#B] | preimages[i] eq 0};
  // (this fails if a is a square)
  j := Min({i: i in [1..#B] | preimages[i] eq 1});
  KernelGenerators join:= {B[i] + B[j]: i in [1..#B] | (preimages[i] eq 1) and (i ne j)};

  return sub<V | KernelGenerators>;
end function;

function G_function(a, V, g, p)
  // use LocalMultiplicativeGroupModSquares to calculate in F*/(O^*)^2
  // (or F^*/(F^*)^2 -- the last component of the result is the p-valuation
  // mod 2, and  F^*/(F^*)^2 =  F*/(O^*)^2 \times C_2.)
  // Also we expect V, g = LocalMultiplicativeGroupModSquares(p);

  // cf. Beli 2003, Def. 4.
  

  U := UnitSquareClassReps(p);
  assert U[1] eq 1;

  //b := a@@g;
  e := Valuation(Order(p)!2, p);
  R := Valuation(a, p);
  d := RelativeQuadraticDefect(-a, p);
  pi := UniformizingElement(p);
 
  if not IsInA(a, p) then     vprintf Genus : "G_function, case F\n"; return N_function(-a, g, p);
  elif ((-1/4)@@g eq a@@g) then     vprintf Genus : "G_function, case G\n"; return sub<V|[V.i: i in [1..Dimension(V)-1]]>;
  elif (R gt 4*e) then     vprintf Genus : "G_function, case H\n"; return sub<V| a@@g>;
  elif (2*e lt R) and (R le 4*e) then
    if (d le 2*e - R/2) then
    vprintf Genus : "G_function, case A\n";  
      return (sub<V| (a@@g)> + sub<V|(1 + pi^(R + d - 2*e))@@g>) meet N_function(-a, g, p);
    else
      vprintf Genus : "G_function, case B\n";
      assert GF(2)!R eq 0;
      return (sub<V| (a@@g)> + sub<V|(1 + pi^(R div 2))@@g>);
    end if;
  elif (R le 2*e) then 
    // printf "e = %o, R = %o, d = %o\n", e, R, d;
    if (d le e - R/2) then
    vprintf Genus : "G_function, case C\n";

    return N_function(-a, g, p);
    elif (e - R/2 lt d) and (d le 3*e/2 - R/4) then
      assert R mod 2 eq 0;
    vprintf Genus : "G_function, case D\n";

      return N_function(-a, g, p) meet sub<V| (1 + pi^(IntegerRing()!(R/2) + d - e))@@g>;
    else
    vprintf Genus : "G_function, case E\n";

    //assert R mod 4 eq 0;
      //assert e mod 2 eq 0;
      // attention! Use the Floor function wherever Beli writes stuff in square brackets. This is due to his citing Hsia's papers, which have this convention.
      return sub<V| (1 + pi^(e - (Floor(e/2 - R/4))))@@g>;
    end if;
  else error "this never happens.";
  end if;
end function;


intrinsic IsGoodBONG(BONG::[], p::RngOrdIdl) -> BoolElt
{Returns true iff BONG is a good BONG at p, as defined by C. Beli.}
  // Given: a BONG of a LatticeModule L at prime p.
  v := &and[Valuation(BONG[i], p) le Valuation(BONG[i+2], p): i in [1..#BONG-2]]; 
  return v;
end intrinsic;

procedure beli_correction(L, ~G, ~JJ, Steps, i, j, p);
  // this is helper procedure for MakeGoodBONG().
  // Correct the Jordan components with number i and j of L from the Jordan decomposition given in G and JJ,
  // in such a way that the new component with number i has maximal norm. 
  assert #Steps[i] eq 2;

  /*
  // for debugging:
  NU := NormUpscaled(G, p);
  // assert orthogonality of the vectors in JJ[Steps[i]] and those in JJ[Steps[j]], i.e. those that make up G[i] and G[j]:
  if Nrows(G[j]) eq 2 then assert IsZero( C*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[i][1], Steps[i][2]]) where C is  RowSubmatrix(JJ, [Steps[j][2]])); end if;

  assert NU[i] lt Valuation(Norm(LatticeModule(G[i])), p); // otherwise, the norm condition is violated
  assert NU[j] le Valuation(Norm(LatticeModule(G[j])), p); // "<=" must always hold, otherwise something about the lattice is broken
  */

  // we will need the original Gram matrix for re-orthogonalization:
  GI := G[i]^(-1); // over the quotient field
  assert #Steps[j] in [1,2];
  
  
  // if G[j] is two-dimensional, make sure that a norm generator is on the [1,1] position:
  if (Ncols(G[j]) eq 2) and (Valuation(G[j][1,1], p) gt Valuation(G[j][2,2], p)) then
    temp := JJ[Steps[j][1]];
    JJ[Steps[j][1]] := JJ[Steps[j][2]];
    JJ[Steps[j][2]] := temp;
    G[j] := B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[j][1]..Steps[j][#Steps[j]]]);
  end if;

  JJ[Steps[i][1]] +:= JJ[Steps[j][1]];
  
  // update Gram matrix for component #i:
  G[i] := B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[i][1], Steps[i][2]]);
  x := JJ[Steps[i][1]];
  y := JJ[Steps[i][2]];
  // Ansatz: v = JJ[Steps[j][1]] + lambda*x + mu*y
  
  // re-orthogonalize the first basis vector of G[j]:
  v := Vector([-G[j][1,1], 0]) * GI; // G[j][1,1] is always non-zero (the lattice is definite)

  // assert integrality at p:
  assert &and[Valuation(v[k], p) ge 0: k in [1..Ncols(v)]];
  JJ[Steps[j][1]] +:= v[1]*JJ[Steps[i][1]] + v[2]*JJ[Steps[i][2]];
  
  // if applicable, re-orthogonalize the second basis vector of G[j]:
  if Ncols(G[j]) eq 2 then
    w := Vector([-G[j][1,2], 0]) * GI; // G[j][1,2] is always non-zero (or else the lattice would have been further decomposed into two 1x1-lattices here)
    assert &and[Valuation(v[k], p) ge 0: k in [1..Ncols(w)]];
    JJ[Steps[j][2]] +:= w[1]*JJ[Steps[i][1]] + w[2]*JJ[Steps[i][2]];
  end if;

  // update Gram matrix for component #j:
  G[j] := B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[j][1]..Steps[j][#Steps[j]]]);

  /*
  assert Valuation(Norm(LatticeModule(G[i])), p) eq NU[i];
  */
end procedure;

function IsMaximalNormSplittingInternal(GramMatrices, Scales, Norms, Dimensions, p)
  // Scales: list of valuations of scales
  // Norms: list of generators of norms
  // Dimensions: list of dimensions of Jordan components
  // occurring in a Genus symbol of L at p, as calculated
  // by GenusSymbol(L, p). 

  Norms := [Valuation(Norms[i], p): i in [1..#Norms]];

  // test if each component is either unary or binary:
  fail := exists(i){i: i in [1..#Dimensions] | not(Dimensions[i] in {1,2})};
  if fail then 
    vprintf Genus : "IsMaximalNormSplittingInternal failed: components are not all unary or binary\n";
    return false, -i, {}; 
  end if;

  // test if binary components are modular:
  fail := exists(i){i: i in [1..#GramMatrices] | #E ne 1 where _,_,E is JordanDecomposition(LatticeModule(GramMatrices[i]), p)};
  if fail then
    vprintf Genus : "IsMaximalNormSplittingInternal failed: at least one of the binary components is not modular\n";
    return false, -i, {}; 
  end if;
  

  // test if sL[1] \supseteq sL[2] \supseteq ... \supseteq sL[#sL]:
  for i in [1..#Scales-1] do if (Scales[i] gt Scales[i+1]) then 
    error Sprintf("scale condition at components %o/%o not met, perhaps something is wrong with your lattice?\nFor your reference, the scale valuations are %o", i, i+1, Scales);
    return false, 0, {}; end if; 
  end for;

  // test if nL[i] = n(L^{sL[i]}):
  NU, _ := NormUpscaled(GramMatrices, p);
  failset := {}; // for testing purposes, we count the components that do not have maximal norm
  for i in [1..#Scales] do
    assert NU[i] le Norms[i];
    if NU[i] lt Norms[i] then 
      vprintf Genus : "IsMaximalNormSplittingInternal failed: norm condition at component %o\n", i;
      Include(~failset, i);
    end if;
  end for;
  if #failset gt 0 then return false, SetToSequence(failset)[1], failset; end if;
  return true, 0, {};  
end function;

intrinsic IsMaximalNormSplitting(G::[], p::RngOrdIdl) -> BoolElt
{returns true iff the given list G of Gram matrices corresponds to a maximal norm splitting at p.}
  sL, aL, _, _ := ScalesAndNorms(G, p, UniformizingElement(p));   
  Dimensions := [Nrows(a): a in G];
  return IsMaximalNormSplittingInternal(G, sL, aL, Dimensions, p);
end intrinsic;

intrinsic MaximalNormSplitting(L::LatMod, p::RngOrdIdl) -> [], []
{A maximal norm splitting of L at a dyadic prime p, as defined by C. Beli. Returns a Jordan decomposition into 1x1 and 2x2 components, and the corresponding list of basis vectors.}
  // we follow Beli, "Representation of Integral Quadratic Forms over Dyadic Local Fields", section 7

  require Order(p) cmpeq BaseRing(L) : "Incompatible arguments";
  require Minimum(p) eq 2 : "p must be a dyadic prime.";

  R := Order(p);
  e := RamificationIndex(p);
  Uniformizer := UniformizingElement(p);

  // we can use the additional values returned by DetailedJordanDecomposition to 
  // further decompose a Jordan decomposition into unary and binary components:

  //  dims, sL, wL, aL, fL, _ := GenusSymbol(L, p);
  JJJ, JG, JE, Blocks, Steps := DetailedJordanDecomposition(L, p);
  // join JJJ into one matrix of base vectors:
  JJ := Matrix(&cat[[J[i]: i in [1..Nrows(J)]]: J in JJJ]);


  // join the individual Gram matrices:
  A := DiagonalJoin(<g: g in JG>);

  // read the finer decomposition:
  G := [* *];
  for i in [1..#Steps] do
    // are we in a 1x1 or in a 2x2 block?
    if #Steps[i] eq 1 then
      AA := SubmatrixRange(A, j,j,j,j) where j is Steps[i][1];
      Append(~G, AA);
    elif #Steps[i] eq 2 then
      AA := SubmatrixRange(A, j,j,j+1,j+1) where j is Steps[i][1];
      assert Steps[i][2] eq Steps[i][1] + 1;
      Append(~G, AA);
    else 
      assert false; // something is wrong with DetailedJordanDecomposition
    end if;
  end for;
  Dimensions := [Nrows(a): a in G];
  assert SequenceToSet(Dimensions) subset {1,2};
  assert &+Dimensions eq Rank(L);
  
  // now turn this decomposition into unary/binary components into a maximal norm splitting:
  failset := {};
  while true do
    // Scale valuations sL, Norm generators aL, Norm valuations uL,
    // weight valuations wL of the decomposition:
    sL, aL, uL, wL := ScalesAndNorms(G, p, Uniformizer); 

    failset_old := failset;
    b, i, failset := IsMaximalNormSplittingInternal(G, sL, aL, Dimensions, p);
    assert (failset_old eq {}) or (#failset_old gt #failset);
    if b then break; end if; // maximal norm splitting reached!
    assert i gt 0; // if b is false and i=0, something is wrong with our lattice

    // here comes the interesting stuff:
    assert #Steps[i] eq 2; // unary components already have maximal norm.
  
    // the maximal norm splitting condition is violated at index i.
    find_j := [Valuation(aL[k], p): k in [i+1..#aL]] cat [2*(sL[i] - sL[k]) + Valuation(aL[k], p): k in [1..i-1]];
    assert #find_j gt 0;
    min, j := Minimum(find_j);

    if j le #aL - i then
      j := j + i;  // we want j to correspond to a Jordan component, not an index in find_j
      assert j gt i;

      // This is case 1 in Beli.
      beli_correction(L, ~G, ~JJ, Steps, i, j, p); 
    else 
      j := j - (#aL - i); // we want j to correspond to a Jordan component, not an index in find_j
      assert j lt i;
      
      // This is case 2 in Beli.
      // switch to a basis of the dual lattice L^#: 
      for k in [1..#Steps] do 
        for l in [1..#Steps[k]] do 
          JJ[Steps[k][l]] *:= Uniformizer^(-sL[k]); 
        end for; 
        assert Valuation(Scale(LatticeModule(B*L`Form*Transpose(B))), p) eq -sL[k] where B is RowSubmatrix(JJ, [Steps[k][1]..Steps[k][#Steps[k]]]);
      end for;

      // apply case 1 to the reversed orthogonal sum of the dual lattices:
      
      Steps := Reverse(Steps);
      // update Gram matrices for the reversed, dualized lattice:
      for k in [1..#G] do 
        G[k] :=  B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[k][1]..Steps[k][#Steps[k]]]);
        assert Nrows(G[k]) in {1,2};
      end for;
      
      assert #Steps[#aL - i + 1] eq 2; // component i is now at position #aL-i+1
      
      beli_correction(L, ~G, ~JJ, Steps, #aL - i + 1, #aL - j + 1, p);

      Steps := Reverse(Steps);
      G := Reverse(G); 

      // update norms/scales from the updated Gram matrices:
      sL, aL, uL, wL := ScalesAndNorms(G, p, Uniformizer); 
      // dualize again:
      for k in [1..#Steps] do for l in [1..#Steps[k]] do 
        JJ[Steps[k][l]] *:= Uniformizer^(-sL[k]); 
      end for; end for;

      // update Gram matrices:
      for k in [1..#G] do 
        G[k] :=  B*L`Form*Transpose(B) where B is RowSubmatrix(JJ, [Steps[k][1]..Steps[k][#Steps[k]]]);
        assert Nrows(G[k]) in {1,2};
      end for;
    end if;
  end while;

  assert &and[Nrows(G[k]) in {1,2}: k in [1..#G]];
  return G, JJ;
end intrinsic;


function MakeBONGdim2(L, p);
  // cf. Beli, Lemma 3.3
  // return a BONG of a 2-dimensional lattice.

  is_modular, r := IsModular(L, p);
  assert is_modular and (Rank(L) eq 2);
  
  pi := UniformizingElement(p);
  R := Order(p);
  
  A := L`Form;
  
  // is the (1,1) entry a norm generator?
  if Valuation(A[1,1], p) ne Valuation(Norm(L), p) then
    // swap stuff around so that the (1,1) entry will be a norm generator
    b := exists{i: i in [1..Rank(A)] | Valuation(A[i,i], p) eq Valuation(Norm(L), p)};
    if b then
      B := LocalBasis(L`Module, p);
      S := Matrix(B);
      SwapRows(~S, 1, 2);
    else
      B:= LocalBasis(L`Module, p);
      S:= Matrix(B);
      S[1] +:= S[2];
    end if;
    A := S*L`Form*Transpose(S);
    L_old := L;
    L := LatticeModule(A);
    assert Valuation(A[1,1], p) eq Valuation(Norm(L), p);
  end if;
  
  n := A[1,1];
  
  // Some final checks:
  disc := A[1,1]*A[2,2]- A[2,1]*A[1,2];
  val := Valuation(disc, p);
  assert GF(2)!val eq 0;
  assert Valuation(disc, p) eq 2*r;
  assert Valuation(n, p) ge 0;
  BONG := [n, disc / n];
  
  return BONG;
end function;


intrinsic MakeGoodBONG(L::LatMod, p::RngOrdIdl) -> []
{Return a good BONG of L at a dyadic prime p, as defined by C. Beli.}
  // Any BONG obtained from a maximal norm splitting is good, cf. Corollary 4.2 in Beli.
  // If a maximal norm splitting is calculated, a good BONG can be read off by joining
  // together the BONGs of the 1- and 2-dimensional components.
  require Order(p) cmpeq BaseRing(L) : "Incompatible arguments";
  require Minimum(p) eq 2 : "p must be a dyadic prime.";

  G, JJ := MaximalNormSplitting(L, p);
  BONG := [];
  pi := UniformizingElement(p);
  
  for i in [1..#G] do
    GG := G[i];
    if Nrows(GG) eq 2 then
      BONG cat:=MakeBONGdim2(LatticeModule(GG), p);
    elif Nrows(GG) eq 1 then
      Append(~BONG, GG[1,1]);
    else 
      assert false;
    end if;
  end for;  
  return BONG;  
end intrinsic;



// 2. Spinor norms of lattices over number fields:

intrinsic SpinorNorm(L::LatMod, p::RngOrdIdl) -> ModTupFld, BoolElt
{The spinor norm of L at p, as calculated by C. Beli for dyadic p and by M. Kneser for non-dyadic p. Returns a subspace of LocalMultiplicativeGroupModSquares(p), and a boolean which is true iff the spinor norm is exactly the units.}
  require Order(p) cmpeq BaseRing(L) : "Incompatible arguments";
  V, g := LocalMultiplicativeGroupModSquares(p);  
  e:= Valuation(BaseRing(L)! 2, p);

  if Minimum(p) ne 2 then
    // cf. Satz 3 in Kneser, "Klassenzahlen indefiniter quadratischer Formen in drei oder mehr Veränderlichen":
    J, G, E := JordanDecomposition(L, p);
    for i in [1..#G] do 
      if Rank(G[i]) ge 2 then
        vprintf Genus : "This lattice has a 2-dimensional Jordan constituent, and p is odd. Spinor norm is either F^* or O_F^*(F^*)^2, i.e. we will find a vector space of dimension %o or %o.\n", Rank(G[i]), Rank(V)-1, Rank(V);
        // well, which of the two is the case?
        if #{e mod 2: e in E} lt 2 then
          // the units only: 
          return sub<V | [V.i: i in [1..Dimension(V)-1]]>, true;
        else
          // the whole group:
          return V, false;
        end if;
      end if; 
    end for;
    // this is the obscure case where p is odd and each
    // of its Jordan components has dimension 1. In this case, the units
    // need not be contained in the Spinor norm.
   
    // generators of the (principal) norm ideals of the Jordan
    // components: since p is odd, the norms (and the scales) are just the 
    // dimensions of the Jordan components
    NormGens := [g[1,1]: g in G];

    TwoNormGens := [];
    for i in [1..#NormGens] do for j in [i..#NormGens] do 
      Append(~TwoNormGens, NormGens[i]*NormGens[j]);
    end for; end for;
    TwoNormVectors := [x@@g: x in TwoNormGens];
    
    vprintf Genus : "odd p, norm generators of the %o Jordan components are: \n\n%o\n\n Products of two elements are \n\n%o\n\nTwoNormVectors:\n%o\n\n", #G, NormGens, TwoNormGens, TwoNormVectors; 
    
    // cf. Kneser 1956, Satz 3:
    SN := sub<V| TwoNormVectors>;
  else
    SN := sub<V| 0>;
    BONG := MakeGoodBONG(L, p);
    // assert IsGoodBONG(BONG, p);
    if not HasPropertyA(L, p) then
      vprintf Genus : "This lattice does not have property A. Spinor norm is either F^* or O_F^*(F^*)^2, i.e. we will find a vector space of dimension %o or %o.\n", Rank(V)-1, Rank(V);
      // Using section 7, Thm. 3 in Beli 2002, one can decide which of the two cases applies. 
      // This is why we needed to compute a *good* BONG:
      for i in [1..#BONG-1] do
        BG := Basis(G_function(BONG[i+1]/BONG[i], V, g, p));
        for bg in BG do if bg[#Basis(V)] ne 0 then 
          // the whole group:
          return V, false, <Sprintf("G(%o=BONG[%o]/BONG[%o] contains non-units at entry #%o=%o", BONG[i+1]/BONG[i], i+1, i, Index(BG, bg), bg)>;
        end if; end for;
      end for;
      for i in [1..#BONG-2] do
        if Valuation(BONG[i], p) eq Valuation(BONG[i+2], p) then 
          assert GF(2)!(Valuation(BONG[i+1], p) - Valuation(BONG[i], p)) eq 0;
          if GF(2)!((Valuation(BONG[i+1], p) - Valuation(BONG[i], p))/2) ne GF(2)!e then
            // the whole group:
            return V, false, i;
          end if;
        end if;
      end for;
       
      // if all checks have passed, the Spinor norm is exactly the units:
      return sub<V | [V.i: i in [1..Dimension(V)-1]]>, true;
    end if;

    // cf. Beli 2003, Thm. 1
    for i in [2..Rank(L)] do SN +:= G_function(BONG[i]/BONG[i-1], V, g, p); end for;
  
    // for why we can take the Minimum in what follows here, see the remark on p. 161 in Beli 2003:
    if exists{i: i in [1..Rank(L)-2]| GF(2)!(Valuation(BONG[i+2], p) - Valuation(BONG[i], p)) eq 0} 
    then 
      alpha := IntegerRing()!Minimum({(Valuation(BONG[i+2], p) - Valuation(BONG[i], p))/2 : i in [1..Rank(L)-2]|  GF(2)!(Valuation(BONG[i+2], p) - Valuation(BONG[i], p)) eq 0});
      Uniformizer := UniformizingElement(p);
      SN +:= sub<V | (1+Uniformizer^alpha)@@g>;
    end if;
  end if;

  return SN, SN eq sub<V | [V.i: i in [1..Dimension(V)-1]]>;
end intrinsic;

intrinsic SpinorsNotExactlyTheUnits(L::LatMod) -> SeqEnum
{those places of BaseRing(L) where the Spinor norm is not exactly the units.}
  assert Rank(L) ge 3;
  // Rank >= 2 necessary for SpinorNorm().
  // Rank >= 3 necessary for \Phi from O'Meara 102:7 to be surjective.
  PP := Factorization(Discriminant(L));
  S := []; // set of bad primes
  for i in [1..#PP] do
    P := PP[i][1];
    SN, b := SpinorNorm(L, P);
    /*
    // we double-check:
    ExactlyTheUnits := sub<V| [V.i: i in [1..Rank(V) - 1]]> where V is LocalMultiplicativeGroupModSquares(P);  

    assert b eq (SN eq ExactlyTheUnits);
    */
    if not b then Include(~S, [*P, SN*]); end if;
  end for;
  return S;
end intrinsic;



// 3. Enumerating genera of lattices over number fields:

intrinsic MapIdeleIntoClassGroup(L::LatMod, Idele::[]: AtInfinity:=false) -> GrpAbElt
{ map an idele into the class group associated to L. The parameter Idele must be a list of tuples <p, a_p>, where p is a prime of BaseRing(L), and a_p is an element of K^* (interpreted as an element of the completion K^*_p). The parameter AtInfinity can be a list of tuples <v, +1 or -1>, where v is an element of InfinitePlaces(BaseRing(L)). All places, finite or infinite, which are unspecified are interpreted as 1.}

  if not assigned(L`GenusComp) then PrepareClassGroup(L); end if;
  assert L`GenusComp`Status gt 0; // class group must be initialized by now
  
  
  // sanity checks: 
  error if not &and[(Type(i) cmpeq Tup) and (#i eq 2) and ISA(Type(i[1]), RngOrdIdl) and (Order(i[1]) cmpeq BaseRing(L)) and ISA(Type(i[2]), FldNumElt) /*and (Valuation(i[2], i[1]) ge 0)*/: i in Idele], "Idele: must be a list of tuples <p, a_p> with p::RngOrdIdl, a_p::FldNumElt";
  F := NumberField(BaseRing(L));
  error if not IsTotallyReal(F), "The base field is not totally real. Please come back when you have found a lattice over a totally real number field.";
  R := Integers(F);
  IP := InfinitePlaces(F);  
  the_idele_inf := [1: I in IP];
  if AtInfinity cmpne false then 
    error if not Type(AtInfinity) eq SeqEnum, "AtInfinity: must be false, or a list";
    error if not &and[(Type(i) cmpeq Tup) and (#i eq 2) and (i[2] in [-1, 1]) and (i[1] in IP): i in AtInfinity], "AtInfinity: must be a list of tuples <v, +1/-1> with v in InfinitePlaces(BaseRing(L))";
    for i in AtInfinity do
      the_idele_inf[Index(IP, i[1])] := i[2];
    end for;
  end if;
  
  // The finite places we need to make congruent to 1 in order to be able to map into the class group:
  the_idele := [R!1: p in L`GenusComp`RayPrimes]; 
  for i in Idele do 
    j := Index(L`GenusComp`RayPrimes, i[1]);
    if j gt 0 then the_idele[j] := i[2]; end if;
  end for;
  t := the_idele; // this is just used to check an assertion later
  
  // L`GenusComp`hom_RCG is the map from L`GenusComp`RCG into the divisors that Magma returns. 
  // So, only the ideals D that are coprime to the ray  lie in the codomain of the map and can be mapped 
  // into L`GenusComp`RCG via I@@L`GenusComp`hom_RCG.
  // The ideles we are considering here are considered to be representatives of classes of ideles (modulo (F^*_S)*J^L, where F^*_S is the set of principal ideles which are totally positive at the infinite places where the envelopping space V of L is anisotropic, and where J^L is the set of ideles which lie in the spinor norm of L at all finite places. 
  // So, we may modify the ideles in the following two ways without changing their class:
  // 1) multiply every component by an element of F^*
  // 2) multiply any individual component (at p) by an element of the spinor norms of L_p. In particular, any component is only relevant modulo squares.
  // Thus we can achieve that the idele in question is congruent to 1 modulo the ray M that defines L`GenusComp`RCG.
  // The idele is then interpreted as the ideal \prod_p{p^Valuation(idele_p, p)}, which can be mapped into the class group by Magma.


  
  // first, we need to be able to invert modulo the divisor L`GenusComp`M:
  if &or[Valuation(the_idele[i], L`GenusComp`RayPrimes[i]) gt 0: i in [1..#the_idele]] then
    // we need to correct with an element of F^* which is positive at the relevant infinite places.
    // Since InverseMod only works for finite places, we will first correct only the finite places ...
    s := WeakApproximation(L`GenusComp`RayPrimes, [-Valuation(the_idele[i], L`GenusComp`RayPrimes[i]): i in [1..#L`GenusComp`RayPrimes]]); 

    // Note that s can of course be non-trivial at some places which are not in RayPrimes, including infinite places.
    
    // first we take care of the places in RayPrimes:
    // multiply with s at every place, and find nice representatives modulo squares in each component:
    for j in [1..#the_idele] do
      the_idele[j] := R!quo<R| L`GenusComp`MM[j]>!(the_idele[j]*s);
    end for;        

    // check that everything is fine now:
    assert &and[Valuation(the_idele[i], L`GenusComp`RayPrimes[i]) eq 0: i in [1..#L`GenusComp`RayPrimes]];
  else s := F!1;
  end if;

  // the idele we have now:
  // equals the_idele[j] at RayPrimes[j]
  // equals s*Idele[j][2] at a finite place Idele[j][1] which does not divide L`GenusComp`M
  // equals s at all other finite places
  // equals (the IP[j]-embedding of s)*the_idele_inf[j] at IP[j]
  
  // we can now invert our idele modulo L`GenusComp`M:
  // (note that L`GenusComp`M = &*L`GenusComp`MM):
  if #the_idele gt 0 then el := InverseMod(CRT(the_idele, L`GenusComp`MM), L`GenusComp`M); else el := R!1; end if;

  // the CRT with L`GenusComp`MM produces an element that is a square at RayPrimes, as per the definition of L`GenusComp`MM and the local square theorem.
  // The idele is now a square at all the p in RayPrimes. We write a "1" instead at all these places (we may do this since the ideles are really classes of ideles, modulo K^* and modulo ideles which lie in the spinor norms at all finite places. The squares at p always lie in the spinor norms at p.
  
  // the idele we have now:
  // equals 1 at RayPrimes[j] for all j
  // equals el*s*Idele[j][2] at any finite place Idele[j][1] which does not divide L`GenusComp`M  
  // equals el*s at all other finite places
  // equals (the IP[j]-embedding of el*s)*the_idele_inf[j] at IP[j]
    
  // We still need to take care of the infinite places. 
  
  the_idele_inf := [the_idele_inf[j] * Evaluate(el*s,IP[j]) lt 0 select -1 else 1: j in [1..#IP]]; 
    
  // This version of CRT takes a global integer for the finite places, and can correct signs at the infinite places:
  el2 := CRT(L`GenusComp`M, [1..#IP], R!1, the_idele_inf);

  // el2 changes at the p in RayPrimes only by a local square. Since we are only interested in square classes,
  // we view the idele as w.l.o.g. unchanged in the p in RayPrimes.
    
  // the idele we have now:
  // equals 1 at RayPrimes[j] for all j
  // equals el2*el*s*Idele[j][2] at a finite place Idele[j][1] which does not divide L`GenusComp`M  
  // equals el2*el*s at all other finite places
  // equals 1 at IP[j] for all j (it doesn't hurt to make it positive at all infinite places, even if only a subset of these defines our ray class group)
  
  assert &and[IsOne(quo< R | L`GenusComp`MM[k] > ! (el2*el*s * t[k])): k in [1..#t]];    // t is just a copy of the value the_idele had in the beginning
  assert &and[ Evaluate(el2*the_idele_inf[j], IP[j]) gt 0: j in [1..#IP]];

  // This idele can be sent into the ray class group L`GenusComp`RCG via the (inverse of the) homomorphism returned by Magma.
  // We first interpret it as the ideal which will actually have to be mapped:
  // i.e., we just collect the p-valuations at the noncritical places (p notin RayPrimes):
  temp := [Idele[j][1]^Valuation(Idele[j][2], Idele[j][1]): j in [1..#Idele]]; // why are we including all the primes here, even those in L`RayPrimes?: -- The correction factor s is still included below, it may have negative valuations at the p in RayPrimes. So we need the original values as well.
  if #temp gt 0 then ideal := &*temp; else ideal := R; end if;  
  
  // then multiply by the correcting scalars, and map the resulting ideal into the ray class group:
  // the second return value is the element in the factor group of the ray class group actually corresponding to a proper spinor genus in Genus(L).
  if L`GenusComp`Status gt 1 then
    return g, g@L`GenusComp`hom_RCG_FactorGroup where g is (ideal*el*el2*s)@@L`GenusComp`hom_RCG;
  else
    return g where g is (ideal*el*el2*s)@@L`GenusComp`hom_RCG;
  end if;
end intrinsic;


intrinsic MapPrincipalIdeleIntoClassGroup(L::LatMod, a::FldNumElt) -> GrpAbElt
{Map the principal idele defined by the element 'a'  into the ray class group identified with the proper spinor genera in Genus(L) } 
  require Parent(a) cmpeq NumberField(BaseRing(L)) : "Incompatible arguments";
  
  F := NumberField(BaseRing(L));
  decomp := Factorization(a*Integers(F));
  primes := {x[1]: x in decomp};
  
  // Finite places:
  the_idele := [];
  for p in primes do 
    Append(~the_idele, <p, F!a>); // a@@g
  end for;
  // Infinite places:
  IP := InfinitePlaces(F);
  the_idele_infinite := [];
  for j in [1..#IP] do
    Append(~the_idele_infinite, <IP[j], Evaluate(a, IP[j]) gt 0 select 1 else -1>);
  end for;
  
  return MapIdeleIntoClassGroup(L, the_idele: AtInfinity:=the_idele_infinite);
end intrinsic;

  
function GetCriticalPrimes(L: FullList:=false)
  // Set FullList:=true to calculate one prime for every element of L`GenusComp`RCG_FactorGroup,
  // instead of just one prime for every generator. This is needed for indefinite lattices.
  
  if (not assigned L`GenusComp) or (L`GenusComp`Status lt 2) then PrepareClassGroup(L); end if;
  
  CriticalPrimes := [];
  F := NumberField(BaseRing(L));
  MyRCG_FactorGroup := FullList select {Identity(L`GenusComp`RCG_FactorGroup) } else sub<L`GenusComp`RCG_FactorGroup | [] >; // set to the trivial subgroup and not to [], otherwise a trivial element might be picked as generator below
  
  maxnorm := 50; // this is really arbitrary -- we just grab some primes to start with and will calculate more later if needed.
  GoodPrimes := [p: p in PrimesUpTo(maxnorm, F: coprime_to:=&*BadPrimes(L))];
  p_ind := 1;
  RCG_FactorGroupGens := []; 
  
  
  if FullList then vprintf Genus : "Looking for %o group elements.\n", #L`GenusComp`RCG; end if;

  while #MyRCG_FactorGroup lt #L`GenusComp`RCG_FactorGroup do
    while p_ind gt #GoodPrimes do
      // calculate more good primes!
      maxnorm := Floor(maxnorm * 1.2); 
      GoodPrimes := [p: p in PrimesUpTo(maxnorm, F: coprime_to:=&*BadPrimes(L))];
    end while;
    p := GoodPrimes[p_ind]; 
    
    
    // map the (1,...,1,pi,1,...,1) idele into the class group, where pi is a uniformizing element of p:
    _, h := MapIdeleIntoClassGroup(L, [<p, F!UniformizingElement(p)>]); // @L`GenusComp`hom_RCG_FactorGroup;
    
    if not(h in Keys(L`GenusComp`GroupElementsToPrimes)) then L`GenusComp`GroupElementsToPrimes[h] := p; end if;
    
    oldsize := #MyRCG_FactorGroup;
    Append(~RCG_FactorGroupGens, h);    
    if not FullList then 
      MyRCG_FactorGroup := sub<L`GenusComp`RCG_FactorGroup | RCG_FactorGroupGens>;
    else
      MyRCG_FactorGroup join:= {h};
    end if;
    if #MyRCG_FactorGroup gt oldsize then 
      vprintf Genus : " -- subgroup size increased from %o to %o.\n", oldsize, #MyRCG_FactorGroup;
      // this prime ideal is relevant for neighbour generation, store it:
      Append(~CriticalPrimes, p);
    end if;  
    p_ind +:= 1;
  end while;
  return CriticalPrimes;
end function;
  
  
  
intrinsic PrepareClassGroup(L::LatMod)
{internal use.}
  deletedata(L);
  
  F := NumberField(BaseRing(L));
  assert IsIdentical(BaseRing(L), Integers(F)); // debug 20150617: this should never fail!
  BP := SetToSequence(BadPrimes(L)); // the primes dividing 2*volume(L)  
  
  // we only need to carry around those finite places where the Spinor norm is not exactly the units:
  RayPrimes := []; 

  Gens := [];
  for p in BP do
    V, g := LocalMultiplicativeGroupModSquares(p);  
    Spinors, ExactlyTheUnits := SpinorNorm(L, p);
    
    if not ExactlyTheUnits then 
      vprintf Genus : "Found a prime over %o where the Spinor norm is not exactly the units of the order. dim(Spinors)=%o, dim(LocalMultGrpModSq)=%o\n", Minimum(p), Dimension(Spinors), Dimension(V);
      Append(~RayPrimes, p);
      // a basis of the Spinor norm of L at p, when viewed (modulo squares) as an F_2-vector space  (cf. LocalMultiplicativeGroupModSquares)
      b := Basis(Spinors); 
      Gens cat:= [<F!(a@g), p>: a in b];
      assert &and[Valuation(gg, p) in {0,1}: gg in [a@g: a in b]];  // we rely on LocalMultiplicativeGroupModSquares mapping to 0- and 1-valued elements, and on the last base vector of V to represent a uniformizer of p
    // else
      // printf "At prime over %o, the Spinor norm is exactly the units of the order.\n", Minimum(p);    
    end if;
  end for;

  // calculate a ray M, and for convenience store its factors as well:
  L`GenusComp`M := 1*Integers(F);
  // the individual powers of prime ideals composing the ray:
  L`GenusComp`MM := [];  
  // the places used for the ray:
  L`GenusComp`RayPrimes := RayPrimes; 

  for i in [1..#RayPrimes] do
    p := RayPrimes[i];
    e_p := Valuation(Order(p)!2, p);
    L`GenusComp`M *:= p^(1 + 2*e_p);
    Append(~L`GenusComp`MM, (p)^(1 + 2*e_p));
  end for;
  // the ray we define includes all infinite places:
  L`GenusComp`RCG, L`GenusComp`hom_RCG := RayClassGroup(L`GenusComp`M, [1..Signature(F)]);
  
  // Step 1: map the generators into the class group to create the factor group.
  RCG_SubgroupGens := [];
  RCG_Subgroup := sub<L`GenusComp`RCG | 2*L`GenusComp`RCG>;
  
  L`GenusComp`Status := 1; // meaning: "initialized to the point that MapIdeleIntoClassGroup can be called"

  for g in Gens do
    // h := MapIntoClassGroupInternal(L, g[1], g[2]); 
    h := MapIdeleIntoClassGroup(L, [<g[2], g[1]>]);
    
    Append(~RCG_SubgroupGens, h); 
    RCG_Subgroup := sub<L`GenusComp`RCG | RCG_SubgroupGens, 2*L`GenusComp`RCG>;
  end for;

  L`GenusComp`RCG_Subgroup := RCG_Subgroup;
  L`GenusComp`RCG_FactorGroup, L`GenusComp`hom_RCG_FactorGroup := quo<L`GenusComp`RCG | RCG_Subgroup>;
  
 vprintf Genus : "*** PrepareClassGroup: ray class group: size = %o\n", #L`GenusComp`RCG;
 vprintf Genus : "*** PrepareClassGroup: subgroup of ray class group: size = %o\n", #RCG_Subgroup;


  // Step 2: find good generators (that is: generators which are not dyadic, and not in BadPrimes(L) -- so that neighbour generation is always possible), 
  // of this factor group.
  // Only pick ideles of the form (1,...,1,p,1,...,1) so that:
  // a)  nothing has to be corrected before they can be mapped down into the factor group
  // b)  we can simply store the prime ideal and know that it corresponds to the (1,...,1,p,1,...,1) 
  // These prime ideals will be stored in L`GenusComp`CriticalPrimes.
  
  // the primes which generate the spinor genera in Genus(L):
  // if L is not definite, we need one prime for each factor group element.
  // if L is definite, we only need a generating set of primes for the factor group.
  
  L`GenusComp`Status := 2; // meaning: "you can calculate a set of critical primes now"
  
  L`GenusComp`CriticalPrimes := GetCriticalPrimes(L: FullList := true); // FullList := not IsDefinite(L));
  
  L`GenusComp`Status := 3; // meaning: "ready, go ahead and do something useful with it"

    
// vprintf Genus : "*** good primes over %o (together with the squares) generate the subgroup.\n", [Minimum(q): q in L`GenusComp`CriticalPrimes];
// vprintf Genus : "*** PrepareClassGroup: factor group of ray class group: size = %o (this is the number of Spinor `+` genera in Genus(L))\n", #L`GenusComp`RCG_FactorGroup;
end intrinsic;


function SmallestNormGoodPrime(L)
  // return a smallest-norm prime which is not contained in BadPrimes(L),
  // i.e. not dyadic and does not divide Discriminant(L) -- making neighbor generation possible (and hopefully fast!).
  // Note that this prime can be contained in L`GenusComp`CriticalPrimes, and that this is not a problem!
  
  if (not assigned(L`GenusComp)) or (L`GenusComp`Status lt 3) then PrepareClassGroup(L); end if; // must be fully initialized by now
  
  if #L`GenusComp`CriticalPrimes gt 0 then maxnorm := Maximum({Norm(a): a in L`GenusComp`CriticalPrimes}); else maxnorm := 10; end if;
  repeat
    PP := [q: q in PrimesUpTo(maxnorm, NumberField(BaseRing(L)))| not(q in BadPrimes(L))];
    maxnorm := NextPrime(maxnorm);
  until #PP gt 0;
  
  return PP[1];
end function;

// for use in for Genus enumeration:
intrinsic GetPrimes(L::LatMod) -> []
{A complete set of primes of BaseRing(L) needed to enumerate the genus of L.}
  require Rank(L) ge 3 : "The rank of L must be at least 3."; // otherwise we are unsure about the Spinor norms at unimodular primes.
  if (not assigned(L`GenusComp)) or (L`GenusComp`Status lt 3) then PrepareClassGroup(L); end if; // must be fully initialized by now
  
  return L`GenusComp`CriticalPrimes;
end intrinsic;

intrinsic IteratedNeighbours(L::LatMod, p::RngOrdIdl: UseAuto:= true) -> []
{ iterated neighbours at the prime p. }
  require IsDefinite(L) : "L must be definite"; // otherwise we cannot test for isometry
  
  todolist := [L];
  Result := [];
  while #todolist gt 0 do
    K := todolist[#todolist];
    Remove(~todolist, #todolist);
    Append(~Result, K);

    // keep if not isometric, continue until the whole graph has been exhausted
    call := func < M | not (&or[IsIsometric(M, x): x in (Result cat todolist)]), true >;

    N := Neighbours(K, p: CallBack:=call, AutoOrbits:=UseAuto); 
    // note that none of the lattices returned in N are isometric to a lattice in todolist or in Result,
    // but the individual lattices in N may be isometric to each other.
    N2 := []; for n in N do if not &or[IsIsometric(m, n): m in N2] then Append(~N2, n); end if; end for; 
    todolist cat:= N2;
  end while;  
  return Result;
end intrinsic;


function SpinorGeneraModuloSomething(L, g: CalculateCriticalPrimes:=true, GeneratorsOnly:=true);
  error if (not assigned(L`GenusComp)) or (L`GenusComp`Status lt 3), "The genus computation record of L has not been fully initialized.";
  error if not (Parent(g) eq L`GenusComp`RCG), "incompatible arguments, g is not an element of the class group";
  
  // g is an element of L`GenusComp`RCG. Calculate the subgroup H generated by < L`GenusComp`RCG_Subgroup , g > .
  // Then the quotient Q := L`GenusComp`RCG / H is a quotient of L`GenusComp`RCG_FactorGroup. We return this quotient, and optionally return a set of critical primes for this quotient Q as well.
  F := NumberField(BaseRing(L));

  NewSubgroup := sub<L`GenusComp`RCG | L`GenusComp`RCG_Subgroup, g>;
  NewRCG_FactorGroup, hom_NewRCG_FactorGroup := quo<L`GenusComp`RCG | NewSubgroup>;
  
  if not CalculateCriticalPrimes then return NewRCG_FactorGroup, hom_NewRCG_FactorGroup, _; end if;
  
  // TODO 20160111: oh god, redo this in a sane way ...
  // list of group elements which represent the factor group:
  to_find  := [(x@@hom_NewRCG_FactorGroup)@L`GenusComp`hom_RCG_FactorGroup : x in (GeneratorsOnly select Generators(NewRCG_FactorGroup) else NewRCG_FactorGroup)]; 
  Crit := [];

  // we use the mapping { Class group elements } -> { prime ideals } that we calculated in PrepareClassGroup:
  Crit := [L`GenusComp`GroupElementsToPrimes[g]: g in to_find];

  return NewRCG_FactorGroup, hom_NewRCG_FactorGroup, Crit;  
end function;

procedure AddLatticeAndItsIteratedNeighbours(L, p, ~Reps, ~Markers: UseAuto, Max)
  // for use in GenusRepresentatives
  if not &or[IsIsometric(r, L): r in Reps] then
    Append(~Reps, L);
    Append(~Markers, false);
  end if;
  N:= IteratedNeighbours(L, p: UseAuto:=UseAuto);
  for n in N do
    if forall{ r: r in Reps | &and[not IsIsometric(r, n): r in Reps]} then 
      Append(~Reps, n);
      Append(~Markers, true);
      if #Reps ge Max then break; end if;
    end if;
  end for;
end procedure;



intrinsic GenusRepresentatives(L::LatMod : Max:= Infinity(), UseAuto:= true, Extra:= 0, NoBreakOnError:=false, IndefiniteClassNumberOnly:=false) -> []
{A list of representatives of the isometry classes in the genus of L. Optional parameters: Max, the maximum number of isometry classes to calculate; Extra, the number of extra primes for safety}

  require Rank(L) ge 3 : "Only implemented for lattices of rank at least 3."; // otherwise the isomorphism to the class group fails, cf. §102 in O'Meara.
  requirege Extra, 0;
  require (Max ge 1): "Max should be an integer or Infinity().";
  
  
  // There are infinitely many primes generating our ray class group for L,
  // so there is no need to choose dyadic primes. Neighbour generation
  // for dyadic primes is a pain and should be avoided.
  
  F := NumberField(BaseRing(L));
  //// Primes join:= {good};
  
  
  ///Reps:= [ L ]; // the todo-list, and result

  if IsDefinite(L) then
      // our favorite prime of small norm:
      good := SmallestNormGoodPrime(L);
      ///  if good in Primes then Primes := {x: x in Primes | x ne good}; end if; --> if not commented out, this line causes lattices to be found in the extra rounds. Since it doesn't do anything bad, we can just leave it in here ...
      
      // Let's see if the favorite prime switches spinor genera too:
      g, h := MapIdeleIntoClassGroup(L, [<good, F!UniformizingElement(good)>]);
      if not IsIdentity(h) then
        _,_,Primes := SpinorGeneraModuloSomething(L, g);
        // at this point, 'good' will not be an element of Primes.
        vprintf Genus : "We have %o critical primes (modulo our good prime).\n", #Primes;
      else
        Primes := {x: x in GetPrimes(L)};
        // here we are not sure whether 'good' is also an element of primes, and we can remove it to save some time:
        if good in Primes then
          /// error "why is this -- we thought h was the identity of our class group quotient?";
          // TODO 20160111: das liegt daran dass wir die CriticalPrimes mit FullList:=true berechnen und deswegen auch einen Repräsentanten für die 1 berechnen ._.
          Primes diff:= {good};
          vprintf Genus : "We have %o critical primes, after removing the good prime from this set.\n", #Primes; 
        else
          vprintf Genus : "We have %o critical primes.\n", #Primes;
        end if;
      end if;
      error if not &and[Minimum(q) gt 2: q in Primes], "PrepareClassGroup should not pick dyadic primes.";
      error if &or[p in BadPrimes(L):  p in Primes], "PrepareClassGroup should not pick any primes where the lattice is not modular.";
      
      
///      Markers := [ false ]; 
      extra:= false; // are we in the extra rounds (for safety?)
      ExtraFail := false; // will be true if an error has happened and NoBreakOnError set to true

      // We begin with the (proper) spinor genus of the given latice:
      Reps := [  ]; // representatives of the classes in Genus(L)
      Markers := [  ]; // marker for "do we have to visit this node?"
      vprintf Genus : "now going to do all neighbours at good prime. Old number of lattices was %o.\n", #Reps;
      AddLatticeAndItsIteratedNeighbours(L, good, ~Reps, ~Markers : UseAuto:=UseAuto, Max:=Max);
      
      while true do
        i:= 1;
        while i le #Reps and #Reps lt Max do
          // If Markers[i] is set, then a different representative of Reps[i]'s spinor genus has been encountered before.
          // We neither need to calculate its neighbourhood graph at our favorite prime, nor do we need to apply the critical primes.
          if Markers[i] then i+:=1; continue; end if;

          Markers[i] := true;
        
          for p in Primes do
            // one neighbour each at all of the critical primes:
            // set AutoOrbits to false unless we wish wo calculate all of the neighbors!
            /// call_one:= func< LL | not iso, iso where iso:= IsIsometric(Reps[i], LL) >;
            vprintf Genus : "One neighbour at prime:\n%o\n\n", p;
            N:= Neighbours(Reps[i], p: Max:=1, /*CallBack:=call_one, */ AutoOrbits:= false);
            assert #N eq 1;

            for n in N do
              if forall{ r: r in Reps | not IsIsometric(r, n) } then 
                // add the remainder of the (proper) spinor genus of n:
                vprintf Genus : "now going to do all neighbours at good prime. Old number of lattices was %o.\n", #Reps;
                AddLatticeAndItsIteratedNeighbours(n, good, ~Reps, ~Markers : UseAuto:=UseAuto, Max:=Max);
                
                // are we in the extra rounds? Then we should NOT find anything new!
                if extra then
                  error if not NoBreakOnError, "Oops, found new rep in the extra rounds! --- Please report";
                  ExtraFail := true;
                end if;
                if #Reps ge Max then break p; end if;
              end if;
            end for;
          end for;
          i +:= 1;
        end while;

        // when we arrive here, then nothing new can be found with the current set of primes.
        // If we were told to use extra primes, then now is the time to start using them:
        if (#Reps ge Max) or extra or (Extra cmpeq 0) then 
          assert (#Reps ge Max) or &and[b: b in Markers];
          if (Extra gt 0) then
            return Reps, ExtraFail; // ExtraFail = true if there were new reps in the extra rounds
          else 
            return Reps;
          end if;
        end if;

        // Add some more extra rounds for safety and repeat
        Primes:= GetExtra(L, Primes, Extra);
        //if Extra gt 0 then printf "Double-checking with %o extra primes over %o...\n", Extra, [Minimum(a): a in GetExtra(L, Primes, Extra)]; 
        // end if;
        extra:= true;
        Markers := [false: i in [1..#Reps]];
      end while;


  
  else
      // The indefinite case.
      // What makes this case easier: proper spinor genera and proper classes coincide (when rank(L) >= 3, cf. 104:5 in O'Meara). 
      // What makes this case harder: Aut(L) is infinite. We cannot test for isometry of global lattices.
      // How do we still achieve a system of representatives for the classes in Genus(L), and avoid hitting a class twice if it decomposes into two proper classes?: We explicitly calculate the group element corresponding to the improper automorphism in O(V) that switches two proper classes, and kick it out. The group element can be calculated by explicitly mapping the spinor norm of the reflection at any basis vector of L into the class group.
      
      G := GramMatrix(L);
      A, T := OrthogonalizeGram(G); // now T*G*Transpose(T) = A
      //M :=  DiagonalMatrix(BaseRing(T), [1: i in [1..Nrows(T)]]);
      // M[1,1] := -1; 

      ind := [i: i in [1..Nrows(G)] | not IsZero(G[i][i])];
      if #ind eq 0 then
        ind := [i: i in [2..Ncols(G)] | not IsZero(G[1][i])];
        error if #ind eq 0, "Lattice does not have full rank.";
        // now Phi(v_1+v_i, v_1+v_i) = 2*Phi(v_1, v_i), for any i in ind.
        spinornorm := 2*G[i][ind[1]];
      else
        spinornorm := G[1,1];
      end if;
      // now spinornorm is a representative of the spinor norm of the reflection at x_i for some i, or at (x_1+x_i) for some i.
            
      
      h, g := MapPrincipalIdeleIntoClassGroup(L, F!spinornorm);
      // g := h@L`GenusComp`hom_RCG_FactorGroup;
      
      /// assert IsIdentity(g@L`GenusComp`hom_RCG_FactorGroup) eq (g in L`GenusComp`RCG_Subgroup);
      
     
      if not IsIdentity(g) then // equally: if not IsIdentity(g@L`GenusComp`hom_RCG_FactorGroup)
        // We factor out a group element that switches proper spinor genera (= proper classes in the indefinite case);
        // but only switches to a different proper spinor genus inside the same isometry class.
        printf "Classes and proper classes do not coincide in Genus(L), i.e. there is no automorphism of L with determinant -1.\n";
        printf "Adding an element to RCG_Subgroup ...";
        MyRCG_FactorGroup, Myhom_RCG_FactorGroup, Crit := SpinorGeneraModuloSomething(L, h: CalculateCriticalPrimes := not IndefiniteClassNumberOnly, GeneratorsOnly:=false);
        
        if IndefiniteClassNumberOnly then return #MyRCG_FactorGroup; end if;
      else
        printf "Classes and proper classes coincide in Genus(L), i.e. there is an automorphism of L with determinant -1.\n";
        MySubgroup := L`GenusComp`RCG;
        MyRCG_FactorGroup := L`GenusComp`RCG_FactorGroup; Myhom_RCG_FactorGroup := L`GenusComp`hom_RCG_FactorGroup;
        if IndefiniteClassNumberOnly then return #MyRCG_FactorGroup; end if;
      end if;  

      // remove the representative of the trivial element:
      Crit := [p: p in L`GenusComp`CriticalPrimes | not IsIdentity(MapIdeleIntoClassGroup(L, [<p, F!UniformizingElement(p)>])@L`GenusComp`hom_RCG_FactorGroup)];
      
      vprintf Genus : "Number of elements: RCG: %o, Subgroup: %o, RCG_FactorGroup: %o, MySubgroup: %o, MyRCG_FactorGroup: %o, Crit: %o\n", #L`GenusComp`RCG, #L`GenusComp`RCG_Subgroup, #L`GenusComp`RCG_FactorGroup, #MySubgroup, #MyRCG_FactorGroup, #Crit;
      
      ClassGroupElements := [x: x in MyRCG_FactorGroup];
      
      for p in Crit do
        // one neighbour each at all of the critical primes:
        // We now only need to calculate one neighbour for each critical
        // prime and have a complete set of representatives of the improper(!) isometry classes in Genus(L).
        /// call_one:= func< LL | true, false >;
        N:= Neighbours(L, p: /*CallBack:=call_one, */ Max:=1, AutoOrbits:= false);
        assert #N gt 0;
        Append(~Reps, N[1]);
      end for;
      return Reps;
  end if;
end intrinsic;

intrinsic IsOneClassGenus(L::LatMod: UseAuto:=false) -> BoolElt
{returns true iff #GenusRepresentatives(L) eq 1, but is much faster than calling GenusRepresentatives.}
  require Rank(L) ge 3 : "Only implemented for lattices of rank at least 3.";
  
  if not IsDefinite(L) then return GenusRepresentatives(L: IndefiniteClassNumberOnly) eq 1; end if;
  
  F := NumberField(BaseRing(L));
 
  assert IsIdentical(BaseRing(L), Integers(F)); // debug 20150617: this should never fail! There were some Magma errors that made this fail, I still don't know why ...
  
  if (not assigned(L`GenusComp)) or (L`GenusComp`Status lt 3) then PrepareClassGroup(L); end if; // must be fully initialized by now
  
  // There are infinitely many primes generating our ray class group for L,
  // so there is no need to choose dyadic primes. Neighbour generation
  // for dyadic primes is a pain and should be avoided.
  error if not &and[Minimum(q) gt 2: q in L`GenusComp`CriticalPrimes], "PrepareClassGroup should not pick dyadic primes.";
  
  good:= SmallestNormGoodPrime(L);

  
  vprintf Genus : "IsOneClassGenus: working with %o critical primes over %o and an auxiliary prime over %o.\n", #L`GenusComp`CriticalPrimes, [Minimum(a): a in L`GenusComp`CriticalPrimes] , Minimum(good); 

  // A callback function for the neighbors algorithm which keeps only those lattices not isometric to L, and stops once one non-isometric lattice has been found:
  call:= func< LL | not iso, true where iso:= IsIsometric(L, LL) >;
  call_one:= func< LL | not iso, iso where iso:= IsIsometric(L, LL) >;
  
  // One neighbour each at all the critical primes
  // (calculating a neighbour here *will* change the spinor-"+"-genus, but the resulting lattice might lie in the same spinor genus.)
  found := false;
  for p in L`GenusComp`CriticalPrimes do 
    // instead of the IteratedNeighbours function, we can use Neighbours because we just need to check
    // whether there is any neighbour at all which is not isometric to L:
    found := found or (#Neighbours(L, p: CallBack:= call_one, AutoOrbits:=false) gt 0);
    if not found then 
      // this can only occur if all automorphisms of L have determinant 1 (hence the dimension must be even):
      assert &and[IsOne(Determinant(x)): x in Generators(AutomorphismGroup(L))];
    end if;
    if found then break; end if;
  end for;

  // All the neighbours at our favorite auxiliary low-norm prime:
  // (calculating a neighbour here might or might not change the spinor genus. Depending on the structure of the RCG, 
  // there will be one or two spinor genera covered if we calculate all the neighbours here.)
  if not found then 
    found := found or (#Neighbours(L, good: CallBack:=call, AutoOrbits:=UseAuto) gt 0); 
  end if;
  if (not found) then
    if GF(2)!Rank(L) eq 0 then assert #L`GenusComp`RCG_FactorGroup in {1,2};
    else assert #L`GenusComp`RCG_FactorGroup in {1}; end if;
  end if;
  return not found;  
end intrinsic;
