//freeze;

declare type ProcPL;

declare attributes ProcPL:
  a,
  v,
  dim,
  depth;

declare attributes RngOrd: 
  SquareClasses;

intrinsic Print(PL::ProcPL)
{internal}
  printf "A projective line process for %m", Parent(PL`v);
end intrinsic;

intrinsic ProjectiveLineProcess(V::ModTupFld[FldFin]) -> ProcPL
{Creates a projective line process for V}
  PL:= New(ProcPL);
  PL`a:= PrimitiveElement(BaseField(V));
  PL`v:= V ! 0;
  PL`dim:= Dimension(V);
  PL`depth:= PL`dim+1;
  return PL;
end intrinsic;

intrinsic ProjectiveLineProcess(k::FldFin, n::RngIntElt) -> ProcPL
{Creates a projective line process for k^n}
  requirege n, 1;
  return ProjectiveLineProcess(VectorSpace(k, n));
end intrinsic;

intrinsic '#'(PL::ProcPL) -> RngIntElt
{The number of lines in PL}
  q:= #Parent(PL`a);
  return (q^PL`dim-1) div (q-1);
end intrinsic;

intrinsic Next(PL::ProcPL) -> ModTupFldElt
{The next element in the process. Returns the zero vector if no more elements left}
  if PL`depth ne 0 then
    i:= PL`dim;
    while true do
      if i eq PL`depth then
        PL`v[i]:= 0; i -:= 1;
      elif i lt PL`depth then
        PL`depth:= i;
	if i ge 1 then PL`v[i]:= 1; end if;
	break;
      elif PL`v[i] eq 0 then PL`v[i]:= 1; break;
      else
        PL`v[i] *:= PL`a;
	if PL`v[i] eq 1 then PL`v[i]:= 0; i -:= 1; else break; end if;
      end if;
    end while;
  end if;
  return PL`v;
end intrinsic;

intrinsic MyOrbitsOfSpaces(G::GrpMat[FldFin], k::RngIntElt) -> []
{Orbit representatives of the k-dimensional subspaces of the natural vector space for G}
  n:= Degree(G);
  requirerange k, 0, n;

  F:= BaseRing(G);
  if IsPrimeField(F) then
    O:= OrbitsOfSpaces(G, k);
    return [ BasisMatrix(o[2]): o in O ];
//  elif 2*k gt n then
//    O:= MyOrbitsOfSpaces(G, n-k);
//    return [ KernelMatrix(Transpose(b)) : b in O ];
  end if;

  V:= VectorSpace(F, n);
  M:= MatrixRing(F,n);
  One:= M ! 1;
  L:= {@@};
  for S in Subsets({1..n}, k) do
    T:= Sort(Setseq(S));
    A:= M ! 0;
    for i in [1..k] do A[i, T[i]]:= 1; end for;
    Cols:= [ Sort(Setseq({T[i]+1..n} diff S)) : i in [1..k] ];
    c:= Integers() ! (k*(n-k/2+1/2)) - &+T;
    if c eq 0 then 
      Include(~L, RowSpace(A));
    else
      C:= CartesianProduct([F^^c]);
      for x in C do
        B:= A; r:= 1;
        for i in [1..k], j in Cols[i] do
          B[i,j]:= x[r]; r +:= 1;
        end for; 
        Include(~L, Rowspace(B));
      end for;
    end if;
  end for;

  if forall{ i : i in [1..Ngens(G)] | IsScalar(G.i) } then
    return [ BasisMatrix(x) : x in L ];
  end if;

  O:= [];
  while not IsEmpty(L) do
    o:= Rep(L);
    Append(~O, BasisMatrix(o));
    L:= L diff Orbit(G, o);
  end while;
  return O;
end intrinsic;

// Computing with K_p^* / (K_p^*)^2.

// Helper function to ensure x is a unit.
function forceunit(x)
  error if x eq 0, "element must be non-zero";
  return x;
end function;

// Helper function to find y such that xy^2 is integral and Valuation(xy^2, p) = 0
function SquareRepNice(x, p, piinv)
  x *:= Denominator(x)^2;
  v:= Valuation(x, p);
  assert v ge 0 and IsEven(v);
  if v ne 0 then x *:= piinv^v; end if;
  return Order(p) ! x;
end function;

// todo: Abelian group would be better....
intrinsic LocalMultiplicativeGroupModSquares(p::RngOrdIdl) -> ModFld, Map
{Given a prime ideal over some number field K, this returns a vectorspace
 over GF(2) isomorphic to K_p^/(K_p^*)^2 and a map representing the isomorphism}
  require IsPrime(p): "The ideal must be prime";
  R:= Order(p);
  K:= NumberField(R);
  if Minimum(p) ne 2 then
    pi:= PrimitiveElement(p);
    k, h:= ResidueClassField(p);
    e:= Nonsquare(k) @@ h;
    V:= VectorSpace(GF(2), 2);
    m:= map< V -> K |
      x :-> (x[1] eq 0 select 1 else e) * (x[2] eq 0 select 1 else pi),
      y :-> [ IsSquare((y/pi^v) @ h) select 0 else 1, v] where v:= Valuation(forceunit(y), p)
    >;
    return V, m;
  else
    if not assigned R`SquareClasses then R`SquareClasses:= AssociativeArray(); end if;
    ok, m:= IsDefined(R`SquareClasses, p);
    if not ok then
      if not IsAbsoluteOrder(R) then
        p:= AbsoluteOrder(R) !! p;
      end if;
      pi:= PrimitiveElement(p);
      e:= RamificationDegree(p);
      dim:= Valuation(Norm(p), 2)*e+2;
      V:= VectorSpace(GF(2), dim);
      I:= p^(2*e+1);
      Q, h:= quo< Order(I) | I >;
      U, g:= UnitGroup(Q);	// fails for non-absolute orders
      M, i:= quo< U | 2*U >;
      assert #M eq 2^(dim-1);
      piinv:= WeakApproximation([p], [-1]);
      m:= map< V -> K | x :-> ((M ! ChangeUniverse(Eltseq(x)[1..dim-1], Integers())) @@ i @ g @@ h) * (x[dim] eq 0 select 1 else pi),
                        y :-> Append(Eltseq(SquareRepNice(y * pi^v, p, piinv) @@ g @ i), v) where v:= Valuation(y, p) mod 2 >;
      R`SquareClasses[p]:= m;
    end if;
    return Domain(m), m;
  end if;
end intrinsic;

intrinsic LocalMultiplicativeGroupModSquares(p::RngIntElt) -> ModFld, Map
{"} //"
  require IsPrime(p: Proof:= false) : "The argument must be prime";
  if p ne 2 then
    e:= Integers() ! Nonsquare(GF(p));
    V:= VectorSpace(GF(2), 2);
    m:= map< V -> Rationals() |
      x :-> (x[1] eq 0 select 1 else e) * (x[2] eq 0 select 1 else p),
      y :-> [ KroneckerSymbol( Integers(p) ! (y/p^v) ) eq 1 select 0 else 1, v] where v:= Valuation(forceunit(y), p)
    >;
  else
    V:= VectorSpace(GF(2), 3);
    R:= Integers(8);
    m:= map< V -> Rationals() | 
      x :-> (x[1] eq 0 select 1 else 3) * (x[2] eq 0 select 1 else -1) * (x[3] eq 0 select 1 else 2),
      y :-> case< R ! (y/p^v) | R!1: [0,0], R!3: [1,0], R!-1:[0,1], default: [1,1] > cat [ v ] where v:= Valuation(forceunit(y), p)
    >;
  end if;
  return V, m;
end intrinsic;

intrinsic UnitSquareClassReps(p::RngOrdIdl) -> []
{Given a prime ideal of some number ring return a system
 of representatives of the units R_p^* modulo squares}
  V, h:= LocalMultiplicativeGroupModSquares(p);
  U:= sub< V | [V.i: i in [1..Dimension(V)-1]] >;
  return [ u @ h : u in U ];
end intrinsic;

intrinsic NiceUnitSquareClassRepresentative(u::RngElt, p::RngOrdIdl) -> RngElt
{Given a unit at the prime p, find a local unit v in the same square class such that v-1 generates the quadratic defect of u}
  require IsPrime(p) : "The second argument must be prime";
  R:= Order(p);
  K:= NumberField(R);
  ok, u:= IsCoercible(K, u);
  require ok : "Incompatible arguments";
  require Valuation(u, p) eq 0 : "The element is not a local unit";

  k,h:= ResidueClassField(p);
  ok, s := IsSquare(h(u));
  assert ok;
  u := u / (s @@ h)^2;
  assert IsOne(h(u));
  e:= Valuation(R ! 2, p);
  pi:= PrimitiveElement(p);

  val := Valuation(u-1, p);
  while val lt 2*e do
    if IsOdd(val) then return u; end if;
    ok, s:= IsSquare(h((u - 1) / pi^val)); assert ok;
    u /:= (1+(s@@h)*pi^(val div 2))^2;
    val := Valuation(u-1, p);
  end while;

  return val eq 2*e and IsIrreducible(Polynomial([k | ((u-1)/4) @ h, 1, 1])) select u else R ! 1;
end intrinsic;

// End of unit group computations


// The even HilbertSymbol code

// We first start with our implementation of Alg. 6.2 -- 6.5 of
// John Voight, Characterizing quaternion rings over an arbitrary base, J. Reine Angew. Math. 657 (2011), 113-134.
function SolveCongruence(a,b,p)
  //assert Valuation(a,p) eq 0;
  //assert Valuation(b,p) eq 1;
  pi:= PrimitiveElement(p);
  k, h:= ResidueClassField(p);
  y:= 1/(SquareRoot(a @ h) @@ h); z:= 0;
  ee:= 2*RamificationIndex(p);
  told:= -1;
  while true do
    N:= 1 - a*y^2 - b*z^2;
    t:= Valuation(N, p); assert t gt told; told:= t;
    if t ge ee then return y,z; end if;
    w:= pi^(t div 2);
    if IsEven(t) then
      y +:= SquareRoot((N/(a*w^2)) @ h) @@ h * w;
    else
      z +:= SquareRoot((N/(b*w^2)) @ h) @@ h * w;
    end if;
  end while;
end function;

// Decide if a*x^2=1 mod p^e
// This is a variant of O'Meara ???. The code assumes that 0 <= e <= 2*RamIdx(p)
function CanLiftSquare(a,p,e)
  assert Valuation(a,p) eq 0;
  if Valuation(a-1, p) ge e then return true, Order(p) ! 1, e; end if;
  aold:= a;
  pi:= PrimitiveElement(p);
  k, h:= ResidueClassField(p);
  x:= SquareRoot( (a @ h)^-1 ) @@ h;
  a *:= x^2;
  t:= Valuation(a-1, p);
  while t lt e and IsEven(t) do
    s:= SquareRoot(h((a-1)/pi^t));
    s:= 1+(s@@h)*pi^(t div 2);
    a /:= s^2;
    x := x/s;
    tt:= Valuation(a-1, p);
    assert tt gt t;
    t:= tt;
  end while;
  t:= Min(t,e);
  assert Valuation(aold*x^2-1, p) ge t;
  return t eq e, x, t;
end function;

function SolveCongruence2(a,b,p)
  if Valuation(b, p) eq 1 then 
    y,z:= SolveCongruence(a,b,p);
    return y,z,0;
  end if;
  
  e:= RamificationIndex(p);
  pi:= PrimitiveElement(p);
  oka, a0, t:= CanLiftSquare(a, p, e);
  if oka then
    okb, b0, t:= CanLiftSquare(b, p, e);
    if okb then return a0, b0, a0*b0; end if;
    bt:= (b-(1/b0)^2) /pi^t;
    y, z:= SolveCongruence(b, -pi*bt/a, p);
    w:= pi^(t div 2);
    return w*b0/z, b0, y*w*b0/z;
  else
    at:= (a-(1/a0)^2) /pi^t;
    y, z:= SolveCongruence(a, -pi*at/b, p);
    w:= pi^(t div 2);
    return a0, w*a0/z, y*w*a0/z;
  end if;
end function;

normalize:= function(a, p)
  vala:= Valuation(a, p);
  n:= vala div 2;
  if n ne 0 then a /:= PrimitiveElement(p)^(2*n); vala -:= 2*n; end if;
  assert Valuation(a, p) eq vala and vala in {0,1};
  return a, vala;
end function;

//function MyEvenHilbertSymbol(a,b,p)
intrinsic MyEvenHilbertSymbol(a::FldElt,b::FldElt,p::RngOrdIdl) -> RngIntElt
{}
  a, vala:= normalize(a, p);
  b, valb:= normalize(b, p);
  if vala eq 1 then
    if valb eq 1 then
      a:= (-a*b) / PrimitiveElement(p)^2;
    else
      x:= a; a:= b; b:= x;
    end if;
  end if;

  // lift solution to precision 2e
  y,z,w:= SolveCongruence2(a,b,p);
  nrm:= 1-a*y^2-b*z^2+a*b*w^2;
  assert Valuation(nrm, p) ge 2*RamificationIndex(p);
  tmp:= (b*z)^2*a + (a*y)^2*b;

  if tmp eq 0 or IsEven(Valuation(tmp, p)) then 
    res:= 1;
  else
    k,h:= ResidueClassField(p);
    res:= not IsIrreducible(Polynomial([k| (nrm/4)@h,-1,1])) select 1 else -1; 
  end if;
  //assert res eq HilbertSymbol(a,b,p);
  return res;
end intrinsic;
//end function;

intrinsic MyHilbertSymbol(a::FldAlgElt, b::FldAlgElt, p::RngOrdIdl) -> RngIntElt
{The Hilbert symbol of a and b at p}
  if IsOne(a) or IsOne(b) then return 1; end if;
  require not IsZero(a) and not IsZero(b): "The elements must be non-zero";
  require IsPrime(p): "The ideal must be prime";
  return Minimum(p) ne 2 select HilbertSymbol(a,b,p) else MyEvenHilbertSymbol(a,b,p);
end intrinsic;


// Working with ideal class groups

function CGPrimes(I, S, Generators, CoprimeTo, Minimal, Quotient)
  R:= Order(I);
  if not IsMaximal(R) or not IsAbsoluteOrder(R) then return false, "The order must be absolute and maximal"; end if;
  r1, r2:= Signature(NumberField(R));
  if not IsEmpty(S) then
    T:= Type(Universe(S));
    if T eq PlcNum then
      X:= S; S:= [];
      for s in X do
        ok, i:= IsInfinite(s); 
        if not ok or not IsReal(s) then return false, "The places must be real"; end if;
        Append(~S, i);
      end for;
    elif (T eq RngInt) and Minimum(S) ge 1 and Maximum(S) le r1 then
      ;
    elif (Universe(S) cmpeq PowerStructure(Infty)) and AbsoluteDegree(R) eq 1 then 
      S:= [1];
    else
      return false, "Wrong infinite places";
    end if;
    S:= Sort(S);
  end if;
  if not IsIntegral(I) then
    return false, "The ray is not integral";
  end if;

  if Type(Quotient) eq SetEnum then Quotient:= Setseq(Quotient); end if;
  if #Quotient ne 0 then
    T:= Type(Quotient[1]);
    if T in {RngIntElt, FldRatElt} then
      Quotient:= [ ideal< R | x> : x in Quotient ];
    elif T eq RngInt then
      Quotient:= [ ideal< R | Generator(x) > : x in Quotient ];
    elif not ISA(T, RngOrdFracIdl) or Order(Quotient[1]) cmpne R then
      return false, "Incompatible user defined generators";
    end if;
    if exists{ u: u in Quotient | not IsOne(I+u) } then 
      return false, "The user defined generators must be comprime to the ray";
    end if;
  end if;

  if CoprimeTo cmpne 1 then
    if ISA(Type(CoprimeTo), RngElt) then
      ok, CoprimeTo:= IsCoercible(FieldOfFractions(R), CoprimeTo);
      if not ok then return ok, "IscoprimeTo must be an ideal or field element"; end if;
    end if;
  end if;

  L:= [ PowerIdeal(R) | ];
  C, h:= RayClassGroup(I, S);
  C0:= sub< C | {u @@ h : u in Quotient} >;
  if #C0 ne 1 then
    C, hh:= quo< C | C0 >;
    h:= hh^-1 * h;
  end if;

  if Generators then
    n:= #AbelianInvariants(C);	// We will end up with this number of gens.
    U:= sub< C | >;
    p:= 1; i:= 1; D:= [];
    while n ne 0 do
      if i ge #D then p:= NextPrime(p); D:= Decomposition(R, p); i:= 1; end if;
      if (CoprimeTo cmpeq 1 or Valuation( CoprimeTo, D[i,1] ) eq 0) and IsOne(I+D[i,1]) then
        g:= D[i,1] @@ h;
        if g notin U then
          nn:= #AbelianInvariants( quo< C | U, g > );
          if nn ne n then
            assert nn eq n-1;
            n:= nn;
            Append(~L, D[i,1]);
            U:= sub< C | U, g >; 
          end if;
        end if;
      end if;
      i +:= 1;
    end while;
  else
    U:= {@ @};
    p:= 1; i:= 1; D:= [];
    while #U ne #C do
      if i ge #D then p:= NextPrime(p); D:= Decomposition(R, p); i:= 1; end if;
      if (CoprimeTo cmpeq 1 or Valuation( CoprimeTo, D[i,1] ) eq 0) and IsOne(I+D[i,1]) then
        g:= D[i,1] @@ h;
        if g notin U then
          Append(~L, D[i,1]);
          Include(~U, g);
        end if;
      end if;
      i +:= 1;
    end while;

    if Minimal then
      B:= Max([Norm(p): p in L]);
      P:= PrimesUpTo(B-1, NumberField(R): coprime_to:= CoprimeTo * I);
    
      for p in P do
        g:= p @@ h;
        i:= Index(U, g);
        if Norm(p) le Norm(L[i]) then 
          L[i]:= p;
        end if;
      end for;
    end if;
  end if;

  Norms:= [ Norm(p): p in L ];
  ParallelSort(~Norms, ~L);

  return true, L;
end function;

// Computing prime ideal generators of (ray) class groups

intrinsic ClassGroupPrimeIdealGenerators(I :: RngOrdIdl, S::[] : CoprimeTo:= 1, Quotient:= []) -> []
{Returns prime ideals that generate the ray class group of I and the infinite places in S}
  ok, L:= CGPrimes(I, S, true, CoprimeTo, true, Quotient);
  require ok: L;
  return L;
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(I :: RngOrdIdl : CoprimeTo:= 1) -> []
{"}; //"
  ok, L:= CGPrimes(I, [], true, CoprimeTo, true, Quotient);
  require ok: L;
  return L;
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(I :: RngIntElt, S::[] : CoprimeTo:= 1, Quotient:= []) -> []
{"} // "
  require I ne 0: "The modulus cannot be zero";
  R:= Integers(QNF());
  I:= ideal< R | I >;
  ok, L:= CGPrimes(I, S, true, CoprimeTo, true, Quotient);
  require ok: L;
  return [ Minimum(l): l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(I :: RngIntElt : CoprimeTo:= 1, Quotient:= []) -> []
{"} // "
  require I ne 0: "The modulus cannot be zero";
  R:= Integers(QNF());
  I:= ideal< R | I >;
  ok, L:= CGPrimes(I, [], true, CoprimeTo, true, Quotient);
  require ok: L;
  return [ Minimum(l): l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(R :: RngInt, S::[]: CoprimeTo:= 1, Quotient:= []) -> []
{"} //"
  g:= Generator(R);
  require g ne 0 : "The ideal must not be zero";
  R:= Integers(QNF());
  I:= ideal< R | g >;
  ok, L:= CGPrimes(I, S, true, CoprimeTo, true, Quotient);
  require ok: L;
  return [ ideal< Integers() | Minimum(l) > : l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(R :: RngInt: CoprimeTo:= 1, Quotient:= []) -> []
{"} //"
  g:= Generator(R);
  require g ne 0 : "The ideal must not be zero";
  S:= Integers(QNF());
  I:= ideal< S | g >;
  ok, L:= CGPrimes(I, [], true, CoprimeTo, true, Quotient);
  require ok: L;
  return [ ideal< Integers() | Minimum(l) > : l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(R :: RngOrd, S::[]: CoprimeTo:= 1, Quotient:= []) -> []
{Returns prime ideals that generate the ideal class group of R}
  ok, L:= CGPrimes(1*R, S, true, CoprimeTo, true, Quotient);
  require ok: L;
  return L;
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(R :: RngOrd: CoprimeTo:= 1, Quotient:= []) -> []
{Returns prime ideals that generate the ideal class group of R}
  ok, L:= CGPrimes(1*R, [], true, CoprimeTo, true, Quotient);
  require ok: L;
  return L;
end intrinsic;

// The same with representatives

intrinsic ClassGroupPrimeIdealRepresentatives(I :: RngOrdIdl, S::[] : CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{Returns prime ideals that generate the ray class group of I and the infinite places in S}
  ok, L:= CGPrimes(I,  S, false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return L;
end intrinsic;

intrinsic ClassGroupPrimeIdealRepresentatives(I :: RngOrdIdl : CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{"}; //"
  ok, L:= CGPrimes(I, [], false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return L;
end intrinsic;

intrinsic ClassGroupPrimeIdealRepresentatives(I :: RngIntElt, S::[] : CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{"} //"
  require I ne 0: "The modulus cannot be zero";
  R:= Integers(QNF());
  I:= ideal< R | I >;
  ok, L:= CGPrimes(I,  S, false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return [ Minimum(l): l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealRepresentatives(I :: RngIntElt: CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{"} //"
  require I ne 0: "The modulus cannot be zero";
  R:= Integers(QNF());
  I:= ideal< R | I >;
  ok, L:= CGPrimes(I,  [], false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return [ Minimum(l): l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealRepresentatives(R :: RngInt, S::[]: CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{"} //"
  g:= Generator(R);
  require g ne 0 : "The ideal must not be zero";
  R:= Integers(QNF());
  I:= ideal< R | g >;
  ok, L:= CGPrimes(I, S, false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return [ ideal< Integers() | Minimum(l) > : l in L ];
end intrinsic;

intrinsic ClassGroupPrimeIdealGenerators(R :: RngInt: CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{"} //"
  g:= Generator(R);
  require g ne 0 : "The ideal must not be zero";
  S:= Integers(QNF());
  I:= ideal< S | g >;
  ok, L:= CGPrimes(I, [], false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return [ ideal< Integers() | Minimum(l) > : l in L ];
end intrinsic;


intrinsic ClassGroupPrimeIdealRepresentatives(R :: RngOrd: CoprimeTo:= 1, Minimal:= true, Quotient:= []) -> []
{Returns prime ideals that represent the ideal classes of R}
  ok, L:= CGPrimes(1*R, [], false, CoprimeTo, Minimal, Quotient);
  require ok: L;
  return L;
end intrinsic;


/* Evaluating Dedekind zeta functions at negative integers exactly */

function ExtensionToHeckeCharacter(E)
  assert Degree(E) eq 2;
  K:= BaseField(E);
//  if not IsAbsoluteField(K) then K:= AbsoluteField(K); end if;
  RE:= Integers(E);

  S:= [];
  for i in RealPlaces(K) do
    if #Decomposition(E, i) ne 2 then
      if Type(i) eq Infty then
        idx:= 1;
      else
        ok, idx:= IsInfinite(i); assert ok;
      end if;
      Append(~S, idx);
    end if;
  end for;
  S:= Sort(S);

//  S:= [1..Degree(K)];
  bad:= Type(K) eq FldRat;
  if bad then 
    DE:= Integers( QNF() ) * Discriminant(RE);
  else
    DE:= Discriminant(RE);
  end if;
  P:= ClassGroupPrimeIdealGenerators(DE, S);
  T:= < < p, IsSplit(bad select Minimum(p) else p, RE) select 1 else -1> : p in P >; 
  h:= HeckeCharacter(DE, S, T);
  assert IsPrimitive(h);
  return h;
end function;

function myEval(K, z, Relative)
  if IsOdd(z) then
    k:= 1-z;
    if Type(K) eq FldRat then
      return BernoulliNumber(k)/-k;
    elif Type(K) eq FldQuad then
      d:= Discriminant(Integers(K));
      if d gt 1 then
        return BernoulliNumber(k) * BernoulliNumber(k, KroneckerCharacter(d, Rationals())) / k^2;
      end if;
    end if;
  end if;

  if Relative then
    H:= ExtensionToHeckeCharacter(K);
    L:= LSeries(H);
//    F:= BaseField(K);
//    K:= OptimizedRepresentation(AbsoluteField(K));
//    L:= LSeries(K : Method:= Degree(F) ge 5 select "Direct" else "Default") / LSeries(F);
  else
    L:= LSeries(K);
  end if;

  i:= 0;
  repeat
    if i ge 1 then 
      LSetPrecision(L, 40 + i*20); 
      "increasing precision", i; 
    end if;
    x:= Evaluate(L, z);
    if Type(x) eq FldComElt and Im(x) le 10^-20 then x:= Re(x); end if;
    X:= Type(x) eq FldReElt select { BestApproximation(x, 10^i) : i in [12, 14, 16, 18] } else [];
    i +:= 1;
  until #X eq 1;
  X:= Rep(X);

//  if Relative then
//    assert Abs(Real(Evaluate(LSeries(H), z)) - X) le 10^-10;
//  end if;
  return X;
end function;

intrinsic DedekindZetaExact(K::FldAlg, z::RngIntElt: Relative:= false) -> FldRatElt
{Evaluates the Dedekind zeta function of K at the negative integer z}
  require (Relative and z eq 0) or z lt 0 : "The argument must be a negative integer";
  return myEval(K, z, Relative);
end intrinsic;

intrinsic DedekindZetaExact(K::FldRat, z::RngIntElt: Relative:= false) -> FldRatElt
{"} //"
  require (Relative and z eq 0) or z lt 0 : "The argument must be a negative integer";
  return myEval(K, z, Relative);
end intrinsic;


intrinsic HasInvariantIdealClassRepresentative(I::RngOrdFracIdl) -> BoolElt, FldAlgElt
{Given an ideal in a quadratic extension, check if xI is fixed under the Galois group for some x. Is so, returns true and x.}
   R:= Order(I);
   require IsMaximal(R) and Degree(R) eq 2: "The order of the ideal must be maximal and of degree 2";
   K:= NumberField(R);
   a:= Automorphisms(K)[2];
   J := ideal< R | [ a(x): x in Generators(I) ] >;
   ok, x:= IsPrincipal( J * I^-1 );
   if not ok then return false, _; end if;
   S:= Integers(BaseField(K));
   ok, n:= IsCoercible(S, Norm(x));
   if not ok or not IsUnit(n) then return false, _; end if;
   ok, o:= NormEquation(R, n : Solutions:= 1, Exact);
   if not ok then return false, _; end if;
   y:= Hilbert90( (K ! x)/(K ! o[1]), a );
   assert J eq ideal< R | [ a(x): x in Generators(J) ] > where  J:= y * I;
   return true, y;
end intrinsic;

//This intrinsic can be called from the debugger to export data
declare attributes RngInt: User1;

intrinsic SaveMyFriend(X) -> RngIntElt
{internal}
  ZZ:= Integers();
  ZZ`User1:= X;
  return 0;
end intrinsic;

intrinsic RestoreMyFriend() -> Any
{internal}
  ZZ:= Integers();
  require assigned ZZ`User1: "No data stored. Call SaveMyFriend() before.";
  return ZZ`User1;
end intrinsic;
