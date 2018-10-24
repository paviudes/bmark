
// See Delfosse, Zemor - 2010 for this construction


//=============================================
//				Construction of the group G
//=============================================


//Compute the minimal polynomial of 2\cos(2\pi/n) with coef in Z
h := function(n)
	Z<x> := PolynomialRing(IntegerRing());
	C<z> := CyclotomicField(n);
	f := MinimalPolynomial(z + z^-1);
	return(f);
end function;


//Compute the normalized Chebyshev polynomial with coef in Z
P := function(k)
	Z<x> := PolynomialRing(IntegerRing());
	S := Coefficients(ChebyshevFirst(k));
	n := #S-1;
	
	//normalization
	for i := 0 to n do
	   temp := 2^(1-i)*S[i+1];
	   Remove(~S, i+1);
	   Insert(~S, i+1, temp);
	end for;

	P := elt< Z	| S>;
	return(P);
end function;


//construct the generating set S
HypGenerators := function(m, p)
	//construction of the coef ring
	F<x> := PolynomialRing(GF(p));
	H := F! h(2*m*m);
	I := ideal< F | H >;
	Ring := quo< F | I >;

	//construction of S
	Pm := Ring!(F!P(m));
	y := Matrix(Ring, 3, 3, [Pm^2-1, 0, Pm, Pm, 1, 0, -Pm, 0, -1]);
	z := Matrix(Ring, 3, 3, [-1,-Pm,0, Pm, Pm^2-1, 0, Pm, Pm^2, 1]);
	Generators := {@ y, z, y*z @};
	return(Generators);
end function;



//==============================================================
//		Construction of the Cosets
//==============================================================


//Construct the list of all the elements of the group generated by S
Group := function(S)
	Group := S;
	card := 0;

	while card ne #Group do
		new := card;
	   card := #Group;
	   for i := new+1 to card do

			//construction of the neighbours of Group[i]
	      for g := 1 to 3 do
				v := Group[i]*S[g];
				j := Index(Group, v);
				if (j eq 0) then
					Include(~Group, v);
					j:= #Group;
				end if;
	      end for;

	   end for;
	end while;

	return(Group);
end function;


//Construct cosets defining the vertices
VertexCosets := function(Group, S);
	V := {@ @};
	m := Order(S[1]);
	for i := 1 to #Group do
		g := Group[i];
		coset := {@ i @};
		for j := 1 to m-1 do
			g := g*S[1];
			Include(~coset, Index(Group, g));
		end for;
		Include(~V, coset);
	end for;
	return(V);
end function;


/*
//test
m := 5; p := 5;
S := HypGenerators(m, p);
V := VertexCosets(Group(S), S);
#V;
*/


//Construct cosets defining edges
EdgeCosets := function(Group, S);

	E := {@ @};
	m := Order(S[1]);
	for i := 1 to #Group do
		g := Group[i];
		coset := {@ i @};
		Include(~coset, Index(Group, g*S[3]));
		Include(~E, coset);
	end for;
	return(E);
end function;

/*
//test
m := 5; p := 5;
S := HypGenerators(m, p);
E := EdgeCosets(Group(S), S);
#E;
*/


//construct cosets defining faces
FaceCosets := function(Group, S)
	F := {@ @};
	m := Order(S[2]);
	for i := 1 to #Group do
		g := Group[i];
		coset := {@ i @};
		for j := 1 to m-1 do
			g := g*S[2];
			Include(~coset, Index(Group, g));
		end for;
		Include(~F, coset);
	end for;
	return(F);
end function;

/*
//test
m := 5; p := 5;
S := HypGenerators(m, p);
F := FaceCosets(Group(S), S);
#F;
*/


//=====================================================
//		Construction of the edge set and the face set
//=====================================================


//return the edge set from the vertices and edges cosets
SelfDualEdges := function(V_cosets, E_cosets)
	E := {@ @};
 	for i := 1 to #E_cosets do
		ei := E_cosets[i];
		v1 := {@ A : A in V_cosets | ei[1] in A @}[1];
		v2 := {@ A : A in V_cosets | ei[2] in A @}[1];
		ei := { Index(V_cosets, v1), Index(V_cosets, v2) };
		Include(~E, ei);
	end for;
	return(E);
end function;

/*
//test
m :=12; p := 3;
S := HypGenerators(m, p);
V_cosets := VertexCosets(Group(S), S);
#V_cosets;
E_cosets := EdgeCosets(Group(S), S);
#E_cosets;
E := SelfDualEdges(V_cosets, E_cosets);
E;
*/


//return the face set from the faces and edges cosets and the degree m
//faces are described by m indices of edges
SelfDualFaces := function(F_cosets, E_cosets, m)
	//construction of the dual edge set
	E_dual := {@ @};
 	for i := 1 to #E_cosets do
		ei := E_cosets[i];
		f1 := {@ A : A in F_cosets | ei[1] in A @}[1];
		f2 := {@ A : A in F_cosets | ei[2] in A @}[1];
		ei := {@ Index(F_cosets, f1), Index(F_cosets, f2) @};
		Include(~E_dual, ei);
	end for;
	if (#E_dual ne #E_cosets) then
		print "denenerated case in SelfDualFaces: choose a larger prime number p";
	end if;

	//construction of the matrix of faces
	F_matrix := ZeroMatrix(Integers(), #F_cosets, m);
	for i := 1 to #E_dual do
		ei := E_dual[i];
		f1 := ei[1];
		j := 1;
		while (F_matrix[f1, j] ne 0) do
			j := j+1;
		end while;
		F_matrix[f1, j] := i;
		f2 := ei[2];
		j := 1;
		while (F_matrix[f2, j] ne 0) do
			j := j+1;
		end while;
		F_matrix[f2, j] := i;
	end for;

	//transformation of the matrix into a set of faces
	F := {@ @};
	for i := 1 to #F_cosets do
		fj := {@ @};
		for j := 1 to m do
			Include(~fj, F_matrix[i, j]);
		end for;
		Include(~F, fj);
	end for;

	return(F);
end function;

/*
//test
m :=6; p := 5;
S := HypGenerators(m, p);
E_cosets := EdgeCosets(Group(S), S);
#E_cosets;
F_cosets := FaceCosets(Group(S), S);
#F_cosets;
F := SelfDualFaces(F_cosets, E_cosets, Order(S[1]));
#F;
*/



//=================================================
//		Write the graph in graph.txt and dual.txt
//=================================================

//The graph has the following format in these files:
//	* v vertices indexed by integers from 0 to v-1
// 	* e edges indexed by integers from 0 to e-1
//	* f faces indexed by integers from 0 to f-1
//	* an egde is a pair of 2 differents vertices, i.e. 2 integers
//	* a face is a m-tuple of edges, i.e. m integers
//Be careful in Magma indices start from 1

//The format of edge.txt is 
//v e
//v1 v2
//v3 v4
//...
//where
//	* v is the number of vertices
//  * e is the number of edges
//	* the i-th line contains the 2 integers represemting the 2 vertices of the (i-1)-th edge

//The format of face.txt is 
//f
//m1 e11 e12 ... e1m1
//m2 e21 e22 ... e2m2
//...
//where
//  * f is the number of faces
//	* the i-th line contains mi+1 integers describing the face indexed by i-1
//	* the first number mi is the size of the face
// 	* it is followed by mi integers in increasing order which are the indices of 
//		the mi edges of this face


//print edges in edges.txt
PrintEdges := function(E, v, filename)
	file := Open(filename, "w");

	//write the type v e of the graph
	nbvertices := IntegerToString(v);
	nbedges := IntegerToString(#E);
	size := &cat [nbvertices, " ", nbedges];
	Puts(file, size);

	//write all the edges
	for k := 1 to #E do
		Ek := SetToIndexedSet(E[k]);
		vi := IntegerToString(Ek[1]-1);
		vj := IntegerToString(Ek[2]-1);
		edgeij := &cat [vi, " ", vj];
		Puts(file, edgeij);
	end for;

	return(filename);
end function;


//print faces in faces.txt
PrintFaces := function(F, filename)
	file := Open(filename, "w");

	//write the type f  of the graph
	nbfaces := IntegerToString(#F);
	Puts(file, nbfaces);

	//write all the edges
	for k := 1 to #F do
		Facek := IntegerToString(#F[k]);
		for i := 1 to #F[k] do
			vi := IntegerToString(F[k][i]-1);
			Facek := &cat [Facek, " ", vi];
		end for;
		Puts(file, Facek);
	end for;

	return(filename);
end function;


//==================================================
//		TESTS					 
//==================================================


m := 6;
p := 79;

//construct the graph type string
i1 := IntegerToString(m);
i2 := IntegerToString(p);
family := "selfdual";
graph_type := &cat [i1, "_" , i2, "_", family, "_"];

//construct the graph
S := HypGenerators(m, p);
Gr := Group(S);
V_cosets := VertexCosets(Gr, S);
v := #V_cosets; v;
E_cosets := EdgeCosets(Gr, S);
#E_cosets;
F_cosets := FaceCosets(Gr, S);
#F_cosets;

Vertices := SetToIndexedSet({1..v});
Edges := SelfDualEdges(V_cosets, E_cosets);
Faces := SelfDualFaces(F_cosets, E_cosets, m);
//G := Graph< Vertices | IndexedSetToSet(Edges) >;

PrintEdges(Edges, v, graph_type cat "edges.txt");
PrintFaces(Faces, graph_type cat "faces.txt");




