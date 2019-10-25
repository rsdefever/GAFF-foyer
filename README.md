# GAFF-foyer

This repository provides a foyer (https://github.com/mosdef-hub/foyer) compatible XML file for the Generalized Amber Force Field (GAFF). 

### Currently in beta. Please always verify the choice of atomtypes.

## Known issues

GAFF has 'paired' atomtypes for conjugated systems (e.g., cc/cd, ce/cf, nc/nd, etc.). The reasoning for this is explained in the GAFF and Antechamber papers (see links below). This implementation does not currently support identifying the latter of the paired types. In other words, all conjugated types will be identified as the first of the pair (e.g., all 'cc' and 'cd' atomtypes will be identified as 'cc').

## Source of parameters

Parameters for GAFF were taken from AmberTools19. The md5sum for the AmberTools19.tar.bz2 file was afffe8a5473a0bd143b98f0396f52f0f. The md5sum for the gaff.dat file in this repo is 2635454817b8f9139d87d4bb397c47b3.

## Implementing GAFF support in foyer

This document describes the development of foyer SMARTS strings for the Generalized Amber Force Field (GAFF). A few related links are provided below:

- [MoSDeF](https://mosdef.org)
- [Foyer GitHub](https://github.com/mosdef-hub/foyer)
- [Foyer paper](https://doi.org/10.1016/j.commatsci.2019.05.026)
- [GAFF paper](https://doi.org/10.1002/jcc.20035)
- [Antechamber paper](https://doi.org/10.1016/j.jmgm.2005.12.005)

The purpose of this document is threefold. (1) to explain the logic used to develop foyer SMARTS strings for GAFF; (2) to document the test cases used to validate the SMARTS strings; and (3) to provide guidance for developing/extending foyer SMARTS strings for GAFF or other forcefields. This document assumes the reader has some familiarity with the concept of atomtyping and the use of SMARTS strings in Foyer to accomplish this end. The reader is directed to the Foyer paper (link above) for further background.

Foyer uses SMARTS strings to create a general atomtyping procedure that can be applied for various forcefields. SMARTS strings are a sequence of characters used to describe the local chemical environment around some atom. Details can be found [here](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html). Foyer does not support the entire SMARTS language. An up-to-date list of supported features is [here](https://github.com/mosdef-hub/foyer/issues/63). In Foyer, the SMARTS strings translate the description of a given atomtype provided by the forcefield authors to a machine-readable graph of a local chemical environment matching the human-readable description. Therefore, implementing new forcefields in Foyer requires developing SMARTS strings that can describe the various atomtypes provided. 

I now describe the development of SMARTS strings for each atomtype in GAFF. The `name` (type name), `class`, `element`, `def` (definition), `overrides`, `desc` (description) fields are taken directly from the forcefield XML file. The atomtype `mass` and `doi` fields have been excluded for brevity. The SMARTS string is located in the `def` field.

I will present the atomtypes in order of increasing complexity of the SMARTS definitions. I first focus on the 35 basic atomtypes in GAFF. 

### Basic atomtypes

Atomtypes for which there is only a single atomtype per element are trivial to define. In the case of GAFF, that includes fluorine, chlorine, bromine, and iodine:

	<Type name="f" class="f" element="F" def="F" desc="Fluorine"/>
	<Type name="cl" class="cl" element="Cl" def="Cl" desc="Chlorine"/>
	<Type name="br" class="br" element="Br" def="Br" desc="Bromine"/>
	<Type name="i" class="i" element="I" def="I" desc="Iodine"/>

Since they are unique for a given element these definitions will always be assigned correctly and will never override nor be overriden by any other atomtypes.

Next, we consider atomtypes that have unique bonding patterns with neighboring elements. In GAFF, those include the following: 

	<Type name="o" class="o" element="O" def="[O;X1]" desc="Oxygen with one connected atom"/>
	<Type name="oh" class="oh" element="O" def="[O;X2]H" desc="Oxygen in hydroxyl group"/>
	<Type name="os" class="os" element="O" def="[O;X2]([!H])[!H]" desc="Ether and ester oxygen"/>

Once again these atomtypes will never be overriden or override other elements; all cases where oxygen or sulfur are bonded to a single other element are type `o` or `s`, respectively. Though the hydroxyl and ether oxygen atoms are both have two connections, they are easily differentiated by the elements with which they are connected. The hydrogen atomtypes are similarly simple to define: 

	<Type name="hc" class="hc" element="H" def="H[C;X4]" desc="H bonded to aliphatic carbon without electrwd. group" />
	<Type name="ha" class="ha" element="H" def="H[C;!X4]" desc="H bonded to aromatic carbon" />
	<Type name="hn" class="hn" element="H" def="HN" desc="H bonded to nitrogen atoms"/>
	<Type name="ho" class="ho" element="H" def="HO" desc="Hydroxyl group" />
	<Type name="hp" class="hp" element="H" def="HP" desc="H bonded to phosphate"/>
	<Type name="hs" class="hs" element="H" def="HS" desc="Hydrogen bonded to sulphur">

Here "aromatic carbon" means any non-sp3 carbon. This is determined from the atomtype definitions distributed with `Antechamber`, in `ATOMTYPE_GFF.DEF`. 

The remainder of the atomtypes begin to either require more complex definitions or are involved in overrides. GAFF often defines atomtypes by hydridization. For carbon, it is easy to infer hybrization from the number of bonds. We can thus define the following: 

	<Type name="c1" class="c1" element="C" def="[C;X2]" desc="Sp C"/>
	<Type name="c2" class="c2" element="C" def="[C;X3]" desc="Sp2 C"/>
	<Type name="c3" class="c3" element="C" def="[C;X4]" desc="Sp3 C"/>

These atomtypes are general and will be overridden by more specific types in certain cases. Once such example follows:

	<Type name="c" class="c" element="C" mass="12.01" def="[C;X3][O&X1,S&X1]" overrides="c2,ca,cc_r5,cc_r6,ce,cv" desc="Sp2 C in C=O or C=S"/>

Note that although the `c` type overrides the `c2`, the order of the definitions in the forcefield XML is irrelevant. Also note the use of the `&` symbol as a high priority `AND` statement that will be evaluated before the `OR` (`,`) statement.

We approach the nitrogen atomtypes similarly. An additional challenge with certain nitrogen atomtypes is that the hydridization cannot be directly mapped onto the number of bonds (i.e., a nitrogen could be `sp2` or `sp3` hydridized and be bonded with three other atoms). The initial definitions however, remain simple.

    <Type name="n1" class="n1" element="N" def="[N;X1]" desc="Sp N"/>
    <Type name="n2" class="n2" element="N" def="[N;X2]" desc="aliphatic Sp2 N with two connected atoms"/>
    <Type name="n3" class="n3" element="N" def="[N;X3]" desc="Sp3 N with three connected atoms"/>
    <Type name="n4" class="n4" element="N" def="[N;X4]" desc="Sp3 N with four connected atoms"/>

Once again these types will be frequently overridden by more specific definitions. Examples are the nitrogen in nitro groups and the nitrogen in amides: 

	<Type name="no" class="no" element="N" def="[N;X3]([O;X1])[O;X1]" overrides="n3,na,na_r5,na_r6,n3,nh_r5,nh_r6,nh" desc="Nitro N">
	<Type name="n" class="n" element="N" def="[N;X3][C;X3][O&X1,S&X1]" overrides="n3,nh,nh_r5,nh_r6,na,na_r5,na_r6" desc="Sp2 nitrogen in amide groups"/>

These atomtypes are well-suited to SMARTS strings since the local chemical environment around the atomtypes can be mapped onto elements and their connectivity.

Next we consider the sulfur atomtypes:

    <Type name="s2" class="s2" element="S" def="[S;X2][!H]" desc="sp2 sulfur (P=S,C=S,etc)"/>
    <Type name="s4" class="s4" element="S" def="[S;X3]" desc="S with three connected atoms"/>
    <Type name="s6" class="s6" element="S" def="[S;X4]" desc="S with four connected atoms"/>
    <Type name="sh" class="sh" element="S" def="[S;X2]H" overrides="s2" desc="Sp3 S connected with hydrogen"/>
    <Type name="ss" class="ss" element="S" def="[S;X2]([C,S,O])[C,S,O]" overrides="s2" desc="Sp3 S in thio-ester and thio-ether"/>

Followed by phosphorus:
    
    <Type name="p2" class="p2" element="P" def="[P;X2]" desc="Phosphate with two connected atoms"/>
    <Type name="p3" class="p3" element="P" def="[P;X3]" desc="Phosphate with three connected atoms, such as PH3"/>
    <Type name="p4" class="p4" element="P" def="[P;X3][%o,%s,%n2]" overrides="p3" desc="Phosphate with three connected atoms, such as O=P(CH3)2"/>
    <Type name="p5" class="p5" element="P" def="[P;X4]" desc="Phosphate with four connected atoms, such as O=P(OH)3"/>

Now we will look at some more difficult types to define; namely, aromatic `sp2` carbon and aromatic `sp2` nitrogen. Starting with carbon:

	<Type name="ca" class="ca" element="C" def="[C;X3;r6]1[C&X3&r6,N&X2&r6][C&X3&r6,N&X2&r6][C&X3&r6,N&X2&r6][C&X3&r6,N&X2&r6][C&X3&r6,N&X2&r6]1" overrides="c2,cc_r5,cc_r6,ce,cz" desc="Sp2 C in pure aromatic systems" />

Here we use the `1[][][][][][]1` notation to describe a ring. Each atom in the ring must be a carbon with three connections or a nitrogen with two connections. GAFF defines "pure aromatic rings" as "heavy atoms in benzene and pyridine" (Case, 2006, pg. 256). 

The `na` and `nh` atomtypes are particularly hard to distinguish. The `na` type is a `sp2` nitrogen bonded to three other atoms. The `nh` type is defined as an amine nitrogen bound to one or more aromatic rings. The behavior of Antechamber and inspection of the `ATOMTYPE_GFF.dat` distributed with Amber suggests that the definition of `aromatic ring` adopted here is loose. The behavior suggests they mean an amine nitrogen bound to an `sp2` hybridized atom. We use that definition. 

First we define the `na` type by requiring it is bound to an atom that is also `sp2` hybrdized:

    <Type name="na" class="na" element="N" def="[N;X3][C&X3,C&X2,N&X2][C&X3,C&X2,N&X2]" overrides="n3" desc="Sp2 N with three connected atoms"/>

We then override the `na` type in many cases with the `nh` type:

	<Type name="nh" class="nh" element="N" def="[N;X3][C&X3,N&X2]" overrides="n3,na" desc="Sp2 N with three connected atoms"/>
	
However, this `nh` type definition is far too broad. Since the `nh` type is supposed to be connected to aromatic rings rather than within aromatic rings, we define two new artificial types, `na_r5` and `na_r6`. These are both of type `na` but allow us to create more specific definitions. Since they both have the class `na`, they will have the same bonded paramaters. We will add the types `na_r5` and `na_r6` to the nonbonded parameter section and ensure the same nonbonded parameters as `na`:

    <Type name="na_r5" class="na" element="N" def="[N;X3;r5]1[C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2]1" overrides="n3,nh,na,na_r5" desc="Sp2 N with three connected atoms"/>
    <Type name="na_r6" class="na" element="N" def="[N;X3;r6]1[C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2][C&X3,N&X3,N&X2,O&X2]1" overrides="n3,nh,na,na_r5" desc="Sp2 N with three connected atoms"/>
    
These two atomtypes both override `nh`. Note that `na_r6` overrides `na_r5`; without this override, if a nitrogen was found at the intersection of two fused rings where one was a six-membered ring and one was a five-membered ring, Foyer would find two matching atomtypes. Since both `na_r5` and `na_r6` represent the same atomtype in GAFF, the directionality of the override is irrelevant. Additionally, it might appear that one could use a definition with `[N;X3;!r5;!r6]` for the `nh` type. However, there are certain molecules for which the `nh` is connected to an aromatic ring and simultaneously in a non-aromatic ring. In those instances the prior definition would fail to properly identify the atom as `nh` even though it is connected to an aromatic ring.

That completes the 35 basic atomtypes in GAFF.

### Special atomtypes

#### Group I  

Group I types are for hydrogens in different chemical environments. 

    <Type name="h1" class="h1" element="H" def="H[C;X4]([N,O,F,Cl,Br,I,S])" overrides="hc" desc="H bonded to aliphatic carbon with 1 electrwd. group" />
    <Type name="h2" class="h2" element="H" def="H[C;X4]([N,O,F,Cl,Br,I,S])[N,O,F,Cl,Br,I,S]" overrides="h1" desc="H bonded to aliphatic carbon with 2 electrwd. group" />
    <Type name="h3" class="h3" element="H" def="H[C;X4]([N,O,F,Cl,Br,I,S])([N,O,F,Cl,Br,I,S])[N,O,F,Cl,Br,I,S]" overrides="h2" desc="H bonded to aliphatic carbon with 3 electrwd. group" />
    <Type name="h4" class="h4" element="H" def="H[C;!X4]([N,O,F,Cl,Br,I,S])" overrides="ha" desc="H bonded to non-sp3 carbon with 1 electrwd. group" />
    <Type name="h5" class="h5" element="H" def="H[C;!X4]([N,O,F,Cl,Br,I,S])([N,O,F,Cl,Br,I,S])" overrides="h4" desc="H bonded to non-sp3 carbon with 2 electrwd. group" />
	<Type name="s" class="s" element="S" def="[S;X1]" desc="S with one connected atom"/>

Notice the use of cascading overrides. Since `h1` overrides `hc`, the definition for `h2` overrides `hc` by overriding `h1`. In some cases it can be useful to include all the overrides explicitly. 

#### Group II 

Group II types contain carbons in three- and four-membered rings.

    <Type name="cx" class="cx" element="C" def="[C;X4;r3]" overrides="c3" desc="Sp3 carbons in triangle systems" />
    <Type name="cy" class="cy" element="C" def="[C;X4;r4]" overrides="c3" desc="Sp3 carbons in square systems" />
    <Type name="cu" class="cu" element="C" def="[C;X3;r3]" overrides="c2" desc="Sp2 carbons in triangle systems" />
    <Type name="cv" class="cv" element="C" def="[C;X3;r4]" overrides="c2" desc="Sp2 carbons in square systems" />

#### Group III 

Group III types are inner sp/sp<sup>2</sup> in conjugated systems. These atomtypes come in pairs: `cc (cd)`, `ce (cf)`, `cg (ch)`, `nc (nd)`, `ne (nf)`, `ne (nf)`, and `pe (pf)`. The first atom of each pair is group (i) whereas the latter is group (ii). Bonds between two atoms of the same group (e.g., `cc-cc`, `cc-ce`,`ce-ne`) are single bonds. Bonds between different groups (e.g., `cc-cd`, `cc-cf`,`ce-nf`) are double bonds. This allows GAFF to account for differences in the bonded parameters of, e.g., `ce-ce` single bond vs. `ce-cf` double bonds. These atomtypes are problematic for the SMARTS string graph matching approach used by Foyer since the surrounding chemical environments are identical for both members of a pair (e.g., `cc` and `cd`). 

For now, I have written SMARTS strings for the first member of each pair. Starting with the conjugated ring systems, we define `cc` as follows: 

	<Type name="cc_r5" class="cc" element="C" def="[C;X3;r5]1[C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2]1" overrides="c2,ce,cz" desc="Sp2 carbons in non-pure aromatic systems" />
    <Type name="cc_r6" class="cc" element="C" def="[C;X3;r6]1[C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2]1" overrides="c2,ce,cz,cc_r5" desc="Sp2 carbons in non-pure aromatic systems" /> 
    
Similar to the `na` type, I once again define two artificial types `cc_r5` and `cc_r6` which are both `cc`. Using two types allows a more precise definition of the other elements allowed to participate in the ring. 

The `nc` type is similarly defined:

	<Type name="nc_r5" class="nc" element="N" def="[N;X2;r5]1[C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2]1" overrides="ne,n2" desc="Sp2 N in non-pure aromatic systems" />
    <Type name="nc_r6" class="nc" element="N" def="[N;X2;r6]1[C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2][C&X3,N&X3,N&X2,O&X2,S&X2]1" overrides="ne,n2,nc_r5" desc="Sp2 N in non-pure aromatic systems" />
    
Moving to sp<sup>2</sup> conjugated chain systems, we define `ce`, `ne`, and `pe` as:

    <Type name="ce" class="ce" element="C" def="[C;X3]([C&X3,C&X2,N&X2])[C&X3,C&X2,N&X2,%pe,%px,%py,%s,%sx,%sy]" overrides="c2" desc="Inner Sp2 carbons in conjugated systems"/>
    <Type name="ne" class="ne" element="N" def="[N;X2]([C&X3,C&X2,N&X2,%pe,%px,%py,%s,%sx,%sy])[C&X2,C&X3,%pe,%px,%py,%sx,%sy][C&X3,C&X2,N&X2,%pe,%s,%sx,%sy,%o]" overrides="n2" desc="Inner Sp2 N in conjugated systems" />
    <Type name="pe" class="pe" element="P" def="[P;X2]([C&X3,C&X2,N&X2])[C&X3,C&X2,N&X2,%pe,%px,%py,%sx,%sy]" overrides="p2" desc="Inner Sp2 P in conjugated systems"/>

The sp-hybridized `cg` is defined as:

	<Type name="cg" class="cg" element="C" def="[C;X2]([C&X2,N&X1])[C&X3,C&X2,N&X2,%pe,%px,%py,%s,%sx,%sy]" overrides="c1" desc="Inner Sp carbons in conjugated systems"/>
	

#### Group IV

Group IV types once again come in pairs and are similarly problematic. `cp (cq)` are the head carbon atoms connecting the two rings in biphenyl systems. 
	
	<Type name="cp" class="cp" element="C" def="[C;X3;r6;R1]([C;X3;r6]1[C&X3,N&X2][C&X3,N&X2][C&X3,N&X2][C&X3,N&X2][C&X3,N&X2]1)([C&X3,N&X2&r6])[C&X3&r6,N&X2&r6]" overrides="ca,cc_r5,cc_r6" desc="Head Sp2 C that connect two rings in biphenyl sys." />

Here we require that the carbon atom is in a six membered ring, (`r6`), but only a single ring (`R1`). Otherwise the carbon atoms located at the vertex of two fused rings would be incorrectly identified as `cp`. 

Atomtypes `nb` and `pb` are nitrogen and phosphorus in pure aromatic systems. Following our previous definitions for atomtypes in pure aromatic systems:

    <Type name="nb" class="nb" element="N" def="[N;X2;r6]1[C&X3,N&X2,P&X2][C&X3,N&X2,P&X2][C&X3,N&X2,P&X2][C&X3,N&X2,P&X2][C&X3,N&X2,P&X2]1" overrides="nc_r6" desc="Sp2 N in pure aromatic systems" />
	<Type name="pb" class="pb" element="P" def="[P;X2;r6]1[C&X3,N&X2,P&X2][C&X3,N&X2,P&X2][C&X3,N&X2,P&X2][C&X3,N&X2,P&X2][C&X3,N&X2,P&X2]1" overrides="pe,pc_r6" desc="Sp2 P in pure aromatic systems" /> 

Next we consider `sx` and `sy`, which are described as sulfur in conjugated sulfoxide and sulfone systems. 

    <Type name="sx" class="sx" element="S" def="[S;X3]([O;X1])[C&X3,N&X2]" overrides="s4" desc="Special s4 in conjugated systems (sulfoxide)"/>
    <Type name="sy" class="sy" element="S" def="[S;X4]([O;X1])([O;X1])[C&X3,N&X2]" overrides="s6" desc="Special s6 in conjugated systems (sulfone)" />


Followed by `px` and `py`, which are phosphite and phosphate in conjugated systems.

    <Type name="px" class="px" element="P" def="[P;X3]([O;X1])[C&X3,N&X2]" overrides="p4" desc="Special p4 in conjugated systems (phosphite)"/>
    <Type name="py" class="py" element="P" def="[P;X4]([O;X1])([O;X2])([O;X2])[C&X3,N&X2,P&X4]" overrides="p5" desc="Special p5 in conjugated systems (phosphate)"/>

### Additional atomtypes

Other types that were not described in the original GAFF paper but included with antechamber include: 

    <Type name="cz" class="cz" element="C" def="[C;X3](N)(N)N" overrides="c2" desc="Sp2 carbon in guanidine group" />
    <Type name="hx" class="hx" element="H" def="H[C][N;%n4]" overrides="h1,h2,h3" desc="H bonded to C next to positively charged group" />
    <Type name="ni" class="ni" element="N" def="[N;X3;r3][C][O&X1,S&X1]" overrides="n,np" desc="like n in RG3" />
    <Type name="nj" class="nj" element="N" def="[N;X3;r4][C][O&X1,S&X1]" overrides="n,nq" desc="like n in RG4" />
    <Type name="nk" class="nk" element="N" def="[N;X4;r3]" overrides="n4" desc="like n4/nx/ny in RG3" />
    <Type name="nl" class="nl" element="N" def="[N;X4;r4]" overrides="n4" desc="like n4/nx/ny in RG4" />
    <Type name="nm" class="nm" element="N" def="[N;X3;r3][C&X3,N&X2]" overrides="nh" desc="like nh in RG3" />
    <Type name="nn" class="nn" element="N" def="[N;X3;r4][C&X3,N&X2]" overrides="nh" desc="like nh in RG4" />
    <Type name="np" class="np" element="N" def="[N;X3;r3]" overrides="n3" desc="like n3 in RG3" />
    <Type name="nq" class="nq" element="N" def="[N;X3;r4]" overrides="n3" desc="like n3 in RG4" />

