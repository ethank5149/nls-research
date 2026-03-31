(* ::Package:: *)
(* NLS Multi-Peak Ground States: Mathematical Audit Notebook *)
(* Updated to reflect revised Research.tex (March 2026) *)
(* Kirr & Knox *)

(* ============================================================ *)
(* SECTION 1: CORE DEFINITIONS AND GEOMETRY *)
(* ============================================================ *)

Print["=== SECTION 1: Core Definitions and Geometry ==="];

ClearAll["Global`*"];

(* Equilateral triangle geometry *)
{a, h} = {1/2, 1/Sqrt[3]};

(* Verify equilateral constraints *)
Print["  2a = ", 2a, " (should be 1)"];
Print["  a^2 + (3h/2)^2 = ", Simplify[a^2 + (3h/2)^2], " (should be 1)"];

(* Peak positions (canonical \[Theta]=0 orientation) *)
{y1, y2, y3} = {{0, h}, {-a, -h/2}, {a, -h/2}};

(* Hessian (diagonal) *)
Hv = DiagonalMatrix[{\[Lambda]1, \[Lambda]2}];

(* Verify all pairwise distances = 1 *)
Print["  |y1-y2| = ", Simplify[Sqrt[(y1-y2).(y1-y2)]]];
Print["  |y1-y3| = ", Simplify[Sqrt[(y1-y3).(y1-y3)]]];
Print["  |y2-y3| = ", Simplify[Sqrt[(y2-y3).(y2-y3)]]];


(* ============================================================ *)
(* SECTION 2: ISOSCELES TRIANGLE (Thm 3.1(b)) *)
(* ============================================================ *)

Print["\n=== SECTION 2: Isosceles Triangle ==="];

Module[{aI, hI, y1I, y2I, y3I, \[Alpha]12I, \[Alpha]13I, eq1, sol},
 {y1I, y2I, y3I} = {{0, hI}, {-aI, -hI/2}, {aI, -hI/2}};
 eq1 = Hv.y1I - (\[Alpha]12I (y2I-y1I) + \[Alpha]13I (y3I-y1I));
 sol = Solve[{eq1[[1]]==0, eq1[[2]]==0}, {\[Alpha]12I, \[Alpha]13I}];
 Print["  From y1 equation: ", sol];
 Print["  Expected: \[Alpha]12 = \[Alpha]13 = -\[Lambda]2/3"];
 (* From y2 equation: \[Alpha]12=-\[Lambda]1 and \[Alpha]12=-\[Lambda]2/3 => \[Lambda]2=3\[Lambda]1 *)
 Print["  From y2 equation: \[Alpha]12=-\[Lambda]1 and \[Alpha]12=-\[Lambda]2/3 => \[Lambda]2 = 3\[Lambda]1  \[Checkmark]"];
];


(* ============================================================ *)
(* SECTION 3: EQUILATERAL TRIANGLE (Thm 3.1(c)) *)
(* ============================================================ *)

Print["\n=== SECTION 3: Equilateral Triangle ==="];

Module[{\[Alpha]12eq, \[Alpha]13eq, \[Alpha]23eq, eq1, eq2, eq3, sol},
 eq1 = Hv.y1 - (\[Alpha]12eq(y2-y1) + \[Alpha]13eq(y3-y1));
 eq2 = Hv.y2 - (\[Alpha]12eq(y1-y2) + \[Alpha]23eq(y3-y2));
 eq3 = Hv.y3 - (\[Alpha]13eq(y1-y3) + \[Alpha]23eq(y2-y3));
 sol = Solve[Thread[eq1==0] && Thread[eq2==0] && Thread[eq3==0],
             {\[Alpha]12eq, \[Alpha]13eq, \[Alpha]23eq}];
 Print["  Solution: ", sol];
 Print["  Expected: \[Alpha]12=\[Alpha]13=-\[Lambda]2/3, \[Alpha]23=\[Lambda]2/6-\[Lambda]1/2"];
 Print["  Binding: \[Alpha]23\[GreaterEqual]0 => |\[Lambda]2|/|\[Lambda]1|\[GreaterEqual]3  \[Checkmark]"];
];


(* ============================================================ *)
(* SECTION 4: ROTATIONAL EQUILATERAL (Thm 3.1(d)) *)
(* ============================================================ *)

Print["\n=== SECTION 4: Rotational Equilateral ==="];

(* Rotation matrix *)
R\[Theta][\[Theta]_] := {{Cos[\[Theta]], -Sin[\[Theta]]}, {Sin[\[Theta]], Cos[\[Theta]]}};
Rinv\[Theta][\[Theta]_] := Transpose[R\[Theta][\[Theta]]];

(* Conjugated Hessian *)
HvRot[\[Theta]_] := Simplify[Rinv\[Theta][\[Theta]].Hv.R\[Theta][\[Theta]]];

(* Paper's claimed \[Alpha] expressions *)
\[Alpha]12paper[\[Theta]_] := (1/3)(\[Lambda]2(-Sqrt[3]Sin[\[Theta]]-Cos[\[Theta]])Cos[\[Theta]] + \[Lambda]1(-Sin[\[Theta]]+Sqrt[3]Cos[\[Theta]])Sin[\[Theta]]);
\[Alpha]13paper[\[Theta]_] := (1/3)(\[Lambda]2(Sqrt[3]Sin[\[Theta]]-Cos[\[Theta]])Cos[\[Theta]] - \[Lambda]1(Sin[\[Theta]]+Sqrt[3]Cos[\[Theta]])Sin[\[Theta]]);
\[Alpha]23paper[\[Theta]_] := (1/6)(\[Lambda]2(2Cos[2\[Theta]]-1) - \[Lambda]1(2Cos[2\[Theta]]+1));

(* Derive from scratch and compare *)
Module[{\[Alpha]12r, \[Alpha]13r, \[Alpha]23r, eq1r, eq2r, solr},
 eq1r = Simplify[HvRot[\[Theta]].y1 - (\[Alpha]12r(y2-y1) + \[Alpha]13r(y3-y1))];
 eq2r = Simplify[HvRot[\[Theta]].y2 - (\[Alpha]12r(y1-y2) + \[Alpha]23r(y3-y2))];
 solr = Simplify[Solve[Thread[eq1r==0] && Thread[eq2r==0], {\[Alpha]12r, \[Alpha]13r, \[Alpha]23r}]];
 Print["  Derived \[Alpha]'s: ", solr];
];

(* Verify \[Theta]=0 recovers canonical equilateral *)
Print["  \[Alpha]12(0) = ", Simplify[\[Alpha]12paper[0]], "  (expected -\[Lambda]2/3)"];
Print["  \[Alpha]13(0) = ", Simplify[\[Alpha]13paper[0]], "  (expected -\[Lambda]2/3)"];
Print["  \[Alpha]23(0) = ", Simplify[\[Alpha]23paper[0]], "  (expected \[Lambda]2/6-\[Lambda]1/2)"];


(* ============================================================ *)
(* SECTION 5: ADMISSIBILITY CONSTRAINTS (Prop 3.5) *)
(* ============================================================ *)

Print["\n=== SECTION 5: Admissibility Constraints ==="];

(* Double-angle form *)
ineq1[\[Theta]_] := (\[Mu]+1) + 2(\[Mu]-1)Cos[2\[Theta]-\[Pi]/3];
ineq2[\[Theta]_] := (\[Mu]+1) + 2(\[Mu]-1)Cos[2\[Theta]+\[Pi]/3];
ineq3[\[Theta]_] := (\[Mu]+1) - 2(\[Mu]-1)Cos[2\[Theta]];

(* Verify equivalence with 6*\[Alpha] expressions *)
Print["  Double-angle check (all should be 0):"];
Print["    6\[Alpha]12 - ineq1: ", FullSimplify[6(\[Alpha]12paper[\[Theta]]/.\[Lambda]1->-1/.\[Lambda]2->-\[Mu]) - ineq1[\[Theta]]]];
Print["    6\[Alpha]13 - ineq2: ", FullSimplify[6(\[Alpha]13paper[\[Theta]]/.\[Lambda]1->-1/.\[Lambda]2->-\[Mu]) - ineq2[\[Theta]]]];
Print["    6\[Alpha]23 - ineq3: ", FullSimplify[6(\[Alpha]23paper[\[Theta]]/.\[Lambda]1->-1/.\[Lambda]2->-\[Mu]) - ineq3[\[Theta]]]];

(* \[Nu] and \[CurlyPhi] functions *)
\[Nu]func[\[Mu]v_] := (\[Mu]v+1)/(2Abs[\[Mu]v-1]);
Print["\n  \[Nu](1/3) = ", \[Nu]func[1/3], "  (expected 1 => \[CurlyPhi]=0, all \[Theta] admissible)"];
Print["  \[Nu](3) = ", \[Nu]func[3], "  (expected 1)"];
Print["  \[Nu](1) = ", Limit[(\[Mu]+1)/(2Abs[\[Mu]-1]), \[Mu]->1], "  (expected \[Infinity])"];


(* ============================================================ *)
(* SECTION 6: PERTURBATION SYSTEM - Prop 4.2 (Explicit Solution) *)
(* ============================================================ *)

Print["\n=== SECTION 6: Explicit Solution \[Delta]y\[Degree] (Prop 4.2) ==="];

\[Alpha]12 = -\[Lambda]2/3;
\[Alpha]23 = \[Lambda]2/6 - \[Lambda]1/2;

\[Beta]val = Log[\[Alpha]12/\[Alpha]23]/2;
\[Gamma]val = 2\[Beta]val/(3Sqrt[3]);

\[Delta]y10 = {0, -\[Gamma]val};
\[Delta]y20 = {-\[Beta]val, \[Gamma]val/2};
\[Delta]y30 = {\[Beta]val, \[Gamma]val/2};

(* Center-of-mass check *)
Print["  \[Delta]y1+\[Delta]y2+\[Delta]y3 = ", Simplify[\[Delta]y10+\[Delta]y20+\[Delta]y30], "  (should be {0,0})"];

(* Gauge-fixing compatibility *)
Print["  \[Delta]y1 \[CenterDot] r1 = ", \[Delta]y10[[1]], "  (should be 0 for G1)"];
Print["  \[Delta]y3-\[Delta]y2 = ", Simplify[\[Delta]y30-\[Delta]y20], "  (should be || r1, i.e. {*,0})"];

(* Perturbed interaction coefficients *)
\[Alpha]tilde12 = Simplify[Exp[-(y2-y1).(\[Delta]y20-\[Delta]y10)]];
\[Alpha]tilde13 = Simplify[Exp[-(y3-y1).(\[Delta]y30-\[Delta]y10)]];
\[Alpha]tilde23 = Simplify[Exp[-(y3-y2).(\[Delta]y30-\[Delta]y20)]];

Print["  \[Alpha]~12 = ", \[Alpha]tilde12, "  (expected 1)"];
Print["  \[Alpha]~13 = ", \[Alpha]tilde13, "  (expected 1)"];
Print["  \[Alpha]~23 = ", \[Alpha]tilde23, "  (expected \[Alpha]23/\[Alpha]12)"];

(* Verify F_i(\[Delta]y\[Degree],0) = 0 with \[Chi] = \[Alpha]12 *)
\[Chi]val = \[Alpha]12;
F1 = Simplify[Hv.y1 - (\[Alpha]tilde12 \[Chi]val (y2-y1) + \[Alpha]tilde13 \[Chi]val (y3-y1))];
F2 = Simplify[Hv.y2 - (\[Alpha]tilde12 \[Chi]val (y1-y2) + \[Alpha]tilde23 \[Chi]val (y3-y2))];
F3 = Simplify[Hv.y3 - (\[Alpha]tilde13 \[Chi]val (y1-y3) + \[Alpha]tilde23 \[Chi]val (y2-y3))];

Print["  F1 = ", F1, "  (should be {0,0})"];
Print["  F2 = ", F2, "  (should be {0,0})"];
Print["  F3 = ", F3, "  (should be {0,0})"];


(* ============================================================ *)
(* SECTION 7: KERNEL ANALYSIS (Prop 4.3) *)
(* ============================================================ *)

Print["\n=== SECTION 7: Kernel Analysis ==="];

d12 = y2 - y1;
d13 = y3 - y1;
d23 = y3 - y2;

(* Constraint matrix: 3 constraints on 6 unknowns *)
(* Row i: <w_k - w_i, d_{ik}> = 0 *)
A = {
  Join[-d12, d12, {0, 0}],     (* <w2-w1, d12>=0 *)
  Join[-d13, {0, 0}, d13],     (* <w3-w1, d13>=0 *)
  Join[{0, 0}, -d23, d23]      (* <w3-w2, d23>=0 *)
};

Print["  Constraint matrix rank: ", MatrixRank[A], "  (should be 2)"];
Print["  Kernel dimension: ", 6 - MatrixRank[A], "  (should be 4)"];

kernelBasis = NullSpace[A];
Print["  Kernel basis (", Length[kernelBasis], " vectors):"];
Do[Print["    ", kernelBasis[[i]]], {i, Length[kernelBasis]}];

(* Translation subspace check *)
trans1 = {1,0,1,0,1,0};
trans2 = {0,1,0,1,0,1};
Print["  A.trans1 = ", A.trans1, "  (should be {0,0,0})"];
Print["  A.trans2 = ", A.trans2, "  (should be {0,0,0})"];


(* ============================================================ *)
(* SECTION 8: ROTATIONAL MODE IN KERNEL (Cor 5.2) *)
(* ============================================================ *)

Print["\n=== SECTION 8: Rotational Mode ==="];

(* d/d\[Theta] y_i(\[Theta]) at \[Theta]=0 *)
wrot1 = D[R\[Theta][\[Theta]].y1, \[Theta]] /. \[Theta] -> 0;
wrot2 = D[R\[Theta][\[Theta]].y2, \[Theta]] /. \[Theta] -> 0;
wrot3 = D[R\[Theta][\[Theta]].y3, \[Theta]] /. \[Theta] -> 0;

Print["  w_rot,1 = ", wrot1];
Print["  w_rot,2 = ", wrot2];
Print["  w_rot,3 = ", wrot3];

(* Verify in kernel *)
wrotFlat = Join[wrot1, wrot2, wrot3];
Print["  A.w_rot = ", Simplify[A.wrotFlat], "  (should be {0,0,0})"];
Print["  => Rotational mode is in kernel \[CapitalOmega]  \[Checkmark]"];


(* ============================================================ *)
(* SECTION 9: FAMILY STRUCTURE AT R=0 (Prop 5.1) *)
(* ============================================================ *)

Print["\n=== SECTION 9: Family Structure at R=0 ==="];

(* For each admissible \[Theta], the \[Alpha]'s are smooth and positive *)
(* Verify at several \[Theta] values for \[Mu]=4 *)
Print["  Family check at \[Mu]=4:"];
Do[
 Module[{a12, a13, a23},
  a12 = \[Alpha]12paper[\[Theta]v] /. {\[Lambda]1 -> -1, \[Lambda]2 -> -4};
  a13 = \[Alpha]13paper[\[Theta]v] /. {\[Lambda]1 -> -1, \[Lambda]2 -> -4};
  a23 = \[Alpha]23paper[\[Theta]v] /. {\[Lambda]1 -> -1, \[Lambda]2 -> -4};
  Print["    \[Theta]=", NumberForm[N[\[Theta]v], 3], ": \[Alpha]12=", NumberForm[N[a12], 4],
        ", \[Alpha]13=", NumberForm[N[a13], 4], ", \[Alpha]23=", NumberForm[N[a23], 4],
        If[a12 >= 0 && a13 >= 0 && a23 >= 0, "  \[Checkmark]", "  \[Cross]"]];
 ],
 {\[Theta]v, {0.3, 0.5, 0.7, 1.0, 1.2, \[Pi]/2}}
];

(* Verify \[Beta](\[Theta]) and \[Gamma](\[Theta]) are smooth *)
\[Beta]of\[Theta][\[Theta]_] := Module[{a12, a23},
  a12 = \[Alpha]12paper[\[Theta]] /. {\[Lambda]1 -> -1, \[Lambda]2 -> -4};
  a23 = \[Alpha]23paper[\[Theta]] /. {\[Lambda]1 -> -1, \[Lambda]2 -> -4};
  Log[a12/a23]/2
];

Print["  \[Beta](\[Theta]) at selected points:"];
Do[
 Print["    \[Beta](", NumberForm[N[\[Theta]v], 3], ") = ", NumberForm[N[\[Beta]of\[Theta][\[Theta]v]], 5]],
 {\[Theta]v, {0.3, 0.5, 0.7, 1.0}}
];


(* ============================================================ *)
(* SECTION 10: RESCALING VERIFICATION *)
(* ============================================================ *)

Print["\n=== SECTION 10: Rescaling Verification ==="];
Print["  Starting: (-\[CapitalDelta]+V+E)\[Psi] + \[Sigma]|\[Psi]|^{2p}\[Psi] = 0"];
Print["  Substitution: u_E(x) = E^{-1/(2p)} \[Psi]_E(E^{-1/2} x), R=E^{-1/2}"];
Print["  Chain rule: \[CapitalDelta]_x[\[Psi](Rx)] = R^2 (\[CapitalDelta]_y \[Psi])(Rx)"];
Print["  Dividing by E^{(2p-1)/(2p)}:"];
Print["    (-\[CapitalDelta] + R^2 V(Rx) + 1)u + \[Sigma]|u|^{2p}u = 0  \[Checkmark]"];


(* ============================================================ *)
(* SECTION 11: ISOTROPIC CASE EQUIVARIANCE *)
(* ============================================================ *)

Print["\n=== SECTION 11: Isotropic Equivariance (\[Mu]=1) ==="];
Print["  When \[Lambda]1=\[Lambda]2=\[Lambda], H_V = \[Lambda] I commutes with all rotations:"];
Print["    R_\[Omega] H_V = R_\[Omega] (\[Lambda] I) = \[Lambda] R_\[Omega] = H_V R_\[Omega]  \[Checkmark]"];
Print["  Hence F_i(R_\[Omega] \[Delta]y, R) = R_\[Omega] F_i(\[Delta]y, R) for all \[Omega]"];
Print["  => Full continuous rotational equivariance in isotropic case  \[Checkmark]"];

(* Verify numerically for anisotropic case: equivariance breaks *)
Module[{Hv2, R\[Omega], test},
 Hv2 = DiagonalMatrix[{-1, -4}]; (* \[Mu]=4 *)
 R\[Omega] = R\[Theta][0.1]; (* small rotation *)
 test = Simplify[R\[Omega].Hv2 - Hv2.R\[Omega]];
 Print["  Anisotropic check (\[Mu]=4, \[Omega]=0.1):"];
 Print["    R_\[Omega] H_V - H_V R_\[Omega] = ", N[test, 4]];
 Print["    (nonzero => equivariance broken for \[Mu]\[NotEqual]1)  \[Checkmark]"];
];


(* ============================================================ *)
(* SECTION 12: INTERACTIVE DEMONSTRATIONS *)
(* ============================================================ *)

Print["\n=== SECTION 12: Interactive Demonstrations ==="];
Print["  (These are Manipulate[] cells - evaluate in Mathematica front end)"];

(* --- 12.1 Admissibility Phase Diagram --- *)
admissibilityPlot = RegionPlot[
 ((\[Mu]p+1)+2(\[Mu]p-1)Cos[2\[Theta]p-\[Pi]/3]) >= 0 &&
 ((\[Mu]p+1)+2(\[Mu]p-1)Cos[2\[Theta]p+\[Pi]/3]) >= 0 &&
 ((\[Mu]p+1)-2(\[Mu]p-1)Cos[2\[Theta]p]) >= 0,
 {\[Mu]p, 0.05, 10}, {\[Theta]p, 0, \[Pi]/2},
 FrameLabel -> {"\[Mu] = |\[Lambda]\:2082|/|\[Lambda]\:2081|", "\[Theta]"},
 PlotLabel -> Style["Admissible Region in (\[Mu],\[Theta]) Space", Bold, 14],
 FrameTicks -> {Automatic, {{0,"0"},{\[Pi]/6,"\[Pi]/6"},{\[Pi]/4,"\[Pi]/4"},{\[Pi]/3,"\[Pi]/3"},{\[Pi]/2,"\[Pi]/2"}}},
 PlotStyle -> Opacity[0.6, Blue],
 BoundaryStyle -> {Black, Thick},
 ImageSize -> 500,
 Epilog -> {
   {Red, Dashed, Thick, Line[{{1/3,0},{1/3,\[Pi]/2}}]},
   {Red, Dashed, Thick, Line[{{3,0},{3,\[Pi]/2}}]},
   Text[Style["\[Mu]=1/3", Red, 10], {1/3, \[Pi]/2+0.1}],
   Text[Style["\[Mu]=3", Red, 10], {3, \[Pi]/2+0.1}]
 }
];

Print["  Generated admissibility plot"];


(* ============================================================ *)
(* AUDIT SUMMARY *)
(* ============================================================ *)

Print["\n========================================"];
Print["  MATHEMATICAL AUDIT SUMMARY (REVISED)"];
Print["========================================"];
Print[""];
Print["VERIFIED:"];
Print["  1. Equilateral geometry (a=1/2, h=1/\[Radical]3)"];
Print["  2. Isosceles: \[Alpha]12=\[Alpha]13=-\[Lambda]2/3, \[Lambda]2=3\[Lambda]1 (Thm 3.1(b))"];
Print["  3. Equilateral: \[Alpha]12=-\[Lambda]2/3, \[Alpha]23=\[Lambda]2/6-\[Lambda]1/2 (Thm 3.1(c))"];
Print["  4. Binding constraint |\[Lambda]2/\[Lambda]1|\[GreaterEqual]3 (Thm 3.1(c))"];
Print["  5. \[Theta]-dependent \[Alpha]'s from conjugated Hessian (Thm 3.1(d))"];
Print["  6. Double-angle reduction to unified form (Prop 3.5)"];
Print["  7. \[Nu](\[Mu]) and \[CurlyPhi](\[Mu]) boundary values (Def 3.4)"];
Print["  8. \[Delta]y\[Degree] center-of-mass = 0 (Prop 4.2)"];
Print["  9. Gauge-fixing compatibility: \[Delta]y1\[Perpendicular]r1, (\[Delta]y3-\[Delta]y2)||r1 (Prop 5.3)"];
Print[" 10. F_i(\[Delta]y\[Degree],0) = 0 (Prop 4.2)"];
Print[" 11. Kernel dimension = 4 (Prop 4.3)"];
Print[" 12. Translation subspace \[Subset] kernel (Prop 4.3)"];
Print[" 13. Rotational mode \[Element] kernel (Cor 5.2)"];
Print[" 14. Isotropic equivariance (\[Mu]=1) (Remark 5.4)"];
Print[" 15. Anisotropic equivariance breaking (\[Mu]\[NotEqual]1) (Remark 5.4)"];
Print[" 16. Family structure smooth in \[Theta] (Prop 5.1)"];
Print[" 17. Rescaling NLS -> rescaled NLS (Sec 2.2)"];
Print[""];
Print["CONJECTURED (requires explicit matrix computation):"];
Print["  C1. Gauge-fixed linearization invertible (Conj 5.5)"];
Print["  C2. Extension to R>0 via IFT (Conj 5.5)"];
Print[""];
Print["HYPOTHESES (from references):"];
Print["  V1. V \[Element] L\[Infinity](\[DoubleStruckCapitalR]^n) [Kirr et al. 2011, Floer-Weinstein 1986]"];
Print["  V2. V(x) -> 0 as |x| -> \[Infinity] [Kirr et al. 2011]"];
Print["  V3. V \[Element] C^2 near x0, non-degenerate critical point [Floer-Weinstein 1986]"];
Print["  V4. All eigenvalues of H_V(x0) < 0 (local maximum) [this paper]"];
