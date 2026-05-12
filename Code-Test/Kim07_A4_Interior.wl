(*
   Kim07 interior A4 optimization via KKT linear system.
   Reproduces Tables 1 and 2 of Kim (2007), AIAA J. 45(1).

   Run with:  wolframscript -file Kim07_A4_Interior.wl

   Scheme (1st derivative, interior):
     beta f'_{i-2} + alpha f'_{i-1} + f'_i + alpha f'_{i+1} + beta f'_{i+2}
         = (1/h) [a1(f_{i+1}-f_{i-1}) + a2(f_{i+2}-f_{i-2}) + a3(f_{i+3}-f_{i-3})]

   4th-order Taylor constraints (matching kappa^1 and kappa^3 in N = kappa*D):
     -2*alpha - 2*beta + 2*a1 + 4*a2 + 6*a3 = 1   (Kim07 Eq 2)
     -3*alpha - 12*beta + a1 + 8*a2 + 27*a3 = 0   (Kim07 Eq 3)

   Spectral error functional (Kim07 Eq 7):
     E = Integrate[((1+delta)*kappa*D(kappa) - N(kappa))^2 * (kappa/r)^n, {kappa, 0, r}]
     D(kappa) = 1 + 2*alpha*Cos[kappa] + 2*beta*Cos[2*kappa]
     N(kappa) = 2*(a1*Sin[kappa] + a2*Sin[2*kappa] + a3*Sin[3*kappa])

   E is quadratic in x = {alpha, beta, a1, a2, a3} and constraints are linear,
   so the constrained minimum solves the 7x7 KKT system:
     {{M, C^T}, {C, 0}} . {x, lambda} = {-f, q}
   where M[i,j] = Integral[Phi_i * Phi_j * w], f[i] = Integral[Phi0 * Phi_i * w],
   w = (kappa/r)^n.

   Basis functions (residual = (1+delta)*kappa*D - N):
     Phi0   = (1+delta)*kappa
     Phi[1] = 2*(1+delta)*kappa*Cos[kappa]   (coeff of alpha)
     Phi[2] = 2*(1+delta)*kappa*Cos[2*kappa] (coeff of beta)
     Phi[3] = -2*Sin[kappa]                   (coeff of a1)
     Phi[4] = -2*Sin[2*kappa]                 (coeff of a2)
     Phi[5] = -2*Sin[3*kappa]                 (coeff of a3)
*)

(* ---- Basis functions ---- *)
Phi0[kap_, delt_] := (1 + delt) * kap;

PhiVec[kap_, delt_] := {
   2 (1 + delt) kap Cos[kap],
   2 (1 + delt) kap Cos[2 kap],
  -2 Sin[kap],
  -2 Sin[2 kap],
  -2 Sin[3 kap]
};

(* ---- KKT solver ---- *)
SolveKim07A4[rVal_?NumericQ, deltaVal_?NumericQ, nVal_?NumericQ] :=
  Module[{wt, gramMat, rhsVec, constrMat, constrRhs, kktMat, kktRhs, sol},

    wt[kap_] := (kap / rVal)^nVal;

    (* Gram matrix  M[i,j] = Integral[Phi_i * Phi_j * w] *)
    gramMat = Table[
      NIntegrate[
        PhiVec[kap, deltaVal][[i]] * PhiVec[kap, deltaVal][[j]] * wt[kap],
        {kap, 0, rVal},
        MaxRecursion -> 30
      ],
      {i, 5}, {j, 5}
    ];

    (* RHS vector  f[i] = Integral[Phi0 * Phi_i * w] *)
    rhsVec = Table[
      NIntegrate[
        Phi0[kap, deltaVal] * PhiVec[kap, deltaVal][[i]] * wt[kap],
        {kap, 0, rVal},
        MaxRecursion -> 30
      ],
      {i, 5}
    ];

    (* Taylor constraint matrix  C x = q,  x = {alpha, beta, a1, a2, a3} *)
    constrMat = {{-2., -2.,  2.,  4.,  6.},
                 {-3., -12., 1.,  8., 27.}};
    constrRhs = {1., 0.};

    (* KKT system  {{M, C^T}, {C, 0}} . {x, lambda} = {-f, q} *)
    kktMat = Join[
      Join[gramMat,   Transpose[constrMat], 2],
      Join[constrMat, ConstantArray[0., {2, 2}], 2]
    ];
    kktRhs = Join[-rhsVec, constrRhs];

    sol = LinearSolve[kktMat, kktRhs];

    <|"r"     -> rVal,   "delta" -> deltaVal, "n"  -> nVal,
      "alpha" -> sol[[1]], "beta"  -> sol[[2]],
      "a1"    -> sol[[3]], "a2"    -> sol[[4]], "a3" -> sol[[5]]|>
  ];

(* ---- Resolution limit kappa_c^{0.1%} / Pi ---- *)
ResolutionLimit[coeff_, tol_: 0.001] :=
  Module[{alp, bet, a1v, a2v, a3v, kapGrid, denom, numer, kapBar, relErr, idx},
    alp  = coeff["alpha"]; bet  = coeff["beta"];
    a1v  = coeff["a1"];    a2v  = coeff["a2"]; a3v = coeff["a3"];
    kapGrid = Range[1*^-8, Pi, Pi/10000];
    denom  = 1 + 2 alp Cos[#] + 2 bet Cos[2 #] & /@ kapGrid;
    numer  = 2 a1v Sin[#] + 2 a2v Sin[2 #] + 2 a3v Sin[3 #] & /@ kapGrid;
    kapBar = numer / denom;
    relErr = Abs[kapBar / kapGrid - 1];
    idx    = SelectFirst[Range[Length[kapGrid]], relErr[[#]] > tol &, Length[kapGrid] + 1];
    If[idx <= Length[kapGrid], kapGrid[[idx]] / Pi, 1.0]
  ];

(* ---- Print helper ---- *)
PrintResult[res_, label_String : ""] := (
  If[label =!= "", Print["\n--- ", label, " ---"]];
  Print["  alpha = ", NumberForm[res["alpha"], 16]];
  Print["  beta  = ", NumberForm[res["beta"],  16]];
  Print["  a1    = ", NumberForm[res["a1"],    16]];
  Print["  a2    = ", NumberForm[res["a2"],    16]];
  Print["  a3    = ", NumberForm[res["a3"],    16]];
  Print["  kappa_c^0.1% = ", NumberForm[ResolutionLimit[res], 4], " pi"]
);

(* ---- Reference values from Kim07 Table 1 ---- *)
table1Ref = <|
  "alpha" -> 0.5862704032801503,
  "beta"  -> 0.0954953355017055,
  "a1"    -> 0.6431406736919156,
  "a2"    -> 0.2586011023495066,
  "a3"    -> 0.007140953479797375
|>;

(* ---- Cases from Kim07 Table 2 ---- *)
table2Cases = {
  {2.661, -0.000210, 15, "r=2.661  delta=-0.000210  n=15"},
  {2.671, -0.000246, 11, "r=2.671  delta=-0.000246  n=11"},
  {2.669, -0.000260, 10, "r=2.669  delta=-0.000260  n=10"},
  {2.672, -0.000233, 10, "r=2.672  delta=-0.000233  n=10"}
};

(* ==================================================================== *)
Print["===================================================="];
Print["Kim07 Interior A4 Optimization  (KKT linear system)"];
Print["===================================================="];

(* Table 1 *)
Print["\n=== TABLE 1: r=2.672, delta=-0.000233, n=10 ==="];
res1 = SolveKim07A4[2.672, -0.000233, 10];
PrintResult[res1, "computed"];

Print["\nKim07 Table 1 reference:"];
Scan[(Print["  ", #, " = ", NumberForm[table1Ref[#], 16]]) &,
     {"alpha", "beta", "a1", "a2", "a3"}];

Print["\nErrors vs Kim07 Table 1:"];
Scan[(Print["  |", #, " - ref| = ",
            ScientificForm[Abs[res1[#] - table1Ref[#]], 2]]) &,
     {"alpha", "beta", "a1", "a2", "a3"}];

(* Table 2 *)
Print["\n=== TABLE 2 ==="];
Scan[
  Function[tc,
    Module[{rv, dv, nv, lbl, res},
      {rv, dv, nv, lbl} = tc;
      res = SolveKim07A4[rv, dv, nv];
      PrintResult[res, lbl]
    ]
  ],
  table2Cases
];
