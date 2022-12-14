(* Patched for use with FeynCalc *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"TwoScalar4",
  Authors -> {"Raj Patil"},
  Emails -> {"raj.patil.physics@gmail.com"}};

FR$ClassesTranslation={Phi1 -> S[1], Phi2 -> S[2]};

FR$InteractionOrderPerturbativeExpansion={};

FR$GoldstoneList={};

(*     Declared indices    *)

IndexRange[ Index[Gen] ] = Range[ 2 ]

(*     Declared particles    *)

M$ClassesDescription = {
S[1] == {
    SelfConjugate -> True,
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> m1,
    Indices -> {},
    PropagatorLabel -> "Phi1" },

S[2] == {
    SelfConjugate -> True,
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> m2,
    Indices -> {},
    PropagatorLabel -> "Phi2" }
}


(*        Definitions       *)

FAGaugeXi[ S[1] ] = 1;

m1[ ___ ] := m1;
m2[ ___ ] := m2;




(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ S[1] , S[1] ] == {{0, (-I)*m1^2*(-1 + Zm1^2*ZP)}, {0, (-I)*(-1 + ZP)}},

C[ S[1] , S[1] , S[1] , S[1] ] == {{(-I)*g40, (-I)*(-g40 - 4*g31*th*Zg31*ZP^2 + g40*Zg40*ZP^2)}},

C[ S[1] , S[1] , S[1] , S[2] ] == {{(-I)*g31, (-I)*(-g31 - 3*g22*th*Zg22*ZP^2 + g31*Zg31*ZP^2 + g40*th*Zg40*ZP^2)}},

C[ S[2] , S[2] ] == {{0, (-I)*m2^2*(-1 + Zm2^2*ZP)}, {0, (-I)*(-1 + ZP)}},

C[ S[1] , S[1] , S[2] , S[2] ] == {{(-I)*g22, (-I)*(-4*g22 - 2*g13*th*Zg13*ZP^2 + g22*Zg22*ZP^2 + 2*g31*th*Zg31*ZP^2)}},

C[ S[1] , S[2] , S[2] , S[2] ] == {{(-I)*g13, (-I)*(-g13 - g04*th*Zg04*ZP^2 + g13*Zg13*ZP^2 + 3*g22*th*Zg22*ZP^2)}},

C[ S[2] , S[2] , S[2] , S[2] ] == {{(-I)*g04, (-I)*(-g04 + g04*Zg04*ZP^2 + 4*g13*th*Zg13*ZP^2)}},

C[ S[1] , S[2] ] == {{0, (-I)*th*(m1*Zm1 - m2*Zm2)*(m1*Zm1 + m2*Zm2)*ZP}, {0, 0}}

}

