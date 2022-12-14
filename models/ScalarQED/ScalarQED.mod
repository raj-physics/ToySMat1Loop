(* Patched for use with FeynCalc *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"ScalarQED",
  Authors -> {"Raj Patil"},
  Emails -> {"raj.patil.physics@gmail.com"}};

FR$ClassesTranslation={};

FR$InteractionOrderPerturbativeExpansion={};

FR$GoldstoneList={};

(*     Declared indices    *)

(*     Declared particles    *)

M$ClassesDescription = {
V[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> "a",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> 0,
    Indices -> {} },

S[1] == {
    PropagatorLabel -> "\[Phi]",
    SelfConjugate -> False,
    QuantumNumbers -> {-Q},
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> mP,
    Indices -> {} }
}


(*        Definitions       *)


mP[ ___ ] := mP;




(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ S[1] , -S[1] ] == {{0, (-I)*mP^2*(-1 + ZmP^2*ZP)}, {0, (-I)*(-1 + ZP)}},

C[ S[1] , S[1] , -S[1] , -S[1] ] == {{(-I)*g4, (-I)*g4*(-1 + ZP^2*ZPPPP)}},

C[ S[1] , -S[1] , V[1] , V[1] ] == {{(2*I)*e^2, (2*I)*e^2*(-1 + ZA*ZAAPP*ZP)}},

C[ S[1] , -S[1] , V[1] ] == {{(-I)*e, (-I)*e*(-1 + Sqrt[ZA]*ZAPP*ZP)}, {I*e, I*e*(-1 + Sqrt[ZA]*ZAPP*ZP)}},

C[ V[1] , V[1] ] == {{0, (-I)*(-1 + ZA)}, {0, I*(-1 + ZA)}}

}

