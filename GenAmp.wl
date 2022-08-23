(* ::Package:: *)

(* ::Title:: *)
(*Definitions*)


(* ::Subsubsection:: *)
(*Description:: *)


(* ::Section:: *)
(*Import*)


$FeynCalcStartupMessages = False;


$LoadAddOns={"FeynArts","FeynHelpers"};
Needs["FeynCalc`"];


SetOptions[FourVector,FeynCalcInternal->False];
$KeepLogDivergentScalelessIntegrals=True;


(* ::Section:: *)
(*temp functions*)


TermByTerm[transfo_,bar_:True][terms__]:=Module[{numterms,count=0,result,time},
numterms=Length[terms];
If[bar,Print["Current diagram: ",Dynamic[count],"      ","Number of diagrams: ",numterms];
Print[ProgressIndicator[Dynamic[count],{0,numterms}]];];
time=Timing[result=Plus@@((count++;transfo[#])&/@terms);]//First;
If[bar,Print["Time: ",time," s"];];
Return[result];
];


(* convert sum of terms to list and vice versa *)
PlusToList[e_Plus] := List @@ e;
PlusToList[e_] := {e};
ListToPlus[l_List] := Plus @@ l;


(* ::Section:: *)
(*Diagrams*)


FAPatch[PatchModelsOnly->True,FAModelsDirectory->StringJoin["./models/",$model]]


top[level_,i_,j_,opts___]:=CreateTopologies[level,i->j,Adjacencies-> {3,4},opts];
topCT[level_,i_,j_,opts___]:=CreateCTTopologies[level,i->j,Adjacencies-> {3,4},opts];
diag[level_,in__,out__]:=Module[{optsD,optsCD,diagrams,mod},
mod=StringJoin["./models/",$model,"/",$model];
If[(level==0)&&(Length[out]>0),optsD="ExcludeTopologies\[Rule]Tadpoles"; optsCD="ExcludeTopologies\[Rule]Tadpoles",optsD="ExcludeTopologies\[Rule]Reducible";optsCD="ExcludeTopologies\[Rule]Reducible" ];
If[Length[out]==0,optsD="TadpolesOnly"; optsCD="TadpoleCTsOnly"];
diagrams=(InsertFields[#,in->out,InsertionLevel->{Particles},Model->mod,GenericModel->mod]&)/@{top[level,Length[in],Length[out],ToExpression[optsD]],topCT[level,Length[in],Length[out],ToExpression[optsCD]]};
Paint[diagrams[[1]],ImageSize->{512,256}];
Paint[diagrams[[2]],ImageSize->{512,256}];
Return[diagrams];
];


(* ::Section:: *)
(*Amplitudes*)


$ampconfigFCFAConvert={ChangeDimension->False,Contract->False,DropSumOver->False,FCFADiracChainJoin->True,FeynAmpDenominatorCombine->True,
FinalSubstitutions->{},InitialSubstitutions->{},List->False,LoopMomenta->{l},LorentzIndexNames->{},Prefactor->1,SMP->False,SUNFIndexNames->{},
SUNIndexNames->{},TransversePolarizationVectors->{},UndoChiralSplittings->False};


$gaugerules={};


amp1L[inmom_,outmom_][diags_]:= Module[{amp,ampproc},
amp=CreateFeynAmp[#,GaugeRules->$gaugerules,Truncated-> True]&/@diags;
Options[FCFAConvert]={IncomingMomenta->inmom,OutgoingMomenta->outmom,$ampconfigFCFAConvert}//Flatten;
ampproc=(FCFAConvert[#])&/@amp 
];


temp1[sub__][exp1_]:=(((exp1/.sub)//Series[#,{Nd,0,1}]&)//Normal)/.{Nd:> 1};
temp2[exp2_]:=FCLoopIsolate[exp2//Contract//OneLoopSimplify[#,l]&,{l}];(*(most time consuming process)*)
temp3[exp3_]:=(ToPaVe[TID[exp3,l],l]//PaVeUVPart[#(*,Prefactor\[Rule]1/(2*Pi)^D*)]&) ;
temp4[exp4_]:=exp4//FCReplaceD[#,D->4-Epsilon]&//FCHideEpsilon[Normal[Series[#,{Epsilon,0,0}]]]& ;

amp1LSimplify[sub__,deltas__,bar_:False][amp__]:=Module[{Amp,AmpCT},
Amp=amp[[1]]//PlusToList//TermByTerm[(#//temp1[sub]//temp2//temp3)&,bar];
AmpCT=amp[[2]]//temp1[sub](*//Total*);
Return[(Amp+AmpCT)//temp4//Expand//SelectNotFree2[#,deltas]&];
];


(* ::Section:: *)
(*Storage*)


dvalues={};


dSave[exp_]:=Module[{expr},
expr=exp//First;
While[(Cases[dvalues,expr//Keys,\[Infinity]]//DeleteDuplicates//Length)>0,
dvalues=dvalues//Delete[Position[dvalues,expr//Keys,\[Infinity]][[1]][[1]]];
];
(*dvalues=dvalues//DeleteCases[#,_First,\[Infinity]]&;*)
AppendTo[dvalues,expr]
];


dSub[not_]:=Delete[dvalues,{Position[dvalues,not,\[Infinity]][[All,1]]}//Transpose]


(* ::Chapter:: *)
(*Package Setup*)


BeginPackage["GenAmp`"];


Begin["`Private`"];


 End[];


EndPackage[];
