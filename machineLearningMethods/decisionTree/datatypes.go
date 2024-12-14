package main

type ProteinData struct {
	PDB_ID              string
	BindingAffinity     float64
	BindingAffinitySD   float64
	Electrostatic       float64
	ElectrostaticSD     float64
	PolarSolvation      float64
	PolarSolvationSD    float64
	NonPolarSolvation   float64
	NonPolarSolvationSD float64
	VdW                 float64
}
type DecisionTree struct {
	value     float64
	feature   int
	threshold float64
	left      *DecisionTree
	right     *DecisionTree
}

type RandomForest struct {
	trees  []*DecisionTree
	nTrees int
}
