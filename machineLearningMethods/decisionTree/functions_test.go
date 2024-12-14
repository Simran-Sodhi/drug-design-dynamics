package main

import (
	"math"
	"testing"
)

func TestFindBestSplit(t *testing.T) {
	// Test data for findBestSplit function
	X := [][]float64{
		{1.0, 2.0},
		{2.0, 3.0},
		{3.0, 4.0},
		{4.0, 5.0},
	}
	y := []float64{1.0, 2.0, 3.0, 4.0}
	features := []int{0, 1} // Only two features to split on

	feature, threshold, value := findBestSplit(X, y, features)

	// Checking expected values
	if feature != 0 {
		t.Errorf("Expected best feature 1, got %d", feature)
	}

	if threshold <= 0.0 {
		t.Errorf("Expected threshold > 0, got %.2f", threshold)
	}

	if math.IsNaN(value) {
		t.Errorf("Expected best value to be a valid number, got NaN")
	}
}

func TestDecisionTreeTrain(t *testing.T) {
	// Test data for DecisionTree training
	X := [][]float64{
		{1.0, 2.0},
		{2.0, 3.0},
		{3.0, 4.0},
		{4.0, 5.0},
	}
	y := []float64{1.0, 2.0, 3.0, 4.0}
	tree := &DecisionTree{}
	tree.train(X, y, 0) // Train with a maximum depth of 0 (root node only)

	// Checking expected results
	if tree.value == 0 {
		t.Errorf("Expected tree value to be set during training")
	}

}

func TestDecisionTreePredict(t *testing.T) {
	// Predefined decision tree for testing prediction
	tree := &DecisionTree{
		value:     3.0,
		feature:   0,
		threshold: 2.5,
		left: &DecisionTree{
			value: 1.0,
		},
		right: &DecisionTree{
			value: 5.0,
		},
	}

	// Test cases with known expected output
	tests := []struct {
		input  []float64
		output float64
	}{
		{[]float64{1.0, 1.0}, 1.0}, // Less than 2.5, so should go left
		{[]float64{3.0, 1.0}, 5.0}, // Greater than 2.5, go right
		{[]float64{2.0, 2.0}, 1.0}, // Less than 2.5, go left
	}

	// Run tests
	for _, tt := range tests {
		t.Run("Test prediction", func(t *testing.T) {
			got := tree.predict(tt.input)
			if got != tt.output {
				t.Errorf("Expected %f, got %f", tt.output, got)
			}
		})
	}
}

func TestRandomForestTrain(t *testing.T) {
	// Test data for RandomForest training
	X := [][]float64{
		{1.0, 2.0},
		{2.0, 3.0},
		{3.0, 4.0},
		{4.0, 5.0},
	}
	y := []float64{1.0, 2.0, 3.0, 4.0}
	rf := NewRandomForest(2) // Using 2 trees instead of 3 for this test

	rf.train(X, y)

	// Verifying the number of trees in the RandomForest
	if len(rf.trees) != 2 {
		t.Errorf("Expected 2 trees, got %d", len(rf.trees))
	}

	// Testing the prediction function
	testX := []float64{3.0, 3.0}
	prediction := rf.predict(testX)
	if math.IsNaN(prediction) {
		t.Errorf("Expected a valid prediction, got NaN")
	}
}
