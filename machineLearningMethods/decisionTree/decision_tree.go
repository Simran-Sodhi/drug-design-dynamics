package main

//this implements the decision tree machine learning model in go

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"time"
)

// findBestSplit finds the best feature and threshold to split the data
func findBestSplit(X [][]float64, y []float64, features []int) (int, float64, float64) {
	bestScore := math.Inf(1)
	bestFeature := 0
	bestThreshold := 0.0
	bestValue := 0.0

	// iterate through each feature
	for _, feature := range features {
		values := make([]float64, len(X))
		for i := range X {
			values[i] = X[i][feature]
		}
		sort.Float64s(values)

		// try different thresholds
		for i := 0; i < len(values)-1; i++ {
			threshold := (values[i] + values[i+1]) / 2
			leftSum, rightSum := 0.0, 0.0
			leftCount, rightCount := 0, 0
			leftSquaredSum, rightSquaredSum := 0.0, 0.0

			// split data based on threshold
			for j := range X {
				if X[j][feature] <= threshold {
					leftSum += y[j]
					leftSquaredSum += y[j] * y[j]
					leftCount++
				} else {
					rightSum += y[j]
					rightSquaredSum += y[j] * y[j]
					rightCount++
				}
			}

			// calculate score for this split
			if leftCount > 0 && rightCount > 0 {
				leftMean := leftSum / float64(leftCount)
				rightMean := rightSum / float64(rightCount)
				leftVar := leftSquaredSum/float64(leftCount) - leftMean*leftMean
				rightVar := rightSquaredSum/float64(rightCount) - rightMean*rightMean
				score := float64(leftCount)*leftVar + float64(rightCount)*rightVar

				// update best split if score is better
				if score < bestScore {
					bestScore = score
					bestFeature = feature
					bestThreshold = threshold
					bestValue = (leftSum + rightSum) / float64(leftCount+rightCount)
				}
			}
		}
	}
	return bestFeature, bestThreshold, bestValue
}

// train recursively builds the decision tree
func (dt *DecisionTree) train(X [][]float64, y []float64, depth int) {
	// stop conditions: no data or max depth reached
	if len(X) == 0 || depth >= 25 {
		dt.value = 0
		return
	}

	// calculate mean of y values
	sum := 0.0
	count := 0
	for _, val := range y {
		if !math.IsNaN(val) && !math.IsInf(val, 0) {
			sum += val
			count++
		}
	}
	if count > 0 {
		dt.value = sum / float64(count)
	}

	// stop if dataset is too small or max depth reached
	if len(X) < 5 || depth >= 25 {
		return
	}

	// select random subset of features
	nFeatures := len(X[0])
	featureCount := int(math.Sqrt(float64(nFeatures))) + 1
	features := rand.Perm(nFeatures)[:featureCount]

	// find best split
	feature, threshold, value := findBestSplit(X, y, features)
	dt.feature = feature
	dt.threshold = threshold
	dt.value = value

	// split data
	leftX := make([][]float64, 0)
	leftY := make([]float64, 0)
	rightX := make([][]float64, 0)
	rightY := make([]float64, 0)
	for i := range X {
		if X[i][feature] <= threshold {
			leftX = append(leftX, X[i])
			leftY = append(leftY, y[i])
		} else {
			rightX = append(rightX, X[i])
			rightY = append(rightY, y[i])
		}
	}

	// recursively train left and right subtrees
	if len(leftX) > 0 && len(rightX) > 0 {
		dt.left = &DecisionTree{}
		dt.right = &DecisionTree{}
		dt.left.train(leftX, leftY, depth+1)
		dt.right.train(rightX, rightY, depth+1)
	}
}

// predict traverses the decision tree to make a prediction
func (dt *DecisionTree) predict(x []float64) float64 {
	if dt.left == nil || dt.right == nil {
		return dt.value
	}
	if x[dt.feature] <= dt.threshold {
		return dt.left.predict(x)
	}
	return dt.right.predict(x)
}

// NewRandomForest creates a new random forest with specified number of trees
func NewRandomForest(nTrees int) *RandomForest {
	return &RandomForest{
		trees:  make([]*DecisionTree, nTrees),
		nTrees: nTrees,
	}
}

// train builds the random forest
func (rf *RandomForest) train(X [][]float64, y []float64) {
	fmt.Println("Training Random Forest...")
	nFolds := 5
	foldSize := len(X) / nFolds

	for i := 0; i < rf.nTrees; i++ {
		// cross-validation fold
		foldIndex := i % nFolds
		startIdx := foldIndex * foldSize
		endIdx := startIdx + foldSize

		validX := X[startIdx:endIdx]
		validY := y[startIdx:endIdx]

		// create training set
		trainX := make([][]float64, 0, len(X)-foldSize)
		trainY := make([]float64, 0, len(y)-foldSize)
		for j := 0; j < len(X); j++ {
			if j < startIdx || j >= endIdx {
				trainX = append(trainX, X[j])
				trainY = append(trainY, y[j])
			}
		}

		// bootstrap sampling
		n := len(trainX)
		sampleSize := int(float64(n) * 0.8)
		sampleX := make([][]float64, sampleSize)
		sampleY := make([]float64, sampleSize)
		for j := 0; j < sampleSize; j++ {
			idx := rand.Intn(n)
			sampleX[j] = trainX[idx]
			sampleY[j] = trainY[idx]
		}

		// train individual tree
		tree := &DecisionTree{}
		tree.train(sampleX, sampleY, 0)
		rf.trees[i] = tree

		// print validation error every 50 trees
		if (i+1)%50 == 0 {
			validError := 0.0
			validCount := 0
			for j := range validX {
				pred := tree.predict(validX[j])
				if !math.IsNaN(pred) && !math.IsInf(pred, 0) {
					diff := math.Abs(pred - validY[j])
					validError += diff
					validCount++
				}
			}
			if validCount > 0 {
				validError /= float64(validCount)
				fmt.Printf("Trained %d trees, Validation MAE: %.4f\n", i+1, validError)
			}
		}
	}
}

// predict makes a prediction using the random forest
func (rf *RandomForest) predict(x []float64) float64 {
	if len(rf.trees) == 0 {
		return 0
	}

	// get predictions from all trees
	predictions := make([]float64, len(rf.trees))
	for i, tree := range rf.trees {
		predictions[i] = tree.predict(x)
	}

	// remove outliers using IQR method
	sort.Float64s(predictions)
	q1 := predictions[len(predictions)/4]
	q3 := predictions[3*len(predictions)/4]
	iqr := q3 - q1
	lowerBound := q1 - 1.5*iqr
	upperBound := q3 + 1.5*iqr

	// calculate mean of non-outlier predictions
	sum := 0.0
	count := 0
	for _, pred := range predictions {
		if pred >= lowerBound && pred <= upperBound {
			sum += pred
			count++
		}
	}

	if count > 0 {
		return sum / float64(count)
	}
	return predictions[len(predictions)/2]
}

// ReadCSV reads protein data from a CSV file
func ReadCSV(filename string) ([]ProteinData, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.FieldsPerRecord = -1 // Allow variable column lengths
	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	var data []ProteinData
	for i, record := range records {
		if i == 0 {
			continue
		}

		// Parse CSV fields
		bindingAffinity, _ := strconv.ParseFloat(record[1], 64)
		bindingAffinitySD, _ := strconv.ParseFloat(record[2], 64)
		electrostatic, _ := strconv.ParseFloat(record[3], 64)
		electrostaticSD, _ := strconv.ParseFloat(record[4], 64)
		polarSolvation, _ := strconv.ParseFloat(record[5], 64)
		polarSolvationSD, _ := strconv.ParseFloat(record[6], 64)
		nonPolarSolvation, _ := strconv.ParseFloat(record[7], 64)
		nonPolarSolvationSD, _ := strconv.ParseFloat(record[8], 64)
		vdw, _ := strconv.ParseFloat(record[9], 64)

		// Create ProteinData struct
		data = append(data, ProteinData{
			PDB_ID:              record[0],
			BindingAffinity:     bindingAffinity,
			BindingAffinitySD:   bindingAffinitySD,
			Electrostatic:       electrostatic,
			ElectrostaticSD:     electrostaticSD,
			PolarSolvation:      polarSolvation,
			PolarSolvationSD:    polarSolvationSD,
			NonPolarSolvation:   nonPolarSolvation,
			NonPolarSolvationSD: nonPolarSolvationSD,
			VdW:                 vdw,
		})
	}

	return data, nil
}

func main() {
	rand.Seed(time.Now().UnixNano())
	fmt.Println("Starting protein binding prediction model...")

	// Read data from CSV
	data, err := ReadCSV("5kdata.csv")
	if err != nil {
		fmt.Printf("Error reading CSV: %v\n", err)
		return
	}

	// Prepare feature matrix X and target vector y
	X := make([][]float64, len(data))
	y := make([]float64, len(data))
	for i, d := range data {
		X[i] = []float64{
			d.Electrostatic,
			d.PolarSolvation,
			d.NonPolarSolvation,
			d.VdW,
		}
		y[i] = d.BindingAffinity
	}

	fmt.Println("Normalizing features...")

	// Split data into train and test sets
	trainSize := int(float64(len(data)) * 0.8)
	indices := rand.Perm(len(data))
	trainX := make([][]float64, trainSize)
	trainY := make([]float64, trainSize)
	testX := make([][]float64, len(data)-trainSize)
	testY := make([]float64, len(data)-trainSize)
	for i := 0; i < trainSize; i++ {
		trainX[i] = X[indices[i]]
		trainY[i] = y[indices[i]]
	}
	for i := trainSize; i < len(data); i++ {
		testX[i-trainSize] = X[indices[i]]
		testY[i-trainSize] = y[indices[i]]
	}

	// Train random forest
	rf := NewRandomForest(500) // Increased number of trees
	rf.train(trainX, trainY)

	// Make predictions on test set
	predictions := make([]float64, len(testX))
	for i, x := range testX {
		predictions[i] = rf.predict(x)
	}

	// Calculate and print performance metrics
	var mse, totalError float64
	errorBuckets := make(map[string]int)
	for i := range predictions {
		diff := math.Abs(predictions[i] - testY[i])
		mse += diff * diff
		totalError += diff
		fmt.Printf("PDB_ID: %s, Predicted: %.3f, Actual: %.3f, Diff: %.3f\n",
			data[indices[i+trainSize]].PDB_ID, predictions[i], testY[i], diff)
	}

	mse /= float64(len(predictions))
	rmse := math.Sqrt(mse)
	mae := totalError / float64(len(predictions))

	fmt.Printf("\nModel Performance:\n")
	fmt.Printf("Mean Squared Error: %.4f\n", mse)
	fmt.Printf("Root Mean Squared Error: %.4f\n", rmse)
	fmt.Printf("Mean Absolute Error: %.4f\n", mae)

	fmt.Printf("\nError Distribution:\n")
	for category, count := range errorBuckets {
		percentage := float64(count) * 100 / float64(len(predictions))
		fmt.Printf("%s: %d predictions (%.1f%%)\n", category, count, percentage)
	}
}
