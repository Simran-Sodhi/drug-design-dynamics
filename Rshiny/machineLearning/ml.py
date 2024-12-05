import subprocess
import sys

def ensure_libraries_installed(packages):
    """
    Ensure the required libraries are installed. If not, install them.

    """
    for package, module_name in packages.items():
        try:
            __import__(module_name)
        except ImportError:
            print(f"Module '{module_name}' not found. Installing '{package}'...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Define the required packages and their corresponding module names
required_packages = {
    "pandas": "pandas",
    "numpy": "numpy",
    "scikit-learn": "sklearn",
    "xgboost": "xgboost",
    "lightgbm": "lightgbm"
}

# Ensure all required libraries are installed
ensure_libraries_installed(required_packages)


import pandas as pd
import numpy as np
import argparse
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import math
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
import xgboost as xgb
import lightgbm as lgb
from sklearn.ensemble import RandomForestRegressor


# Model definitions
available_models = {
    "RandomForest": RandomForestRegressor(n_estimators=100, random_state=42),
    "DecisionTree": DecisionTreeRegressor(random_state=42),
    "XGBoost": xgb.XGBRegressor(random_state=42),
    "LightGBM": lgb.LGBMRegressor(random_state=42),
    "SVM": SVR()
}


# Set your dependent variable
target = 'binding_affinity (kcal/mol)' 

# Set variables that sotore the results
results = []
feature_importance = {}


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Evaluate models and feature importance.")
parser.add_argument("--data", required=True, help="Path to the dataset (CSV).")
parser.add_argument("--features", nargs="+", required=True, help="Selected features.")
parser.add_argument("--models", nargs="+", required=True, help="Selected models.")
parser.add_argument("--output_Dir", nargs="+", required=True, help="Selected models.")
args = parser.parse_args()


# Get input data/parameters
data = args.data
features = args.features
models = args.models
output_dir = args.output_Dir[0] if isinstance(args.output_Dir, list) else args.output_Dir

# check the output directory
print(f"Output Directory: {output_dir}")



# Load the data
df = pd.read_csv(data)

if target not in df.columns:
    raise ValueError(f"Target column '{target}' not found in the dataset.")

# Get X and y
X = df[features]
y = df[target]


# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Scale the features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


# Filter models based on input
selected_models = {name: available_models[name] for name in models}

# Define function to train and evaluate models
def evaluate_model(model, X_train, y_train, X_test, y_test):
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    rmse = math.sqrt(mse)

    mape = np.mean(np.abs((y_test - y_pred) / y_test)) 
    
    return {
        'Model': model.__class__.__name__,
        'MSE': mse,
        'RMSE': rmse,
        'MAE': mae,
        'R2': r2,
        'MAPE': mape
    }

# Evaluate all models
for name, model in selected_models.items():
    model_results = evaluate_model(model, X_train_scaled, y_train, X_test_scaled, y_test)
    model_results['Model'] = name
    results.append(model_results)

# Output evalution results as .csv file for Rshiny
results_df = pd.DataFrame(results)
results_df.to_csv(f"{output_dir}/evaluation_results.csv", index=False)

# Check evalution results
# print("Model Evaluation Results:")
# print(results_df)


# Calculate feature importance
for model_name in models:
    model = selected_models[model_name]
    if hasattr(model, 'feature_importances_'):  # Check if the model has 'feature_importances_' attribute
        fi = pd.DataFrame({
            'feature': features,
            'importance': model.feature_importances_,
            'model': model_name
        })
        fi = fi.sort_values('importance', ascending=False)
        feature_importance[model_name] = fi

# Output feature importance tables as .csv file for Rshiny
if feature_importance:
    for model_name, fi_df in feature_importance.items():
        output_file = f"{output_dir}/feature_importance_{model_name}.csv"
        fi_df.to_csv(output_file, index=False)
        print(f"Feature importance saved for {model_name} to {output_file}")


print("All Results saved to:", f"{output_dir}")