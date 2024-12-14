# we used in built python libraries to compare different machine learning models
# we made use of sklearn library for this task
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import math
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
import xgboost as xgb
import lightgbm as lgb
from sklearn.ensemble import RandomForestRegressor  # make sure to import RandomForestRegressor

# load the data from a csv file into a dataframe
df = pd.read_csv('5kdata.csv')

# select features (independent variables) and the target variable (dependent variable)
features = ['electrostatic (kcal/mol)', 
            'polar_solvation (kcal/mol)', 
            'non_polar_solvation (kcal/mol)', 
            'vdW (kcal/mol)']
target = 'binding_affinity (kcal/mol)'

# extract features and target variable from the dataframe
X = df[features]
y = df[target]

# split the data into training (80%) and testing (20%) sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# scale the features using standard scaler (normalizes data to mean 0 and variance 1)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)  # fit and transform the training data
X_test_scaled = scaler.transform(X_test)  # transform the testing data (using training data parameters)

# define the models to train and evaluate
models = {
    "Random Forest": RandomForestRegressor(n_estimators=100, random_state=42),  # random forest with 100 estimators
    "Decision Tree": DecisionTreeRegressor(random_state=42),  # decision tree regressor
    "XGBoost": xgb.XGBRegressor(random_state=42),  # xgboost regressor
    "LightGBM": lgb.LGBMRegressor(random_state=42),  # lightgbm regressor
    "SVM": SVR()  # support vector regressor
}

# function to train a model and evaluate its performance
def evaluate_model(model, X_train, y_train, X_test, y_test):
    model.fit(X_train, y_train)  # train the model using the training data
    y_pred = model.predict(X_test)  # predict the target variable for the test set

    # calculate evaluation metrics
    mse = mean_squared_error(y_test, y_pred)  # mean squared error
    r2 = r2_score(y_test, y_pred)  # r-squared score
    mae = mean_absolute_error(y_test, y_pred)  # mean absolute error
    rmse = math.sqrt(mse)  # root mean squared error

    # mean absolute percentage error (MAPE) expressed as a percentage
    mape = np.mean(np.abs((y_test - y_pred) / y_test)) * 100

    
    return {
        'Model': model.__class__.__name__,
        'MSE': mse,
        'RMSE': rmse,
        'MAE': mae,
        'R²': r2,
        'MAPE (%)': mape
    }


results = []

# train and evaluate each model in the models dictionary
for name, model in models.items():
    model_results = evaluate_model(model, X_train_scaled, y_train, X_test_scaled, y_test)  # evaluate the model
    results.append(model_results)  


results_df = pd.DataFrame(results)


print("Model Evaluation Results:")
print(results_df)


for model_name in ['Random Forest', 'Decision Tree', 'XGBoost', 'LightGBM']:
    model = models[model_name]
    # check if the model has feature importance attribute
    if hasattr(model, 'feature_importances_'):
        feature_importance = pd.DataFrame({'feature': features, 'importance': model.feature_importances_})
        # sort the features by their importance in descending order
        feature_importance = feature_importance.sort_values('importance', ascending=False)
        print(f"\n{model_name} Feature Importance:")
        print(feature_importance)

# make a sample prediction using the random forest model (as an example)
sample = X_test_scaled[0].reshape(1, -1)  
rf_model = models["Random Forest"]  # retrieve the random forest model
prediction = rf_model.predict(sample)  
print(f"\nSample Prediction from Random Forest: {prediction[0]:.4f}")

# print actual vs predicted values for all models
for name, model in models.items():
    y_pred = model.predict(X_test_scaled)  
    
    actual_vs_predicted = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
    print(f"\nActual vs Predicted Values for {name}:")
    print(actual_vs_predicted)
