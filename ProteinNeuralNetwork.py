import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler


df = pd.read_csv('protein_data.csv')


# Used to get the density
def get_average_pairwise_distance(coords):
    distances = []
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            dist = np.linalg.norm(np.array(coords[i]) - np.array(coords[j]))
            distances.append(dist)
    return np.mean(distances)



def extract_final_features(row):
    coords = np.array(ast.literal_eval(row['coords']))
    seq = row['sequence']
    
    # Geometric Features
    compactness = np.var(coords[:, 0]) + np.var(coords[:, 1])
    width = np.max(coords[:, 0]) - np.min(coords[:, 0])
    height = np.max(coords[:, 1]) - np.min(coords[:, 1])
    
    # Physical Features
    h_ratio = seq.count('H') / len(seq)
    avg_dist = get_average_pairwise_distance(coords)
    
    return [compactness, width, height, h_ratio, avg_dist]


X = np.array([extract_final_features(row) for _, row in df.iterrows()])
y = df['energy'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Model Parameters for Stability
nn = MLPRegressor(
    hidden_layer_sizes=(64,), # just 1 layer
    activation='relu',        
    solver='adam',            
    alpha=0.0001, 
    learning_rate_init=0.01,
    max_iter=2000, 
    random_state=42
)

print("Training the optimized Neural Network...")
nn.fit(X_train_scaled, y_train)

accuracy = nn.score(X_test_scaled, y_test)
print(f"\nFinal Neural Network R^2 Accuracy: {accuracy:.4f}")

# Test on a sample
sample_pred = nn.predict(X_test_scaled[:1])
print(f"Actual Energy: {y_test[0]} | AI Predicted Energy: {sample_pred[0]:.2f}")



y_pred = nn.predict(X_test_scaled)
residuals = y_test - y_pred

plt.figure(figsize=(8, 5))
plt.scatter(y_test, residuals, color='crimson')
plt.axhline(y=0, color='black', linestyle='--')
plt.xlabel('Actual Energy')
plt.ylabel('Prediction Error (Residuals)')
plt.title('Residual Analysis')
plt.show()