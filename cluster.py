import pandas as pd
import numpy as np
import ast
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt



df = pd.read_csv('protein_data.csv')



def get_features(coords_str):
    coords = np.array(ast.literal_eval(coords_str))
    
    # Feature 1: Compactness (Variance of coordinates), dense folds have low variance/spread
    compactness = np.var(coords[:, 0]) + np.var(coords[:, 1])
    
    # Feature 2: Bounding Box Area
    width = np.max(coords[:, 0]) - np.min(coords[:, 0])
    height = np.max(coords[:, 1]) - np.min(coords[:, 1])
    area = width * height
    
    return [compactness, area]


features = np.array([get_features(c) for c in df['coords']])
energies = df['energy'].values


# Run K-Means
# look for 3 clusters sticks/ L-shapes/ Dense cores
kmeans = KMeans(n_clusters=3, random_state=42)
df['cluster'] = kmeans.fit_predict(features)


# Visualization
plt.figure(figsize=(10, 6))
scatter = plt.scatter(features[:, 0], energies, c=df['cluster'], cmap='viridis')
plt.colorbar(scatter, label='Cluster ID')
plt.xlabel('Spread/Compactness (Lower is more dense)')
plt.ylabel('Energy (Lower is more stable)')
plt.title('Cluster Analysis of Protein Folds')
plt.show()

print("Clustering Complete")