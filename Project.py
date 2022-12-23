import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt
from sklearn_som.som import SOM
from sklearn.preprocessing import StandardScaler

class Model:
    def __init__(self, data : pd.DataFrame, n_rows : int, n_cols : int, 
                n_iter : int = 1000, random_seed : int = 10):
        self.data = data
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.som = SOM(n_rows, n_cols, data.shape[1] - 2,
                        random_state=random_seed)
        self.n_iter = n_iter
        self.data_scaled = None

    def transform(self):
        numeric_cols = self.data.select_dtypes(include=np.number).columns
        scaled = StandardScaler().fit_transform(self.data[numeric_cols])
        self.data_scaled = pd.concat([self.data.drop(numeric_cols, axis=1), 
                                    pd.DataFrame(scaled, columns=numeric_cols)],
                                    axis=1)

    def train(self):
        if self.data_scaled is None:
            self.transform()
        groupe = self.data_scaled.groupby("ImageNumber")
        for name, group in groupe:
            features = group._get_numeric_data().values
            print(features.shape)
            self.som.fit(features, self.n_iter*len(features))

if __name__ == "__main__":
    data = pd.read_csv("data/full_data_clean.csv")
    som = Model(data, 4, 3)
    # som.transform()
    som.train()
    print(som.data_scaled)

    
