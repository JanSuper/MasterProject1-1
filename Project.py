import pandas as pd
import pickle
import numpy as np
from sklearn_som.som import SOM
from sklearn.preprocessing import StandardScaler
from cellType import CellType
import seaborn as sns
import matplotlib.pyplot as plt

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
        self.trainbool = False

    def transform(self):
        numeric_cols = self.data.select_dtypes(include=np.number).columns
        scaled = StandardScaler().fit_transform(self.data[numeric_cols])
        self.data_scaled = pd.concat([self.data.drop(numeric_cols, axis=1), 
                                    pd.DataFrame(scaled, columns=numeric_cols)],
                                    axis=1)

    def train(self):
        self.trainbool = True
        if self.data_scaled is None:
            self.transform()
        groupe = self.data_scaled.groupby("ImageNumber")
        for name, group in groupe:
            features = group._get_numeric_data().values
            # print(features.shape)
            self.som.fit(features, self.n_iter)


    def plot(self):
        if not self.trainbool:
            raise Exception("Model not trained")
        _, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,7))
        sns.scatterplot(x="cx", y="cy", hue='cellType', data=self.Data, sizes=(15, 15), ax=ax)
        ax.set_title("Cell Type")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Cell Type")
        plt.show()

        pass

    def assigncluster(self):
        if not self.trainbool:
            raise Exception("Model not trained")
        cluster = self.som.cluster_centers_.reshape(self.n_rows * self.n_cols
                                                    , -1)
        ct = CellType(cluster)
        ct.setCellType()
        self.Data = self.data.copy()
        self.Data["cluster"] = self.som.predict(self.data_scaled._get_numeric_data().values)
        self.Data["cellType"] = self.Data["cluster"].apply(lambda x: ct.cellTypeDict[x])

    def save(self, path : str):
        self.trainbool = True
        with open(path, "wb") as f:
            pickle.dump(self.som, f)

    def load(self, path : str):
        self.trainbool = True
        with open(path, "rb") as f:
            self.som = pickle.load(f)


if __name__ == "__main__":
    data = pd.read_csv("data/full_data_clean.csv")
    som = Model(data, 4, 3)
    som.train()
    som.assigncluster()
    som.plot()