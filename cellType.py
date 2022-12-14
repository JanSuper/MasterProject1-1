import numpy as np
import pandas as pd
import pickle


class CellType:
    def __init__(self, cluster_centers : np.ndarray):
        self.cluster_centers = cluster_centers
        self.cellType = ["Fibroblasts", "Epithelium", "Bcells", "Monocytes",
                        "Macrophages", "IL17", "T cells", "Others"]
        self.cellTypeDict = None

    def getCellType(self, cluster : int):
        if self.cellTypeDict is None:
            self.setCellType()
        return self.cellTypeDict(cluster)

    def setCellType(self):
        # Extract from the matlab code and convert to python
        self.cellTypeDict = {}
        save = np.array(self.cluster_centers)
        mask = np.ones(self.cluster_centers.shape[0])

        id = np.where(mask == 1, self.cluster_centers[:, 7], -100).argmax()
        self.cellTypeDict[id] = "Fibroblasts"
        mask[id] = 0

        id = np.where(mask == 1, self.cluster_centers[:, 13], -100).argmax()
        self.cellTypeDict[id] = "Epithelium"
        mask[id] = 0

        id = np.where(mask == 1, self.cluster_centers[:, 21], -100).argmax()
        self.cellTypeDict[id] = "Epithelium"
        mask[id] = 0

        id = np.where(mask == 1, self.cluster_centers[:, 14], -100).argmax()
        self.cellTypeDict[id] = "Bcells"
        mask[id] = 0

        id = np.where(mask == 1, self.cluster_centers[:, 19], -100).argmax()
        self.cellTypeDict[id] = "Monocytes"
        mask[id] = 0

        id = np.where(mask == 1, self.cluster_centers[:, 33], -100).argmax()
        self.cellTypeDict[id] = "Macrophages"
        mask[id] = 0
        
        tr = self.cluster_centers[:, 42] + self.cluster_centers[:, 7]
        id = np.where(mask == 1, tr, -100).argmax()
        self.cellTypeDict[id] = "IL17"
        mask[id] = 0

        id = np.where(mask == 1, self.cluster_centers[:, 46], -100).argmax()
        self.cellTypeDict[id] = "T cells"

        lis = list()
        for i in range(len(self.cluster_centers)):
            if i not in self.cellTypeDict:
                lis.append(i)

        for i in lis:
            self.cellTypeDict[i] = "Others"
        self.cluster_centers = save
                            
if __name__ == "__main__":
    with open("data/som.pkl", "rb") as f:
        som = pickle.load(f)
    cluster_centers = som.cluster_centers_
    cluster_centers = cluster_centers.reshape(16, 56)
    print(cluster_centers.shape)
    cellType = CellType(cluster_centers)
    cellType.setCellType()
    # dfn["pred"].apply(lambda x: cellType.cellType[x])
    for i in range(16):
        print(cellType.cellTypeDict[i], i)
