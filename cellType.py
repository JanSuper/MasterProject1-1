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
        pairs = []
        for i in lis:
            # self.cellTypeDict[i] = "T cells"
            min = -1
            for y in self.cellTypeDict:                
                dist = np.linalg.norm(i-y)
                if min == -1 or dist < min:
                    min = dist
                    id = y
            pairs.append((i, id))
        for i in pairs:
            self.cellTypeDict[i[0]] = self.cellTypeDict[i[1]]
                    
                            
if __name__ == "__main__":
    n_cluster = 9
    with open("data/som.pkl", "rb") as f:
        som = pickle.load(f)
    cluster_centers = som.cluster_centers_
    cluster_centers = cluster_centers.reshape(n_cluster, 56)
    print(cluster_centers.shape)
    cellType = CellType(cluster_centers)
    cellType.setCellType()
    # dfn["pred"].apply(lambda x: cellType.cellType[x])
    for i in range(n_cluster):
        print(cellType.cellTypeDict[i], i)
