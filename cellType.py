import numpy as np
import pandas as pd
import pickle


class CellType:
    def __init__(self, cluster_centers : np.ndarray):
        self.cluster_centers = cluster_centers
        self.cellType = ["Fibroblasts", "Epithelium", "Bcells", "Monocytes",
                        "Macrophages", "IL17", "T cells", "Others"]
        self.cellTypeDict = {}

    def getCellType(self, cluster : int):
        return self.cellTypeDict(cluster)

    def setCellType(self):
        # Extract from the matlab code and convert to python
        self.cellTypeDict[self.cluster_centers[:, 7].argmax()] = "Fibroblasts"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 7].argmax(), axis = 0)
        self.cellTypeDict[self.cluster_centers[:, 13].argmax()] = "Epithelium"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 13].argmax(), axis = 0)
        self.cellTypeDict[self.cluster_centers[:, 21].argmax()] = "Epithelium"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 21].argmax(), axis = 0)
        self.cellTypeDict[self.cluster_centers[:, 14].argmax()] = "Bcells"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 14].argmax(), axis = 0)
        self.cellTypeDict[self.cluster_centers[:, 19].argmax()] = "Monocytes"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 19].argmax(), axis = 0)
        self.cellTypeDict[self.cluster_centers[:, 33].argmax()] = "Macrophages"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 33].argmax(), axis = 0)
        tr = self.cluster_centers[:, 42] + self.cluster_centers[:, 7]
        self.cellTypeDict[tr.argmax()] = "IL17"
        self.celustercenters = np.delete(self.cluster_centers, tr.argmax(), axis = 0)
        self.cellTypeDict[self.cluster_centers[:, 46].argmax()] = "T cells"
        self.celustercenters = np.delete(self.cluster_centers, self.cluster_centers[:, 46].argmax(), axis = 0)

        lis = list()
        for i in range(len(self.cluster_centers)):
            if i not in self.cellTypeDict:
                lis.append(i)

        for i in lis:
            self.cellTypeDict[i] = "Others"
                            
if __name__ == "__main__":
    with open("data/som.pkl", "rb") as f:
        som = pickle.load(f)
    cluster_centers = som.cluster_centers_
    cluster_centers = cluster_centers.reshape(16, 56)
    cellType = CellType(cluster_centers)
    cellType.setCellType()
    # dfn["pred"].apply(lambda x: cellType.cellType[x])
    for i in range(15):
        print(cellType.cellTypeDict[i], i)
